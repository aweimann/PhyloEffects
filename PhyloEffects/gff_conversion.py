#Functions to work with GFF files

from io import StringIO
from Bio import SeqIO
import gffutils as gff
import pyranges
import gzip

#Clean other "##" starting lines from gff file, as it confuses parsers
#Taken from prokka.py in Panaroo
def clean_gff_string(gff_string):
    splitlines = gff_string.splitlines()
    lines_to_delete = []
    for index in range(len(splitlines)):
        if '##sequence-region' in splitlines[index]:
            lines_to_delete.append(index)
    for index in sorted(lines_to_delete, reverse=True):
        del splitlines[index]
    cleaned_gff = "\n".join(splitlines)
    return cleaned_gff

#Takes a GFF file and returns a list of lists
#Each entry in list is a gene with 4 components - gene name, gene start, gene end, strand
def convertGFF(gff_file_name, output_dir):
    if gff_file_name.endswith(".gz"):
        gff_file = gzip.open(gff_file_name, "rt", encoding='utf-8')
    else:
        gff_file = open(gff_file_name, "r", encoding='utf-8')

    #Open file, split into genes and sequence
    lines = gff_file.read().replace(",", "")
    split = lines.split("##FASTA")

    parsed_gff = gff.create_db(clean_gff_string(split[0]),
                                dbfn = ":memory:",
                                force = True,
                                keep_order = True,
                               merge_strategy="create_unique",
                               sort_attribute_values=True,
                                from_string = True)
    
    gene_annotation = {} 
    
    pyr_chr, pyr_start, pyr_stop, pyr_strand, pyr_id = [[] for i in range(5)]
    attributes_header = set()
    for entry in parsed_gff.all_features(featuretype = ()):
        if entry.featuretype == "exon":
            continue
        if entry.featuretype != "gene" and entry.featuretype != "pseudogene" and entry.featuretype != "region":
            attribute_names = []
            if 'locus_tag' in entry.attributes:
                locus_tag = entry.attributes['locus_tag'][0]
            elif 'locus' in entry.attributes:
                locus_tag = entry.attributes['locus'][0]
            else:
                continue
            if 'gene' in entry.attributes:
                attribute_names.append('gene')
            if 'old_locus_tag' in entry.attributes :
                attribute_names.append('old_locus_tag')
            if 'product' in entry.attributes :
                attribute_names.append('product')
            attributes = {}
            for name in attribute_names:
                attributes[name] = entry.attributes[name][0]
            for name in attribute_names:
                attributes_header.add(name)
            # Legionella pneumophila ST1
            # if entry.seqid == "NC_006368.1":
            #     entry.start += 131885
            #     entry.stop += 131885
            gene_annotation[entry.id] = [entry.start, entry.stop, entry.strand, locus_tag,
                                         entry.featuretype, attributes]
            #add stop codon to coordinates? (entry.stop + 1)
            pyr_chr.append(entry.seqid)
            pyr_start.append(entry.start)
            pyr_stop.append(entry.stop)
            pyr_id.append(entry.id)

    gene_ranges = pyranges.from_dict({"Chromosome": pyr_chr, "Start": pyr_start, "End": pyr_stop, "Id": pyr_id})
    write_gene_annotation(gene_annotation, attributes_header, output_dir)
    return gene_annotation, gene_ranges


def write_gene_annotation(gene_coordinates, attribute_header, output_dir):
    with open(output_dir + "gene_annotation.txt", 'w') as f:
        # header
        f.write("start\tend\tstrand\tlocus_tag\tfeature")
        # add attribute header
        for attr in attribute_header:
            f.write("\t" + attr)
        f.write("\n")
        for value_list in gene_coordinates.values():
            f.write("\t".join([str(i) for i in value_list[0:5]]))
            # add attributes
            for attr in attribute_header:
                if attr in value_list[5]:
                    f.write("\t" + value_list[5][attr])
                else:
                    f.write("\t")
            f.write("\n")


def extract_intergenic_regions(gene_annotation, gene_ranges):
    """
    Extract intergenic regions from parsed GFF data.
    Returns a dictionary with intergenic region info: upstream_gene, downstream_gene, start, end.
    """
    intergenic_regions = {}

    # Group genes by chromosome
    genes_by_chr = {}
    for gene_id, (start, stop, strand, locus_tag, featuretype, attributes) in gene_annotation.items():
        chr_info = None
        for feature in gene_ranges.features:
            if feature.id == gene_id:
                chr_info = feature.chromosome
                break
        if chr_info not in genes_by_chr:
            genes_by_chr[chr_info] = []
        genes_by_chr[chr_info].append((start, stop, gene_id, locus_tag))

    # Sort genes by start position within each chromosome
    for chr_id in genes_by_chr:
        genes_by_chr[chr_id].sort(key=lambda x: x[0])

    # Extract intergenic regions
    region_id = 0
    for chr_id, sorted_genes in genes_by_chr.items():
        for i in range(len(sorted_genes) - 1):
            upstream_start, upstream_stop, upstream_id, upstream_locus = sorted_genes[i]
            downstream_start, downstream_stop, downstream_id, downstream_locus = sorted_genes[i + 1]

            # Define intergenic region: from end of upstream gene to start of downstream gene
            if upstream_stop < downstream_start:
                region_start = upstream_stop + 1
                region_end = downstream_start - 1

                intergenic_regions[f"intergenic_{region_id}"] = {
                    "chromosome": chr_id,
                    "start": region_start,
                    "end": region_end,
                    "length": region_end - region_start + 1,
                    "upstream_gene_id": upstream_id,
                    "upstream_locus_tag": upstream_locus,
                    "downstream_gene_id": downstream_id,
                    "downstream_locus_tag": downstream_locus
                }
                region_id += 1

    return intergenic_regions


def write_intergenic_regions(intergenic_regions, output_dir):
    with open(output_dir + "intergenic_regions.txt", 'w') as f:
        f.write("chromosome\tstart\tend\tlength\tupstream_gene_id\tupstream_locus_tag\tdownstream_gene_id\tdownstream_locus_tag\n")
        for region_id in sorted(intergenic_regions.keys(), key=lambda x: int(x.split('_')[1])):
            region = intergenic_regions[region_id]
            f.write("\t".join([
                str(region["chromosome"]),
                str(region["start"]),
                str(region["end"]),
                str(region["length"]),
                str(region["upstream_gene_id"]),
                str(region["upstream_locus_tag"]),
                str(region["downstream_gene_id"]),
                str(region["downstream_locus_tag"])
            ]) + "\n")

