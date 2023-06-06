#Functions to work with GFF files

from io import StringIO
from Bio import SeqIO
import gffutils as gff
import pyranges

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
def convertGFF(gff_file_name):
    gff_file = open(gff_file_name, "r")

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
    for entry in parsed_gff.all_features(featuretype = ()):
        if entry.featuretype != "gene" and entry.featuretype != "pseudogene" and entry.featuretype != "region":
            if 'locus_tag' in entry.attributes :
                gene_annotation[entry.id] = [entry.start, entry.stop, entry.strand, entry.attributes['locus_tag'][0],
                                             entry.featuretype]
            elif 'locus' in entry.attributes :
                gene_annotation[entry.id] = [entry.start, entry.stop, entry.strand, entry.attributes['locus'][0],
                                             entry.featuretype]
            else:
                continue
            #add stop codon to coordinates? (entry.stop + 1)
            pyr_chr.append(entry.seqid)
            pyr_start.append(entry.start)
            pyr_stop.append(entry.stop)
            pyr_id.append(entry.id)

    gene_ranges = pyranges.from_dict({"Chromosome": pyr_chr, "Start": pyr_start, "End": pyr_stop, "Id": pyr_id})

    return(gene_annotation, gene_ranges)
