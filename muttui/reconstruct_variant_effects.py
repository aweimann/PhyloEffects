# Takes the alignment and tree output by treetime and reconstructs the mutational spectrum

from Bio import SeqUtils
from Bio.Seq import Seq
import numpy as np
import array
import pyranges
from Bio import AlignIO
import pandas as pd
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from variant_effect import VariantEffect
import time

# write out effect TODO determine proximity to coding features nearby

translation_table = np.array([[[b'K', b'N', b'K', b'N', b'X'],
                               [b'T', b'T', b'T', b'T', b'T'],
                               [b'R', b'S', b'R', b'S', b'X'],
                               [b'I', b'I', b'M', b'I', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'Q', b'H', b'Q', b'H', b'X'],
                               [b'P', b'P', b'P', b'P', b'P'],
                               [b'R', b'R', b'R', b'R', b'R'],
                               [b'L', b'L', b'L', b'L', b'L'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'E', b'D', b'E', b'D', b'X'],
                               [b'A', b'A', b'A', b'A', b'A'],
                               [b'G', b'G', b'G', b'G', b'G'],
                               [b'V', b'V', b'V', b'V', b'V'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'*', b'Y', b'*', b'Y', b'X'],
                               [b'S', b'S', b'S', b'S', b'S'],
                               [b'*', b'C', b'W', b'C', b'X'],
                               [b'L', b'F', b'L', b'F', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']]])


# Converts the positional translation to a dictionary with alignment positions as keys and genome positions as values
def convertTranslation(positionsFile):
    # Import the translation file
    positions = open(positionsFile.name).readlines()

    conversion = {}

    for eachPosition in positions:
        conversion[int(eachPosition.strip().split("\t")[0])] = int(eachPosition.strip().split("\t")[1])

    return conversion


# Creates a positional translation dictionary with alignment sites as keys and values
# Used when --all_sites is specified so there is no need to convert the sites
def all_sites_translation(alignment):
    conversion = {}

    for eachPosition in range(len(alignment[0].seq)):
        conversion[int(eachPosition) + 1] = int(eachPosition) + 1

    return conversion


# Extracts the mutations in branch_mutations.txt to a dictionary with branch names as keys and mutations as values
def get_branch_mutation_dict(branchFile, translation):
    branchDict = {}

    # Iterate through the mutations and add to branchDict
    with open(branchFile) as fileobject:
        # Skip header
        next(fileobject)

        for line in fileobject:
            node = line.strip().split("\t")[0]

            if node in branchDict:
                branchDict[node].append([line.strip().split("\t")[1], int(line.strip().split("\t")[2]),
                                         translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]])
            else:
                branchDict[node] = [[line.strip().split("\t")[1], int(line.strip().split("\t")[2]),
                                     translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]]]

    return branchDict


def get_mutations_from_alignment(alignment, ref, translation):
    clade2mutations = {}
    clades = []
    for record in SeqIO.parse(alignment, "fasta"):
        variants = []
        name = record.name
        clades.append(name)
        clade2mutations[name] = []
        for i, pos in zip(record.seq, range(1, len(record.seq) + 1)):
            j = ref[translation[pos] - 1]
            if i != "N" and i != "-" and j != "N" and i.upper() != j.upper():
                clade2mutations[name].append((j, pos, translation[pos], i))
    return clade2mutations, clades


# Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the alignment
# are assumed
# and the root sequence from the ancestral reconstruction is used
def get_reference(reference, all_sites, alignment, positionTranslation):
    if all_sites:
        # Get the reference as the root of the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                reference_sequence = sequence.seq

    # Use the provided reference as the reference sequence and substitute in the nucleotides in the inferred root
    # sequence at each variable position. This means the updateReference function uses the root sequence as a
    # starting sequence and updates based on this
    else:
        ref = AlignIO.read(reference.name, "fasta")
        # Extract the sequence of the reference and ensure it is uppercase
        refSeq = ref[0].seq.upper()

        # Reverse the position translation so genome positions are keys and alignment positions are values
        # Can then iterate through the genome positions and check if they are in this
        toChange = {j: i for i, j in positionTranslation.items()}

        # Extract the root sequences from the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                root_seq = sequence.seq

        # Will be the original reference with the root sequence base substituted at each reconstructed position
        reference_sequence = ""

        # Iterate through the length of the reference genome. If the position is in toChange, take the dictionary
        # value position from
        # the root sequence which only contains positions in the translate file. Otherwise, take the position from
        # the reference
        for pos in range(len(refSeq)):
            # Check if the position is variable and needs to be taken from the root
            if (pos + 1) in toChange:
                reference_sequence += root_seq[toChange[pos + 1] - 1]
            else:
                reference_sequence += refSeq[pos]

    return reference_sequence


# Updates the original reference sequence so it includes the mutations acquired up to the given clade
# Converts the reference to an array, iterates through the upstream branches and updates each of the positions
# that mutated along an upstream branch to the mutated base
# If a position has changed multiple times along the upstream branches, this keeps the most recent change as the
# clades are iterated through from the root through to the most recent upstream branch
def updateReference(tree, clade, branchMutationDict, refSeq):
    # Convert the reference sequence to an array
    referenceArray = array.array("u", refSeq)

    # Iterate through the upstream branches leading to the node at the start of the current branch
    for upstreamClade in tree.get_path(clade)[:-1]:
        # Check if there are any mutations along the current upstream branch
        if upstreamClade.name in branchMutationDict:
            # Iterate through the previous mutations and update the reference sequence
            for eachMutation in branchMutationDict[upstreamClade.name]:
                referenceArray[eachMutation[2] - 1] = eachMutation[3]

    return "".join(referenceArray)


# Removes mutations that are at the start or end of the genome or do not involve 2 nucleotides
# Writes the removed mutations to outMutationsNotUsed
# Splits double substitutions into a separate list
# Returns the filtered mutations as a list
def filter_mutations(branch_mutations, clade, nucleotides, reference_length, out_mutations_not_used):
    positions_to_remove = []
    double_substitutions = []

    # Check if there are double substitutions along the branch and remove these mutations
    for mutation1 in range(0, len(branch_mutations)):
        # Check if the following mutation is at the adjacent genome position
        if (mutation1 + 1) != len(branch_mutations):
            if branch_mutations[mutation1 + 1][2] == (branch_mutations[mutation1][2] + 1):
                positions_to_remove.append(mutation1)
                positions_to_remove.append(mutation1 + 1)
                double_substitutions.append(branch_mutations[mutation1])
                double_substitutions.append(branch_mutations[mutation1 + 1])
        # Check if the position is at the start or end of the genome
        if (branch_mutations[mutation1][2] == 1) or (branch_mutations[mutation1][2] == reference_length):
            positions_to_remove.append(mutation1)
            # Write the mutation to the mutations not used file
            out_mutations_not_used.write(
                branch_mutations[mutation1][0] + str(branch_mutations[mutation1][1]) + branch_mutations[mutation1][
                    3] + "," + branch_mutations[mutation1][0] + str(branch_mutations[mutation1][2]) +
                branch_mutations[mutation1][3] + "," + clade.name + ",End_of_genome\n")
        # Check if the mutation does not involve 2 nucleotides
        elif (branch_mutations[mutation1][0] not in nucleotides) or (branch_mutations[mutation1][3] not in nucleotides):
            positions_to_remove.append(mutation1)
            out_mutations_not_used.write(
                branch_mutations[mutation1][0] + str(branch_mutations[mutation1][1]) + branch_mutations[mutation1][
                    3] + "," + branch_mutations[mutation1][0] + str(branch_mutations[mutation1][2]) +
                branch_mutations[mutation1][3] + "," + clade.name + "Mutation_does_not_involve_two_nucleotides\n")

    # If there is a tract of 3 or more substitutions at adjacent positions, some of the mutations
    # will be in doubleSubstitutions twice. Extract the unique double substitution positions
    uniqueDoubleSubstitutions = []
    if len(double_substitutions) != 0:
        dsSet = set()
        for ds in double_substitutions:
            if tuple(ds) not in dsSet:
                uniqueDoubleSubstitutions.append(ds)
                dsSet.add(tuple(ds))
    ####Write the double substitutions
    ####Will be removed once double substitutions are incorporated
    for uds in uniqueDoubleSubstitutions:
        out_mutations_not_used.write(uds[0] + str(uds[1]) + uds[3] + "," + uds[0] + str(uds[2]) + uds[
            3] + "," + clade.name + ",Double_substitution\n")

    # Remove the positions that will not be included in the spectrum
    double_substitution_dict = {}
    if len(positions_to_remove) != 0:
        # Identify the unique set of positions, if there are 3 or more consecutive positions to be removed,
        # at least one of these positions will be in positionsToRemove more than once
        uniquePositionsToRemove = list(set(positions_to_remove))
        # Flag double substitutions
        for ele in sorted(uniquePositionsToRemove, reverse=True):
            double_substitution_dict[ele] = True
    return (double_substitution_dict, uniqueDoubleSubstitutions)


# Translates a given nucleotide sequence to protein
def translate_sequence(sequence, strand, check_cds=False):
    if strand == "+":
        return Seq("".join(sequence)).translate(cds=check_cds, table=11)
    else:
        return Seq("".join(sequence)).reverse_complement().translate(cds=check_cds, table=11)


# Identifies the amino acid position of a nucleotide within a gene
# If the gene is on the positive strand, this is the nucleotide position divided by 3, rounded up
# If the gene is on the negative strand, this is the nucleotide position from the end of the gene divided by 3,
# rounded up
def extract_position(gene_coordinates, position_in_gene):
    if gene_coordinates[2] == "+":
        return (position_in_gene + 1) // 3 + ((position_in_gene + 1) % 3 != 0)
    else:
        reverse_position = gene_coordinates[1] - (gene_coordinates[0] + position_in_gene)
        return (reverse_position + 1) // 3 + ((reverse_position + 1) % 3 != 0)


# Extracts synonymous mutations along a given branch
# Used when --synonymous is specified
def extract_synonymous(clade, branch_mutations, updated_reference, reference_sequence,
                       variant_effect2clades, gene_coordinates, position_gene,
                       output_dir):
    # Gene sequences at the upstream node
    upstream_genes = dict()
    # Gene sequences containing mutations along the branch
    downstream_genes = dict()
    # Reference genes
    reference_genes = dict()

    # Protein sequences
    upstream_aas = dict()
    downstream_aas = dict()
    reference_aas = dict()

    # reverse complement
    rc_dss = {}
    rc_uss = {}
    rc_rss = {}

    positions_to_remove = []

    # open file for writing effect predictions to disk
    effects = open(output_dir + "variant_effect_predictions.txt", "a")

    # intersect genes with variants
    var_chromosome, var_start, var_end, var_ref, var_alt = [[] for i in range(5)]
    chromosome = position_gene.chromosomes[0]
    for i in branch_mutations:
        var_chromosome.append(chromosome)
        var_start.append(i[2])
        var_end.append(i[2] + 1)
        var_ref.append(i[0])
        var_alt.append(i[3])
    var_ranges = pyranges.from_dict(
        {"Chromosome": var_chromosome, "Start": var_start, "End": var_end, "ref": var_ref, "alt": var_alt})
    gene_overlaps = var_ranges.nearest(position_gene).as_df()
    gene_overlaps = gene_overlaps.query('Distance == 0')
    knearest = var_ranges.overlap(position_gene, invert=True).k_nearest(position_gene, k=2).as_df()
    genes_proximity = pd.DataFrame()
    node = str(clade)
    if knearest.shape[0] != 0:
        genes_us = knearest.query('Distance < 0').assign(
            left=lambda x: x.groupby('Start')['Distance'].transform('max')).query('Distance == left')
        genes_ds = knearest.query('Distance > 0').assign(
            right=lambda x: x.groupby('Start')['Distance'].transform('min')).query('Distance == right')
        genes_proximity = pd.concat([genes_us, genes_ds])
    aa_pos2effect = {}
    print("intergenic variants")
    for row_tuple in genes_proximity.itertuples():
        gene_name = row_tuple.Id
        pos = str(row_tuple.Start)
        variant_effect = VariantEffect(pos, row_tuple.ref, row_tuple.alt, node)
        variant_effect.impact = "MODIFIER"
        distance = row_tuple.Distance
        if gene_coordinates[gene_name][2] == "+":
            if distance > 0:
                mutation_type = "upstream_gene_variant"
            else:
                mutation_type = "downstream_gene_variant"
        else:
            if distance > 0:
                mutation_type = "downstream_gene_variant"
            else:
                mutation_type = "upstream_gene_variant"
        variant_effect.mutation_type = mutation_type

        variant_effect.locus_tag = gene_coordinates[gene_name][3]
        variant_effect.upstream_allele = row_tuple.ref
        variant_effect.downstream_allele = row_tuple.alt
        aa_pos2effect["%s_%s" % (variant_effect.position, variant_effect.locus_tag)] = variant_effect
        if variant_effect not in variant_effect2clades:
            variant_effect2clades[variant_effect] = [node]
        else:
            variant_effect2clades[variant_effect].append(node)
    print("first pass")
    for row in gene_overlaps.itertuples():
        # for mutation in branchMutations:
        mutation = [row.ref, 0, row.Start, row.alt]
        gene_name = row.Id
        # Iterate through the genes the position is in, check if they are already present, mutate if so, if not
        # add to dictionaries and mutate
        # Position of the mutation in the gene, zero based
        position_in_gene = mutation[2] - gene_coordinates[gene_name][0]
        if gene_name in downstream_genes:
            downstream_genes[gene_name][position_in_gene] = mutation[3]
        else:
            upstream_genes[gene_name] = array.array("u", updated_reference[
                                                        (gene_coordinates[gene_name][0] - 1):
                                                        gene_coordinates[gene_name][
                                                            1]])
            reference_genes[gene_name] = array.array("u", reference_sequence[(gene_coordinates[gene_name][0] - 1):
                                                                           gene_coordinates[gene_name][1]])
            downstream_genes[gene_name] = array.array("u", updated_reference[(gene_coordinates[gene_name][0] - 1):
                                                                           gene_coordinates[gene_name][
                                                                               1]])  # Mutate the position in the
            # downstream gene downstream_genes[geneName][positionInGene] = mutation[3]
            downstream_genes[gene_name][position_in_gene] = mutation[3]

    # Iterate through the mutations, check if the translated position of the mutation changes along the branch
    # Second iteration needed because 1st and 3rd codon positions may both change along a branch
    # dictionary aa position vs change
    is_pseudogene = {}
    start = time.time()
    print('second pass')
    for row in gene_overlaps.itertuples():
        mutation = [row.ref, 0, row.Start, row.alt]
        gene_name = row.Id
        # initialise fields for effect prediction
        # upstream_aa, downstream_aa, reference_aa, upstream_codon, downstream_codon, reference_codon, impact, \
        pos = str(mutation[2])
        upstream_allele = mutation[0]
        downstream_allele = mutation[3]
        # Check if the mutation is in a CDS, if the position is in an RNA gene assign moderate effect, else annotate
        # by nearby genes
        # Will change to False if the amino acid of the mutation has changed in any gene it is present in
        synonymous = True
        # Iterate through the genes the position is in, check if the position's amino acid changes, if so change
        # synonymous to False
        # for geneName in positionGene[mutation[2]].split("____"):
        feature_type = gene_coordinates[gene_name][4]
        variant_effect = VariantEffect(pos, upstream_allele, downstream_allele, node)
        variant_effect.upstream_allele = upstream_allele
        variant_effect.downstream_allele = downstream_allele
        variant_effect.locus_tag = gene_coordinates[gene_name][3]
        if feature_type != "CDS":
            variant_effect.impact = "MODERATE"
            variant_effect.mutation_type = "non-coding_gene"
            synonymous = False
            aa_pos2effect[pos] = variant_effect
            if variant_effect not in variant_effect2clades:
                variant_effect2clades[variant_effect] = [node]
            else:
                variant_effect2clades[variant_effect].append(node)
            continue
        # Position of the mutation in the gene, zero based
        # determine if this is a proper CDS otherwise set pseudogene flag
        # TODO selenocysteine stop codon
        # translate without the check cds flag
        if gene_name not in upstream_aas:
            try:
                translate_sequence(upstream_genes[gene_name], gene_coordinates[gene_name][2], check_cds=True)
            except Exception as err:
                is_pseudogene[gene_name] = True
                continue
            upstream_aas[gene_name] = translate_sequence(upstream_genes[gene_name], gene_coordinates[gene_name][2])
            reference_aas[gene_name] = translate_sequence(reference_genes[gene_name], gene_coordinates[gene_name][2])
            downstream_aas[gene_name] = translate_sequence(downstream_genes[gene_name], gene_coordinates[gene_name][2])
            # transcribe gene
            if gene_coordinates[gene_name][2] == "+":
                rc_dss[gene_name] = Seq("".join(downstream_genes[gene_name]))
                rc_uss[gene_name] = Seq("".join(upstream_genes[gene_name]))
                rc_rss[gene_name] = Seq("".join(reference_genes[gene_name]))
            else:
                # reverse complement
                rc_dss[gene_name] = Seq("".join(downstream_genes[gene_name])).reverse_complement()
                rc_uss[gene_name] = Seq("".join(upstream_genes[gene_name])).reverse_complement()
                rc_rss[gene_name] = Seq("".join(reference_genes[gene_name])).reverse_complement()
        if gene_name in is_pseudogene:
            variant_effect.pseudogene = "pseudogene"
        # position in gene depending on what strand the cds is on
        position_in_gene = mutation[2] - gene_coordinates[gene_name][0]
        aa_position = extract_position(gene_coordinates[gene_name], position_in_gene) - 1
        # extract codon
        rc_ds = rc_dss[gene_name]
        rc_us = rc_uss[gene_name]
        rc_rs = rc_rss[gene_name]
        variant_effect.reference_codon = rc_rs[(aa_position * 3): (aa_position * 3 + 3)]
        variant_effect.upstream_codon = rc_us[(aa_position * 3): (aa_position * 3 + 3)]
        variant_effect.downstream_codon = rc_ds[(aa_position * 3): (aa_position * 3 + 3)]
        # amino acid at substitution position
        variant_effect.reference_aa = reference_aas[gene_name][aa_position]
        variant_effect.downstream_aa = downstream_aas[gene_name][aa_position]
        variant_effect.upstream_aa = upstream_aas[gene_name][aa_position]
        # check for start codon substitution, as we're using check_cds = False, non ATG start codons are
        # translated into non methionine amino acids
        if aa_position == 0:
            if variant_effect.upstream_codon == "GTG" or variant_effect.upstream_codon == "TTG":
                variant_effect.upstream_aa = "M"
            if variant_effect.downstream_codon == "GTG" or variant_effect.downstream_codon == "TTG":
                variant_effect.downstream_aa = "M"
            if variant_effect.reference_codon == "GTG" or variant_effect.reference_codon == "TTG":
                variant_effect.reference_aa = "M"
        # check if amino acid changed due to the current substitution
        if variant_effect.upstream_aa != variant_effect.downstream_aa:
            # stop codon -> different codon
            if variant_effect.upstream_aa == "*":
                variant_effect.mutation_type = "stop_lost"
                downstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[variant_effect.downstream_aa]
                upstream_aa3 = "*"
                variant_effect.upstream_aa = "*"
                variant_effect.impact = "HIGH"
            # stop codon inserted
            elif variant_effect.downstream_aa == "*":
                variant_effect.mutation_type = "stop_gained"
                upstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[variant_effect.upstream_aa]
                downstream_aa3 = "*"
                variant_effect.downstream_aa = "*"
                variant_effect.impact = "HIGH"
            else:
                # start codon -> different codon
                if aa_position == 0 and variant_effect.upstream_aa == "M" and variant_effect.downstream_aa != "M":
                    variant_effect.impact = "HIGH"
                    variant_effect.mutation_type = "start_lost"
                # regular non-synonymous SNP
                else:
                    variant_effect.impact = "MODERATE"
                    variant_effect.mutation_type = "missense_variant"
                downstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[variant_effect.downstream_aa]
                upstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[variant_effect.upstream_aa]
                aa_change = "p.%s%s%s" % (upstream_aa3, aa_position + 1, downstream_aa3)
                variant_effect.aa_change = aa_change
            synonymous = False
        # synonymous SNP
        else:
            variant_effect.impact = "LOW"
            variant_effect.mutation_type = "synonymous_variant"
        # account for codons mutated at different positions on the same branch
        dict_key = "%s_%s_%s" % (node, variant_effect.locus_tag, aa_position + 1)
        if dict_key in aa_pos2effect:
            variant_effect.multi_codon_substitution = "multi_codon_substitution"
            # correct position to beginning of the codon
            if gene_coordinates[gene_name][2] == "+":
                pos = str(int(pos) - (position_in_gene % 3))
                variant_effect.upstream_allele = variant_effect.upstream_codon
                variant_effect.downstream_allele = variant_effect.downstream_codon
            else:
                pos = str(gene_coordinates[gene_name][1] - aa_position * 3 - 2)
                variant_effect.upstream_allele = variant_effect.upstream_codon.reverse_complement()
                variant_effect.downstream_allele = variant_effect.downstream_codon.reverse_complement()
        aa_pos2effect[dict_key] = variant_effect
        if variant_effect not in variant_effect2clades:
            variant_effect2clades[variant_effect] = [node]
        else:
            variant_effect2clades[variant_effect].append(node)
        # If the mutation is nonsynonymous within any gene, remove it
        if synonymous == False:
            positions_to_remove.append(pos)
    # write effects to effect prediction out file
    # for variant_effect in aa_pos2effect.values():
    #    effects.write(variant_effect.to_string() + "\n")

    # Flag the positions that will not be included in the spectrum
    synonymous_substitution_dict = {}
    if len(positions_to_remove) != 0:
        # Identify the unique set of positions
        uniquePositionsToRemove = list(set(positions_to_remove))
        # Flag synonymous substitutions
        for ele in sorted(uniquePositionsToRemove, reverse=True):
            synonymous_substitution_dict[ele] = False

    return synonymous_substitution_dict
