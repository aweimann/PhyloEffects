#Takes the alignment and tree output by treetime and reconstructs the mutational spectrum

import argparse
from Bio import SeqUtils
from Bio.Seq import Seq
import numpy as np
import array
import pyranges
from Bio import Phylo
import pandas as pd

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


#Converts the positional translation to a dictionary with alignment positions as keys and genome positions as values
def convertTranslation(positionsFile):
    #Import the translation file
    positions = open(positionsFile.name).readlines()

    conversion = {}

    for eachPosition in positions:
        conversion[int(eachPosition.strip().split("\t")[0])] = int(eachPosition.strip().split("\t")[1])
    
    return(conversion)

#Creates a positional translation dictionary with alignment sites as keys and values
#Used when --all_sites is specified so there is no need to convert the sites
def allSitesTranslation(alignment):
    conversion = {}

    for eachPosition in range(len(alignment[0].seq)):
        conversion[int(eachPosition) + 1] = int(eachPosition) + 1
    
    return(conversion)

#Extracts the mutations and positions along a given branch
def getMutations(branch, translation):
    mutations = []

    for mutation in branch.split('="')[1].split('"')[0].split(","):
        mutations.append([mutation[0], int(mutation[1:-1]), translation[int(mutation[1:-1])], mutation[-1]])
    
    return(mutations)

#Gets the name of the upstream clade
def getParentName(tree, clade):
    #Check if the clade is one downstream of the root
    if len(tree.get_path(clade)) >= 2:
        #Nodes with no mutations along the upstream branch have their IDs imported as name, nodes with mutations as confidence
        if tree.get_path(clade)[-2].confidence:
            node = tree.get_path(clade)[-2].confidence
        else:
            node = tree.get_path(clade)[-2].name
    else:
        node = str(tree.root)
    return(node)

#Identifies the name of the branch
def getBranchName(tree, clade):
    parentNode = getParentName(tree, clade)

    if clade.is_terminal():
        return(parentNode + "->" + clade.name)
    elif clade.confidence:
        return(parentNode + "->" + clade.confidence)
    else:
        return(parentNode + "->" + clade.name)

#Extracts the mutations to a dictionary with branch names as keys and mutations as values
def getBranchDict(tree, positionTranslation):
    branchDict = {}

    for clade in tree.find_clades():
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:
            #Check if there are mutations along the branch
            if (clade.comment) and (clade.comment != '[&mutations=""]'):
                branchDict[getBranchName(tree, clade)] = getMutations(clade.comment, positionTranslation)
            else:
                branchDict[getBranchName(tree, clade)] = "None"
    
    return(branchDict)

#Extracts the mutations in branch_mutations.txt to a dictionary with branch names as keys and mutations as values
def getBranchMutationDict(branchFile, translation):
    branchDict = {}

    #Iterate through the mutations and add to branchDict
    with open(branchFile) as fileobject:
        #Skip header
        next(fileobject)

        for line in fileobject:
            node = line.strip().split("\t")[0]

            if node in branchDict:
                branchDict[node].append([line.strip().split("\t")[1], int(line.strip().split("\t")[2]), translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]])
            else:
                branchDict[node] = [[line.strip().split("\t")[1], int(line.strip().split("\t")[2]), translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]]]
    
    return(branchDict)


#Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the alignment are assumed
#and the root sequence from the ancestral reconstruction is used
def getReference(reference, all_sites, alignment, positionTranslation):
    if all_sites:
        #Get the reference as the root of the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                referenceSequence = sequence.seq
    
    #Use the provided reference as the reference sequence and substitute in the nucleotides in the inferred root
    #sequence at each variable position. This means the updateReference function uses the root sequence as a
    #starting sequence and updates based on this
    else:
        ref = Phylo.AlignIO.read(reference.name, "fasta")
        #Extract the sequence of the reference and ensure it is uppercase
        refSeq = ref[0].seq.upper()

        #Reverse the position translation so genome positions are keys and alignment positions are values
        #Can then iterate through the genome positions and check if they are in this
        toChange = {j:i for i, j in positionTranslation.items()}

        #Extract the root sequences from the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                root_seq = sequence.seq
        
        #Will be the original reference with the root sequence base substituted at each reconstructed position
        referenceSequence = ""

        #Iterate through the length of the reference genome. If the position is in toChange, take the dictionary value position from
        #the root sequence which only contains positions in the translate file. Otherwise, take the position from the reference
        for pos in range(len(refSeq)):
            #Check if the position is variable and needs to be taken from the root
            if (pos + 1) in toChange:
                referenceSequence += root_seq[toChange[pos + 1] - 1]
            else:
                referenceSequence += refSeq[pos]
    
    return(referenceSequence)


#Updates the original reference sequence so it includes the mutations acquired up to the given clade
#Converts the reference to an array, iterates through the upstream branches and updates each of the positions
#that mutated along an upstream branch to the mutated base
#If a position has changed multiple times along the upstream branches, this keeps the most recent change as the
#clades are iterated through from the root through to the most recent upstream branch
def updateReference(tree, clade, branchMutationDict, refSeq):
    #Convert the reference sequence to an array
    referenceArray = array.array("u", refSeq)

    #Iterate through the upstream branches leading to the node at the start of the current branch
    for upstreamClade in tree.get_path(clade)[:-1]:
        #Check if there are any mutations along the current upstream branch
        if upstreamClade.name in branchMutationDict:
            #Iterate through the previous mutations and update the reference sequence
            for eachMutation in branchMutationDict[upstreamClade.name]:
                referenceArray[eachMutation[2] - 1] = eachMutation[3]
    
    return("".join(referenceArray))

#Identifies the label category of a given branch
#By default returns None if the branch is a transition between labels
#If --include_all_branches is used, no branches are labelled None and transition branches are
#labelled with the downstream node's label
def getBranchCategory(tree, clade, include_all_branches):
    #Identify the label of the upstream node
    #Check if the upstream node is the root
    if len(tree.get_path(clade)) >= 2:
        #Check if the label of the upstream clade is the same as the clade
        if clade.clade_label == tree.get_path(clade)[-2].clade_label:
            return(clade.clade_label)
        elif include_all_branches:
            return(clade.clade_label)
        else:
            return(None)
    
    else:
        if clade.clade_label == tree.root.clade_label:
            return(clade.clade_label)
        elif include_all_branches:
            return(clade.clade_label)
        else:
            return(None)

#Removes mutations that are at the start or end of the genome or do not involve 2 nucleotides
#Writes the removed mutations to outMutationsNotUsed
#Splits double substitutions into a separate list
#Returns the filtered mutations as a list
def filterMutations(branchMutations, clade, nucleotides, referenceLength, outMutationsNotUsed):
    positionsToRemove = []
    doubleSubstitutions = []

    #Check if there are double substitutions along the branch and remove these mutations
    for mutation1 in range(0, len(branchMutations)):
        #Check if the following mutation is at the adjacent genome position
        if (mutation1 + 1) != len(branchMutations):
            if branchMutations[mutation1 + 1][2] == (branchMutations[mutation1][2] + 1):
                positionsToRemove.append(mutation1)
                positionsToRemove.append(mutation1 + 1)
                doubleSubstitutions.append(branchMutations[mutation1])
                doubleSubstitutions.append(branchMutations[mutation1 + 1])
        #Check if the position is at the start or end of the genome
        if (branchMutations[mutation1][2] == 1) or (branchMutations[mutation1][2] == referenceLength):
            positionsToRemove.append(mutation1)
            #Write the mutation to the mutations not used file
            outMutationsNotUsed.write(branchMutations[mutation1][0] + str(branchMutations[mutation1][1]) + branchMutations[mutation1][3] + "," + branchMutations[mutation1][0] + str(branchMutations[mutation1][2]) + branchMutations[mutation1][3] + "," + clade.name + ",End_of_genome\n")
        #Check if the mutation does not involve 2 nucleotides
        elif (branchMutations[mutation1][0] not in nucleotides) or (branchMutations[mutation1][3] not in nucleotides):
            positionsToRemove.append(mutation1)
            outMutationsNotUsed.write(branchMutations[mutation1][0] + str(branchMutations[mutation1][1]) + branchMutations[mutation1][3] + "," + branchMutations[mutation1][0] + str(branchMutations[mutation1][2]) + branchMutations[mutation1][3] + "," + clade.name + "Mutation_does_not_involve_two_nucleotides\n")
    
    #If there is a tract of 3 or more substitutions at adjacent positions, some of the mutations
    #will be in doubleSubstitutions twice. Extract the unique double substitution positions
    uniqueDoubleSubstitutions = []
    if len(doubleSubstitutions) != 0:
        dsSet = set()
        for ds in doubleSubstitutions:
            if tuple(ds) not in dsSet:
                uniqueDoubleSubstitutions.append(ds)
                dsSet.add(tuple(ds))
    ####Write the double substitutions
    ####Will be removed once double substitutions are incorporated
    for uds in uniqueDoubleSubstitutions:
        outMutationsNotUsed.write(uds[0] + str(uds[1]) + uds[3] + "," + uds[0] + str(uds[2]) + uds[3] + "," + clade.name + ",Double_substitution\n")
    
    #Remove the positions that will not be included in the spectrum
    double_substitution_dict = {}
    if len(positionsToRemove) != 0:
        #Identify the unique set of positions, if there are 3 or more consecutive positions to be removed,
        #at least one of these positions will be in positionsToRemove more than once
        uniquePositionsToRemove = list(set(positionsToRemove))
        #Flag double substitutions
        for ele in sorted(uniquePositionsToRemove, reverse = True):
            double_substitution_dict[ele] = True 
    return(double_substitution_dict, uniqueDoubleSubstitutions)

#Translates a given nucleotide sequence to protein
def translateSequence(sequence, strand, check_cds = False):
    if strand == "+":
        return(Seq("".join(sequence)).translate(cds = check_cds, table = 11))
    else:
        return(Seq("".join(sequence)).reverse_complement().translate(cds = check_cds, table = 11))

#Identifies the amino acid position of a nucleotide within a gene
#If the gene is on the positive strand, this is the nucleotide position divided by 3, rounded up
#If the gene is on the negative strand, this is the nucleotide position from the end of the gene divided by 3, rounded up
def extractPosition(geneCoordinates, positionInGene):
    if geneCoordinates[2] == "+":
        return((positionInGene + 1)//3 + ((positionInGene + 1) % 3 != 0)) 
    else:
        reversePosition = geneCoordinates[1] - (geneCoordinates[0] + positionInGene)
        return((reversePosition + 1)//3 + ((reversePosition + 1) % 3 != 0))

#Extracts synonymous mutations along a given branch
#Used when --synonymous is specified
def extractSynonymous(clade, branchMutations, updatedReference, referenceSequence, geneCoordinates, positionGene, output_dir):
    #Gene sequences at the upstream node
    upstreamGenes = dict()
    #Gene sequences containing mutations along the branch
    downstreamGenes = dict()
    referenceGenes = dict()

    positionsToRemove = []


    
    #open file for writing effect predictions to disk
    effects = open(output_dir + "variant_effect_predictions.txt", "a")

    
    #intersect genes with variants
    var_chromosome, var_start, var_end, var_ref, var_alt  = [[] for i in range(5)]
    for i in branchMutations:
        #var_chromosome.append("chromosome") 
        var_chromosome.append(positionGene.chromosomes[0]) 
        var_start.append(i[2]) 
        var_end.append(i[2] + 1) 
        var_ref.append(i[0])
        var_alt.append(i[3])
    var_ranges = pyranges.from_dict({"Chromosome": var_chromosome, "Start": var_start, "End": var_end, "ref": var_ref, "alt": var_alt})
    gene_overlaps = var_ranges.nearest(positionGene).as_df()
    gene_overlaps = gene_overlaps.query('Distance == 0')
    knearest = var_ranges.overlap(positionGene, invert = True).k_nearest(positionGene, k = 2).as_df()
    genes_proximity = pd.DataFrame()
    node = str(clade)
    if knearest.shape[0] != 0:
        genes_us = knearest.query('Distance < 0').assign(left = lambda x: x.groupby('Start').transform('max')['Distance']).query('Distance == left')
        genes_ds = knearest.query('Distance > 0').assign(right = lambda x: x.groupby('Start').transform('min')['Distance']).query('Distance == right')
        genes_proximity = pd.concat([genes_us, genes_ds])
    aa_pos2effect = {}
    for i in genes_proximity.index:
        upstream_aa, downstream_aa, reference_aa, upstream_codon, downstream_codon, reference_codon, impact, aa_change, multi_codon_substitution, locus_tag, pseudogene = ["" for i in range(11)]
        impact = "MODIFIER"
        geneName = genes_proximity.loc[i, "Id"]
        pos = str(genes_proximity.loc[i, "Start"])
        upstream_allele = genes_proximity.loc[i, "ref"]
        downstream_allele = genes_proximity.loc[i, "alt"]
        distance = genes_proximity.loc[i, "Distance"]
        if geneCoordinates[geneName][2]  == "+":
            if distance > 0:
                mutation_type = "upstream_gene_variant"
            else:
                mutation_type = "downstream_gene_variant"
        else:
            if distance > 0:
                mutation_type = "downstream_gene_variant"
            else:
                mutation_type = "upstream_gene_variant"

        locus_tag = geneCoordinates[geneName][3]
        #write out effect TODO determine proximity to coding features nearby
        aa_pos2effect["%s_%s" % (pos, locus_tag)]= [node, pos, upstream_allele, downstream_allele, mutation_type, upstream_aa, downstream_aa, reference_aa, str(downstream_codon), str(upstream_codon), str(reference_codon), impact, aa_change, multi_codon_substitution, locus_tag, pseudogene]
    #1) group by position to get all genes affected by this position
    #Iterate through the mutations to create dictionaries of each gene that mutates along the branch
    #with their nucleotide sequences pre and post mutations
    for pos_group in gene_overlaps.groupby("Start"):
    #for mutation in branchMutations:
        #Check if the mutation is in a gene, if the position isn't in positionGene, the position is intergenic
        mutation = [pos_group[1].iloc[0].loc["ref"], 0, pos_group[0], pos_group[1].iloc[0].loc["alt"]]
        #if mutation[2] not in positionGene:
        #    continue
        #Iterate through the genes the position is in, check if they are already present, mutate if so, if not
        #add to dictionaries and mutate
        #for geneName in positionGene[mutation[2]].split("____"):
        for match in pos_group[1].index:
            geneName = pos_group[1].loc[match, "Id"]
            #Position of the mutation in the gene, zero based
            positionInGene = mutation[2] - geneCoordinates[geneName][0]
            if geneName in downstreamGenes:
                downstreamGenes[geneName][positionInGene] = mutation[3]
            else:
                upstreamGenes[geneName] = array.array("u", updatedReference[(geneCoordinates[geneName][0] - 1):geneCoordinates[geneName][1]])
                referenceGenes[geneName] = array.array("u", referenceSequence[(geneCoordinates[geneName][0] - 1):geneCoordinates[geneName][1]])
                downstreamGenes[geneName] = array.array("u", updatedReference[(geneCoordinates[geneName][0] - 1):geneCoordinates[geneName][1]]) #Mutate the position in the downstream gene downstreamGenes[geneName][positionInGene] = mutation[3]
                downstreamGenes[geneName][positionInGene] = mutation[3]
    
    #Iterate through the mutations , check if the translated position of the mutation changes along the branch
    #Second iteration needed because 1st and 3rd codon positions may both change along a branch
    #dictionary aa position vs change
    for pos_group in gene_overlaps.groupby("Start"):
    #for i, mutation in enumerate(branchMutations):
        mutation = [pos_group[1].iloc[0].loc["ref"], 0, pos_group[0], pos_group[1].iloc[0].loc["alt"]]
        #initialise fields for effect prediction
        upstream_aa, downstream_aa, reference_aa, upstream_codon, downstream_codon, reference_codon, impact, aa_change, multi_codon_substitution, locus_tag, pseudogene = ["" for i in range(11)]
        pos = str(mutation[2])
        upstream_allele = mutation[0]
        downstream_allele = mutation[3]
        #Check if the mutation is in a CDS, if the position is in an RNA gene assign moderate effect, else annotate by nearby genes 
        #Will change to False if the amino acid of the mutation has changed in any gene it is present in
        synonymous = True
        #Iterate through the genes the position is in, check if the position's amino acid changes, if so change synonymous to False
        #for geneName in positionGene[mutation[2]].split("____"):
        for match in pos_group[1].index:
            geneName = pos_group[1].loc[match, "Id"]
            feature_type = geneCoordinates[geneName][4]
            if feature_type != "CDS":
                for match in pos_group[1].index:
                    impact = "MODERATE"
                    mutation_type = "non-coding_gene"
                    geneName = pos_group[1].loc[match, "Id"]
                    locus_tag = geneCoordinates[geneName][3]
                    synonymous = False 
                    #write out effect TODO determine proximity to coding features nearby
                    aa_pos2effect[pos]= [node, pos, upstream_allele, downstream_allele, mutation_type, upstream_aa, downstream_aa, reference_aa, str(downstream_codon), str(upstream_codon), str(reference_codon), impact, aa_change, multi_codon_substitution, locus_tag, pseudogene]
                continue
            #Position of the mutation in the gene, zero based
            locus_tag = geneCoordinates[geneName][3]
            aa_change = ""
            #determine if this is a proper CDS otherwise set pseudogene flag 
            #TODO selenocysteine stop codon 
            try:
                translateSequence(upstreamGenes[geneName], geneCoordinates[geneName][2], check_cds = True)
            except Exception as err:
                pseudogene = "pseudogene"
            #translate without the check cds flag
            upstreamAAs = translateSequence(upstreamGenes[geneName], geneCoordinates[geneName][2])
            referenceAAs = translateSequence(referenceGenes[geneName], geneCoordinates[geneName][2])
            downstreamAAs = translateSequence(downstreamGenes[geneName], geneCoordinates[geneName][2])
            #position in gene depending on what strand the cds is on
            positionInGene = mutation[2] - geneCoordinates[geneName][0]
            aaPosition = extractPosition(geneCoordinates[geneName], positionInGene) - 1
            #extract codon 
            #transcribe gene
            if geneCoordinates[geneName][2] == "+":
                rc_ds = Seq("".join(downstreamGenes[geneName]))
                rc_us = Seq("".join(upstreamGenes[geneName]))
                rc_rs = Seq("".join(referenceGenes[geneName]))
            else:
                #reverse complement
                rc_ds = Seq("".join(downstreamGenes[geneName])).reverse_complement()
                rc_us = Seq("".join(upstreamGenes[geneName])).reverse_complement()
                rc_rs = Seq("".join(referenceGenes[geneName])).reverse_complement()
            reference_codon = rc_rs[(aaPosition * 3): (aaPosition * 3 + 3)]
            upstream_codon = rc_us[(aaPosition * 3): (aaPosition * 3 + 3)]
            downstream_codon = rc_ds[(aaPosition * 3): (aaPosition * 3 + 3)]
            #amino acid at substitution position
            reference_aa = referenceAAs[aaPosition]
            downstream_aa = downstreamAAs[aaPosition]
            upstream_aa = upstreamAAs[aaPosition]
            #check for start codon substitution, as we're using check_cds = False, non ATG start codons are translated into non methionine amino acids
            if aaPosition == 0:
                if upstream_codon == "GTG" or upstream_codon == "TTG":
                    upstream_aa = "M"
                if downstream_codon == "GTG" or downstream_codon == "TTG":
                    downstream_aa = "M"
                if reference_codon == "GTG" or reference_codon == "TTG":
                    reference_aa = "M"
            #check if amino acid changed due to the current substitution
            if upstream_aa != downstream_aa:
                #stop codon -> different codon
                if upstream_aa == "*":
                    mutation_type = "stop_lost"
                    downstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[downstream_aa]
                    upstream_aa3 = "*"
                    upstream_aa = "*"
                    impact = "HIGH"
                #stop codon inserted
                elif downstream_aa == "*":
                    mutation_type = "stop_gained"
                    upstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[upstream_aa]
                    downstream_aa3 = "*"
                    downstream_aa = "*"
                    impact = "HIGH"
                else:
                    #start codon -> different codon
                    if aaPosition == 0 and upstream_aa == "M" and downstream_aa != "M":
                        impact = "HIGH"
                        mutation_type = "start_lost"
                    #regular non-synonymous SNP
                    else:
                        impact = "MODERATE"
                        mutation_type = "missense_variant"
                    downstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[downstream_aa]
                    upstream_aa3 = SeqUtils.IUPACData.protein_letters_1to3[upstream_aa]
                    aa_change = "p.%s%s%s" % (upstream_aa3, aaPosition + 1, downstream_aa3)
                synonymous = False
            #synonymous SNP
            else:
                impact = "LOW"
                mutation_type = "synonymous_variant"
            #account for codons mutated at different positions on the same branch
            dict_key = "%s_%s_%s" % (clade, locus_tag, aaPosition + 1)
            if dict_key in aa_pos2effect:
                multi_codon_substitution = "multi_codon_substitution"
                #correct position to beginning of the codon
                if geneCoordinates[geneName][2] == "+":
                    pos = str(int(pos) - (positionInGene  % 3))
                    upstream_allele = upstream_codon
                    downstream_allele = downstream_codon
                else:
                    pos = str(geneCoordinates[geneName][1] - aaPosition * 3 - 2)
                    upstream_allele = upstream_codon.reverse_complement()
                    downstream_allele = downstream_codon.reverse_complement()
                aa_pos2effect[dict_key] = [node, pos, str(upstream_allele), str(downstream_allele), mutation_type, upstream_aa, downstream_aa, reference_aa, str(upstream_codon), str(downstream_codon),  str(reference_codon), impact, aa_change, multi_codon_substitution, locus_tag, pseudogene]
            else:
                aa_pos2effect[dict_key] = [node, pos, upstream_allele, downstream_allele, mutation_type, upstream_aa, downstream_aa, reference_aa, str(upstream_codon), str(downstream_codon), str(reference_codon), impact, aa_change, multi_codon_substitution, locus_tag, pseudogene]
            #If the mutation is nonsynonymous within any gene, remove it
            if synonymous == False:
                positionsToRemove.append(pos)
    #write effects to effect prediction out file
    for substitution in aa_pos2effect.values():
        effects.write("\t".join(substitution) + "\n")
    
    #Flag the positions that will not be included in the spectrum 
    synonymous_substitution_dict = {}
    if len(positionsToRemove) != 0:
        #Identify the unique set of positions
        uniquePositionsToRemove = list(set(positionsToRemove))
        #Flag synonymous substitutions
        for ele in sorted(uniquePositionsToRemove, reverse = True):
            synonymous_substitution_dict[ele] = False 
    
    return(synonymous_substitution_dict)

#Identifies the context of a mutation
def getContext(mutation, updatedReference):
    return(updatedReference[mutation[2] - 2], updatedReference[mutation[2]])

#Takes a base and returns the complement base
def complement(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"


