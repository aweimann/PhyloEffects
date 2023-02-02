#Runs the MutTui pipeline on an input tree and alignment

import os
import argparse
import isvalid
from treetime import run_treetime
from Bio import AlignIO, Phylo
from branch_labelling import labelAllBranches
import reconstruct_spectrum as rs
import gff_conversion

from __init__ import __version__

#Parse command line options
def get_options():
    description = "Run the MutTui pipeline on a given alignment and tree"

    parser = argparse.ArgumentParser(description = description)

    #Options for input and output files
    io_opts = parser.add_argument_group("Input/output")
    io_opts.add_argument("-a",
                        "--alignment",
                        dest = "alignment",
                        required = True,
                        help = "Input fasta alignment",
                        type = argparse.FileType("r"))
    io_opts.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick tree file",
                        type = argparse.FileType("r"))
    io_opts.add_argument("-o",
                        "--out_dir",
                        dest = "output_dir",
                        required = True,
                        help = "Location of output directory, should already be created and ideally be empty",
                        type = lambda x: isvalid.is_valid_folder(parser, x))
    io_opts.add_argument("-r",
                        "--reference",
                        dest = "reference",
                        help = "Reference genome sequence containing all sites, used to get context of mutations, not required if using --all_sites",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-c",
                        "--conversion",
                        dest = "conversion",
                        help = "Conversion file for alignment position to genome position, used to get context of mutations, not required if using --all_sites",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-g",
                        "--gff",
                        dest = "gff",
                        help = "GFF reference containing gene coordinates in reference sequence. Used to split mutations into transcription strands " + 
                        "and identify synonymous mutations when --synonymous is used",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-to",
                        "--treetime_out",
                        dest = "treetime_out",
                        help = "The location of the directory containing treetime output files from ancestral reconstruction. Only " + 
                        "used if option --start_from_treetime is specified, in which case the output files in this directory will " + 
                        "be used to calculate the mutational spectrum",
                        default = None)

    #Options to include in treetime reconstruction
    treetime = parser.add_argument_group("treetime command")
    treetime.add_argument("--add_treetime_cmds",
                        dest = "add_treetime_cmds",
                        help = "Additional options to supply to treetime (these are not checked). Supply these together in quotes",
                        type = str,
                        default = None)
    
    #Other options
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen, MutTui will output an RNA mutational spectrum",
                        action = "store_true",
                        default = False)
    parser.add_argument("--all_sites",
                        dest = "all_sites",
                        help = "Specify that the alignment contains all sites, in which case a reference genome does not need to be provided",
                        action = "store_true",
                        default = False)
    parser.add_argument("--filter",
                        dest = "filter",
                        help = "Converts gaps to Ns in the sequence alignment. treetime will reconstruct gaps onto the tree. " + 
                        "This is fine when there's not many gap but if the alignment contains many gaps, the annotated tree from " + 
                        "treetime becomes very large and time consuming to read into python. Ns are not reconstructed by treetime. " + 
                        "By converting gaps to Ns, it reduces the number of reconstructed mutations and greatly speeds up run time. " + 
                        "Use this option if your alignment contains many gaps",
                        action = "store_true",
                        default = False)
    parser.add_argument("--start_from_treetime",
                        dest = "start_from_treetime",
                        help = "Use this option to start with treetime output and so skip inference of ancestral mutations. Use this " + 
                        "if you have already run treetime. The directory containing the treetime output files needs to be provided with -to",
                        action = "store_true",
                        default = False)
    parser.add_argument("--synonymous",
                        dest = "synonymous",
                        help = "Use non-coding and synonymous mutations only to calculate the mutational spectrum. A GFF file will need to be " + 
                        "provided with option -g which will be used to identify genes",
                        action = "store_true",
                        default = False)
    parser.add_argument("--version",
                        action = "version",
                        version = "%(prog)s " + __version__)
    
    args = parser.parse_args()
    return(args)

def main():
    args = get_options()

    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #Open output files
    outMutationsNotUsed = open(args.output_dir + "mutations_not_included.csv", "w")
    outMutationsNotUsed.write("Mutation_in_alignment,Mutation_in_genome,Branch,Reason_not_included\n")

    outAllMutations = open(args.output_dir + "all_included_mutations.csv", "w")
    outAllMutations.write("Mutation_in_alignment,Mutation_in_genome,Substitution,Branch\n")

    if not args.start_from_treetime:
        print("Running treetime ancestral reconstruction to identify mutations")

        #Check if the alignment is to be converted so gaps become Ns. If so, run the conversion
        #and run treetime on the new alignment
        if args.filter:
            rs.change_gaps_to_Ns(args.alignment, args.output_dir)
            run_treetime(open(args.output_dir + "gaps_to_N_alignment.fasta"), args.tree, args.output_dir, args.add_treetime_cmds)
        else:
            pass
            #Run treetime on the input alignment and tree with any provided options
            run_treetime(args.alignment, args.tree, args.output_dir, args.add_treetime_cmds)
    
        print("treetime reconstruction complete. Importing alignment from reconstruction and tree")

        #Import the alignment from treetime
        alignment = AlignIO.read(args.output_dir + "ancestral_sequences.fasta", "fasta")
        
    else:
        #Import the already calculated alignment from treetime
        print("Importing alignment from reconstruction and tree")
        #Make sure the trailing forward slash is present in output directory
        args.treetime_out = os.path.join(args.treetime_out, "")

        #Import the alignment from treetime
        alignment = AlignIO.read(args.treetime_out + "ancestral_sequences.fasta", "fasta")
    
    #Import the original unlabelled tree
    tree = Phylo.read(args.tree.name, "newick")
    #Ladderize the tree so the branches are in the same order as the treetime tree
    tree.ladderize()

    print("Alignment and tree imported. Reconstructing spectrum")

    #Check if a GFF file is needed, if so read it in and process it
    if args.synonymous:
        if not args.gff:
            raise RuntimeError("GFF file needs to be provided with -g when using --strand_bias or --synonymous")
        else:
            geneCoordinates, positionGene = gff_conversion.convertGFF(args.gff.name)
            with open(args.output_dir + "gene_annotation.txt", 'w') as f:
                #header
                f.write("start\tend\tstrand\tlocus_tag\tfeature\n")
                for value_list in geneCoordinates.values():
                    f.write("\t".join([str(i) for i in value_list]) + "\n")

    labelledTree, treeLabels = labelAllBranches(tree)

    #The 4 nucleotides, used to check if mutated, upstream and downstream bases are nucleotides
    nucleotides = ["A","C","G","T"]

    #Convert the positions in the alignment to genome positions, if --all_sites specified the positions will be the same
    if args.all_sites:
        positionTranslation = rs.allSitesTranslation(alignment)
    else:
        positionTranslation = rs.convertTranslation(args.conversion)
    
    #Extracts mutations to a dictionary from the branch_mutations.txt file
    if not args.start_from_treetime:
        branchMutationDict = rs.getBranchMutationDict(args.output_dir + "branch_mutations.txt", positionTranslation)
    else:
        branchMutationDict = rs.getBranchMutationDict(args.treetime_out + "branch_mutations.txt", positionTranslation)

    #Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the alignment are assumed
    #and the root sequence from the ancestral reconstruction is used
    referenceSequence = rs.getReference(args.reference, args.all_sites, alignment, positionTranslation)
    referenceLength = len(referenceSequence)
    
    #Iterate through the branches, get the category of the branch, identify the contextual mutations, add to the corresponding spectrum
    if args.synonymous:
        effects = open(args.output_dir + "variant_effect_predictions.txt", "w")
        effects.write("\t".join(["node", "pos", "upstream_allele", "downstream_allele", "mutation_type", "upstream_aa", "downstream_aa", "reference_aa", "upstream_codon", "downstream_codon", "reference_codon", "impact", "aa_change",  "multi_codon_substitution", "locus_tag", "pseudogene"]) + "\n")
        effects.close()
        #get reference
        ref = AlignIO.read(args.reference.name, "fasta")
        #Extract the sequence of the reference and ensure it is uppercase
        refSeq = ref[0].seq.upper()
    for clade in labelledTree.find_clades():
        #Check if there are mutations along the current branch, only need to analyse branches with mutations
        if clade.name in branchMutationDict:

            #Extract the mutations along the branch. This will be None if there are no mutations but treetime has still added the mutation comment
            branchMutations = branchMutationDict[clade.name]

            #If using --branch_specific, check if the branch contains at least -bm mutations, if so
            #add the branch to spectraDict. If not, set branchCategory to None so it won't be analysed

            #Check if the branch has a category, will be None if the branch is a transition between categories
            #and option --include_all_branches is not specified
            #Extract double substitutions, remove mutations at the ends of the genome or not involving 2 nucleotides
            double_substitution_dict, doubleSubstitutions = rs.filterMutations(branchMutations, clade, nucleotides, referenceLength, outMutationsNotUsed)

            #Update the reference sequence to get the current context
            updatedReference = rs.updateReference(tree, clade, branchMutationDict, referenceSequence)

            #Check if only synonymous mutations should be included, if so filter the mutations
            if args.synonymous:
                synonymous_substitution_dict = rs.extractSynonymous(clade, branchMutations, updatedReference, refSeq, geneCoordinates, positionGene, args.output_dir)







if __name__ == "__main__":
    main()
