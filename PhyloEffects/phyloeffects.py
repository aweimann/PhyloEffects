# Runs the PhyloEffects pipeline on an input tree and alignment

import argparse
import os
import sys

from Bio import SeqIO, AlignIO, Phylo

import gff_conversion
import isvalid
import reconstruct_variant_effects as rs
from __init__ import __version__
from run_treetime import run_treetime


# Parse command line options
def get_options(argv=None):
    description = "Run the PhyloEffects pipeline on a given alignment and tree"

    parser = argparse.ArgumentParser(description=description)

    # Options for input and output files
    io_opts = parser.add_argument_group("Input/output")
    io_opts.add_argument("-a",
                         "--alignment",
                         dest="alignment",
                         help="Input fasta alignment",
                         type=argparse.FileType("r"))
    io_opts.add_argument("-v",
                         "--vcf",
                         dest="vcf",
                         required=False,
                         help="Input VCF file",
                         type=argparse.FileType("r"))
    io_opts.add_argument("-vt",
                         "--variant-table",
                         dest="variant_table",
                         required=False,
                         help="variant table file",
                         type=argparse.FileType("r"))
    io_opts.add_argument("-t",
                         "--tree",
                         dest="tree",
                         required=False,
                         help="Newick tree file",
                         type=argparse.FileType("r"))
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="Location of output directory, should already be created and ideally be empty",
                         type=lambda x: isvalid.is_valid_folder(parser, x))
    io_opts.add_argument("--output_prefix",
                         dest="output_prefix",
                         required=False,
                         default="",
                         help="Optional prefix for variant_effect_predictions.txt (becomes <prefix>.variant_effect_predictions.txt). If not specified, defaults to empty string so output file is variant_effect_predictions.txt")
    io_opts.add_argument("-r",
                         "--reference",
                         dest="reference",
                         help="Reference genome sequence containing all sites, used to get context of mutations, "
                              "not required if using --all_sites",
                         type=argparse.FileType("r"),
                         default=None)
    io_opts.add_argument("-c",
                         "--conversion",
                         dest="conversion",
                         help="Conversion file for alignment position to genome position, used to get context of "
                              "mutations, not required if using --all_sites",
                         type=argparse.FileType("r"),
                         default=None)
    io_opts.add_argument("-g",
                         "--gff",
                         required=True,
                         dest="gff",
                         help="GFF reference containing gene coordinates in reference sequence. Used to split "
                              "mutations into transcription strands " +
                              "and identify synonymous mutations when --synonymous is used",
                         type=argparse.FileType("r"),
                         default=None)
    io_opts.add_argument("-to",
                         "--treetime_out",
                         dest="treetime_out",
                         help="The location of the directory containing treetime output files from ancestral "
                              "reconstruction. Only " +
                              "used if option --start_from_treetime is specified, in which case the output files in "
                              "this directory will " +
                              "be used to calculate the mutational spectrum",
                         default=None)

    # Options to include in treetime reconstruction
    treetime = parser.add_argument_group("treetime command")
    treetime.add_argument("--add_treetime_cmds",
                          dest="add_treetime_cmds",
                          help="Additional options to supply to treetime (these are not checked). Supply these "
                               "together in quotes",
                          type=str,
                          default=None)

    # Other options
    parser.add_argument("--all_sites",
                        dest="all_sites",
                        help="Specify that the alignment contains all sites, in which case a reference genome does "
                             "not need to be provided",
                        action="store_true",
                        default=False)
    parser.add_argument("--filter",
                        dest="filter",
                        help="Converts gaps to Ns in the sequence alignment. treetime will reconstruct gaps onto the "
                             "tree. " +
                             "This is fine when there's not many gap but if the alignment contains many gaps, "
                             "the annotated tree from " +
                             "treetime becomes very large and time consuming to read into python. Ns are not "
                             "reconstructed by treetime. " +
                             "By converting gaps to Ns, it reduces the number of reconstructed mutations and greatly "
                             "speeds up run time. " +
                             "Use this option if your alignment contains many gaps",
                        action="store_true",
                        default=False)
    parser.add_argument("--start_from_treetime",
                        dest="start_from_treetime",
                        help="Use this option to start with treetime output and so skip inference of ancestral "
                             "mutations. Use this " +
                             "if you have already run treetime. The directory containing the treetime output files "
                             "needs to be provided with -to",
                        action="store_true",
                        default=False)
    parser.add_argument("--version",
                        action="version",
                        version="%(prog)s " + __version__)

    args = parser.parse_args(argv)
    # if treetime directory is not set; set to output directory
    if args.treetime_out is None:
        args.treetime_out = args.output_dir
    return args


def run_with_args(args):

    # Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")


    if not args.start_from_treetime:
        print("Running treetime ancestral reconstruction to identify mutations")

        if args.filter:
            raise NotImplementedError("The --filter option has been deprecated. Gaps are ignored by default")
        else:
            # Run treetime on the input alignment and tree with any provided options
            if args.tree:
                run_treetime(args.alignment, args.tree, args.output_dir, args.add_treetime_cmds)
                print("treetime reconstruction complete. Importing alignment from reconstruction and tree")

                # Import the alignment from treetime
                alignment = AlignIO.read(args.output_dir + "ancestral_sequences.fasta", "fasta")

    else:
        if args.tree:
            # Import the already calculated alignment from treetime
            print("Importing alignment from reconstruction and tree")
            # Make sure the trailing forward slash is present in output directory
            args.treetime_out = os.path.join(args.treetime_out, "")

            # Import the alignment from treetime
            alignment = AlignIO.read(args.treetime_out + "ancestral_sequences.fasta", "fasta")
        else:
            alignment = AlignIO.read(args.reference.name, format='fasta')

    # Import the original unlabelled tree
    if args.tree:
        tree = Phylo.read(os.path.join(args.treetime_out, "annotated_tree.nexus"), "nexus")
        for node in tree.find_clades():
            if node.name is None and node.confidence is not None:
                node.name = node.confidence
        # Ladderize the tree so the branches are in the same order as the treetime tree
        tree.ladderize()

    print("Alignment and tree imported. Reconstructing variant effects")

    gene_coordinates, position_gene, attribute_header = gff_conversion.convertGFF(args.gff.name)
    with open(args.output_dir + "gene_annotation.txt", 'w') as f:
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


    #labelled_tree, tree_labels = labelAllBranches(tree)

    # The 4 nucleotides, used to check if mutated, upstream and downstream bases are nucleotides
    nucleotides = ["A", "C", "G", "T"]

    # get reference
    ref = next(SeqIO.parse(args.reference.name, "fasta"))
    # Convert the positions in the alignment to genome positions, if --all_sites specified the positions will be the
    # same
    if args.all_sites:
        alignment = AlignIO.read(args.treetime_out + "ancestral_sequences.fasta", "fasta")
        position_translation = rs.all_sites_translation(alignment)
    elif args.alignment:
        position_translation = rs.convert_translation(args.conversion)
    else:
        position_translation = rs.all_sites_translation([ref])

    # Extract the sequence of the reference and ensure it is uppercase
    ref_seq = ref.seq.upper()
    # Extracts mutations to a dictionary from the branch_mutations.txt file
    if args.tree:
        if args.start_from_treetime:
            branch_mutation_dict = rs.get_branch_mutation_nexus_dict(args.treetime_out + "annotated_tree.nexus",
                                                                     position_translation)
        else:
            branch_mutation_dict = rs.get_branch_mutation_nexus_dict(args.output_dir + "annotated_tree.nexus",
                                                                     position_translation)
    elif args.alignment:
        branch_mutation_dict, clades = rs.get_mutations_from_alignment(args.alignment, ref_seq, position_translation)
    elif args.variant_table:
        branch_mutation_dict, clades = rs.get_mutations_from_table(args.variant_table)
    else:
        if not args.vcf:
            raise Exception()
        else:
            branch_mutation_dict, clades = rs.get_mutations_from_vcf(args.vcf)
    print("Alignment/variants processed")


    # Iterate through the branches
    if args.output_prefix:
        output_filename = args.output_prefix + ".variant_effect_predictions.txt"
    else:
        output_filename = "variant_effect_predictions.txt"
    effects = open(args.output_dir + output_filename, "w")
    effects.write("\t".join(
        ["pos", "upstream_allele", "downstream_allele", "mutation_type", "upstream_aa", "downstream_aa",
         "reference_aa", "upstream_codon", "downstream_codon", "reference_codon", "impact", "aa_change",
         "multi_codon_substitution", "locus_tag", "pseudogene", "samples"]) + "\n")
    # Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the
    # alignment are assumed
    # and the root sequence from the ancestral reconstruction is used
    if not args.tree:
        # If a tree is not provided, the reference sequence is assumed to be the first sequence in the alignment
        if not args.vcf and not args.variant_table:
            args.alignment.seek(0, 0)
        reference_sequence = next(SeqIO.parse(args.reference.name, "fasta"))
        reference_sequence = reference_sequence.seq.upper()
    else:
        if args.reference:
            reference_sequence = rs.get_reference(args.reference, args.all_sites, alignment, position_translation)
        else:
            sys.exit("A reference sequence must be provided with -r when using --tree")
    reference_length = len(reference_sequence)
    clade2name = {}
    if args.tree:
        clades = tree.find_clades()
        for clade in clades:
            clade2name[clade] = clade.name
        clades = tree.find_clades()
    else:
        for clade in clades:
            clade2name[clade] = clade

    variant_effect2clades = {}
    for clade in clades:
        # Check if there are mutations along the current branch, only need to analyse branches with mutations
        clade_name = clade2name[clade]
        if clade_name in branch_mutation_dict and len(branch_mutation_dict[clade_name]) != 0:
            # Extract the mutations along the branch. This will be None if there are no mutations but treetime has
            # still added the mutation comment
            branch_mutations = branch_mutation_dict[clade2name[clade]]
            # Update the reference sequence to get the current context
            if args.tree:
                updated_reference = rs.update_reference(tree, clade, branch_mutation_dict, reference_sequence)
            else:
                updated_reference = reference_sequence
            # infer variant effects
            rs.reconstruct_effects(clade, branch_mutations, updated_reference, ref_seq, variant_effect2clades,
                                   gene_coordinates, position_gene, args.output_dir, args.output_prefix)


    for variant_effect, clades in variant_effect2clades.items():
        effects.write(variant_effect.to_string() + "\t" + ",".join(clades) + "\n")


def main(argv=None):
    args = get_options(argv)
    run_with_args(args)


if __name__ == "__main__":
    main()
