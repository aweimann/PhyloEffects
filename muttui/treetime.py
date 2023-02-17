# Runs treetime ancestral reconstruction on an alignment and tree

import os
import subprocess
from Bio import AlignIO


# Runs treetime ancestral reconstruction
def run_treetime(alignment, tree, output_dir, add_treetime_cmds):
    ####This is the command for normal treetime
    ####Revert to this once the branch_mutations.txt file is added to treetime
    # cmd = "treetime ancestral "

    ####Used to run a specific version of treetime until branch_mutations.txt
    ####is added into treetime
    ####Remove this once the file is incorporated
    cmd = "treetime ancestral "

    # Check if a model has been specified, if not use --gtr infer
    if add_treetime_cmds == None:
        cmd += "--gtr infer"
    elif "--gtr" not in add_treetime_cmds:
        cmd += "--gtr infer"
    else:
        cmd += add_treetime_cmds

    cmd += " --aln " + alignment.name
    cmd += " --tree " + tree.name
    cmd += " --outdir " + output_dir

    subprocess.run(cmd, shell=True, check=True, executable='/bin/sh')

    # Check if treetime completed successfully
    if (os.stat(output_dir + "ancestral_sequences.fasta").st_size == 0) or (
            os.stat(output_dir + "annotated_tree.nexus") == 0):
        raise Exception("treetime did not complete successfully")


# treetime reconstructs mutations at gap sites, which are later removed by MutTui
# This can be a problem if there are large numbers of gap sites as the resulting annotated nexus tree
# takes a very long time to read into python
# Ns are not reconstructed by treetime. This function replaces gaps with Ns in the alignment
# so they will not be reconstructed. As gaps and Ns are both treated as missing data
# in phylogenetic methods, this will not alter any inferences but will remove large numbers
# of gaps from the infered mutations
def change_gaps_to_Ns(alignmentFile, output_dir):
    alignment = AlignIO.read(alignmentFile.name, "fasta")

    new_a = open(output_dir + "gaps_to_N_alignment.fasta", "w")

    for seq in alignment:
        new_a.write(">" + seq.id + "\n" + str(seq.seq).replace("-", "N") + "\n")

    new_a.close()
