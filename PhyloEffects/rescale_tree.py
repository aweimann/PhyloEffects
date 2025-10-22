import argparse
import ete3 as ete 
from Bio import SeqIO

def parse():
    parser = argparse.ArgumentParser("rescale, prune and (midpoint) root gubbins tree")
    parser.add_argument("input_tree", help="input tree")
    parser.add_argument("--alignment", help="alignment")
    parser.add_argument("output_tree", help="output tree")
    parser.add_argument("--alignment_out", help="alignment pruned")
    parser.add_argument("--outgroup", help="outgroup to root tree")
    parser.add_argument("--rescale", help="scale by alignment length (if Gubbins tree)")
    parser.add_argument("--midpoint", help="if set, midpoint root tree", action = "store_true")
    parser.add_argument("--relabel", help="name internal nodes", action = "store_true")
    args = parser.parse_args()
    prune_tree(**vars(args))

def prune_tree(input_tree, output_tree, alignment, alignment_out, outgroup, midpoint, rescale, relabel):
    tree = ete.Tree(input_tree, format = 1)
    i = 1
    if relabel:
        for node in tree.traverse():
            if node.is_leaf():
                continue
            else:
                node.name = "internal{}".format(i)
                i+=1
    if rescale:
        with open(alignment_out, 'a') as f:
            with open(alignment) as aln:
                record_iter = SeqIO.parse(aln, 'fasta')
                for entry in record_iter:
                    if outgroup and entry.id != outgroup:
                        SeqIO.write(entry, f, "fasta")
                    else:
                        SeqIO.write(entry, f, "fasta")
        with open(alignment) as f:
            f.readline()
            scaling_factor = len(f.readline())
        for node in tree.traverse():
            node.dist = node.dist/scaling_factor
    if outgroup is not None:
        for node in tree.traverse():
            if node.name == outgroup:
                tree.set_outgroup(node)
        for node in tree.get_children():
            if node.name.startswith("internal"):
                tree.set_outgroup(node)
                node.name = "root_branch"
        terminal_nodes = []
        for node in tree.iter_leaves():
            if node.name != outgroup:
                terminal_nodes.append(node)
        #prune outgroup
        tree.prune(terminal_nodes)
    else:
         if midpoint:
            mp = tree.get_midpoint_outgroup()
            tree.set_outgroup(mp)

    tree.write(outfile = output_tree, format = 1)

if __name__ == "__main__":
    parse()
