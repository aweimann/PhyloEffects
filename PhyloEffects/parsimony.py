import ete3 as ete
import pandas as pd
import argparse
import os

def parse():
    parser = argparse.ArgumentParser("reconstruct gain/losses using max parsimony")
    parser.add_argument("tree", help="input phylogenetic tree")
    parser.add_argument("table", help="input binary table (samples as rows, characters as columns)")
    parser.add_argument("outdir", help="isolate group specific nodes out")
    parser.add_argument("--tree_is_named", action = "store_true", help = "set if tree has internal node names")
    parser.add_argument("--is_transposed", action = "store_true", help = "set to indicate input is genes (rows) X samples (columns) like in Panaroo presence/absence")
    parser.add_argument("--prefix", type=str, required=True, help="prefix for output files")
    args = parser.parse_args()
    parsimony(**vars(args))

def parsimony(tree, table, outdir, tree_is_named, is_transposed, prefix):
    t = ete.Tree(tree, format=1)
    i = 0
    if not tree_is_named:
        for n in t.traverse("preorder"):
            if not n.is_leaf():
                i += 1
                n.name = "N%s" % i
        t.write(outfile = f"{outdir}/{prefix}.named_tree.nwk", format = 1)

    table = pd.read_csv(table, index_col = 0, sep = "\t")
    if is_transposed:
        table = table.T
    table.index = [i.replace(".velvet", "") for i in table.index]
    table.index = [i.replace(".spades", "") for i in table.index]
    table.fillna(float('inf'), inplace=True)
    node_names = [n.name for n in t.traverse("preorder")]
    out_table = pd.DataFrame([[float('inf')] * table.shape[1] for _ in range(len(node_names))],
                             index=node_names, columns=table.columns)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    events = open(f"{outdir}/{prefix}.indel_events.txt", 'w')
    events.write("variant_id\tnode\tparent_node\tnode_state\tparent_node_state\n")
    for c in table.columns:
        if all(pd.isnull(table.loc[:, c])) | all(table.loc[:, c] == float('inf')):
            out_table.loc[:, c] = np.nan
        else:
            node2state = down_pass(t, table.loc[:, c])
            reconstruction = up_pass(node2state, t, c, events)
            out_table.loc[reconstruction.keys(), c] = list(reconstruction.values())
    events.close()

    out_table.to_csv(f"{outdir}/{prefix}.ancestral_states.txt", sep = "\t")

def up_pass(node2state, tree, c, events):
    reconstruction = {}
    for n in tree.traverse("preorder"):
        if n == tree.get_tree_root():
            states = list(node2state[n.name])
            if len(states) > 1:
                reconstruction[n.name] = states[0]
            else:
                reconstruction[n.name] = states.pop()
        else:
            up_state = reconstruction[n.up.name]
            n_states = list(node2state[n.name])
            if up_state in node2state[n.name]:
                reconstruction[n.name] = up_state
            else:
                if float('inf') in node2state[n.name]:
                    reconstruction[n.name] = float('nan')
                elif len(n_states) == 0:
                    reconstruction[n.name] = float('nan')
                else:
                    reconstruction[n.name] = n_states[0]
                    events.write(f"{c}\t{n.name}\t{n.up.name}\t{reconstruction[n.up.name]}\t{reconstruction[n.name]}\n")
    return reconstruction


def down_pass(tree, annotation):
    node2state = {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            node2state[n.name] = set([annotation.loc[n.name]])
        else:
            node2state[n.name] = set.intersection(*[node2state[i.name] for i in n.children])
            if len(node2state[n.name]) == 0:
                node2state[n.name] = set.union(*[node2state[i.name] for i in n.children])
            node2state[n.name].discard(float('inf'))
    return node2state

        

    
if __name__ == "__main__":
    parse()
