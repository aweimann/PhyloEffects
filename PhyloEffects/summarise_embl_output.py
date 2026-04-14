import argparse

def parse():
    parser = argparse.ArgumentParser("parse embl gubbisn output and generate table of substitutions")
    parser.add_argument("embl",  help="gubbins embl output listing substitutions on the tree")
    parser.add_argument("out_table", help="output summary table")
    args = parser.parse_args()
    extract(**vars(args))


def extract(embl, out_table):
    out = open(out_table, 'w')
    out.write("pos\tnode\tparent\tchild\n")
    with open(embl, 'r') as f:
        for l in f:
            if "variation" in l:
                pos = l[21:].strip()
                out.write(pos + "\t")
            elif "/node" in l:
                branch = l.split("=")[1].strip('\n"').split("->")[1]
                out.write(branch + "\t")
            elif "/parent_base" in l:
                parent = l.split("=")[1].strip('\n"')
                out.write(parent + "\t")
            elif "/replace" in l:
                child = l.split("=")[1].strip('\n"')
                out.write(child + "\n")

    out.close()

if __name__ == "__main__":
    parse()