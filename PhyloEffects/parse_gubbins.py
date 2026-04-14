import argparse

def parse():
    parser = argparse.ArgumentParser("extract table with recombination interval and tree node from Gubbins output")
    parser.add_argument("gubbins",  help="gubbins output")
    parser.add_argument("out_table", help="output summary table")
    args = parser.parse_args()
    extract(**vars(args))

def extract(gubbins, out_table):
    out_table = open(out_table, 'w')
    out_table.write("start\tstop\tnode\n")
    with open(gubbins, 'r') as f:
        for l in f:
            if not l.startswith("#"):
                fields = l.split("\t")
                start = fields[3]
                end = fields[4]
                attributes = dict([s.split("=") for s in fields[8].strip("\n;").split(";")])
                node = attributes["node"]
                node = node.strip('"').split("->")[1]
                out_table.write("%s\t%s\t%s\n" % (start, end, node))

if __name__ == "__main__":
    parse()