import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
import argparse

def parse():
    parser = argparse.ArgumentParser("extract reconstructed variants from input pseudoalignment")
    parser.add_argument("variants",  help="variant annotation")
    parser.add_argument("pseudoalignment",  help="pseudoalignemtn will all target samples")
    parser.add_argument("sample2variant_out",  help="output table for variants per sample")
    parser.add_argument("--is_one_hot", action = "store_true",  help="output table for variants per sample")
    parser.add_argument("--group_impute", help="impute missing calls by majority allele per group")
    parser.add_argument("--is_coding", action = "store_true",  help="do coding variants at the amino acid level")
    args = parser.parse_args()
    extract(**vars(args))


#TODO exclude non-coding genes (done, needs testing)
#TODO ignore indels (done, needs testing)
#TODO N as character for ambigious amino acid not ideal as also used for amino acid 
#TODO join annotated variants and annotated variants
#TODO codon start position rather than nucleotide position (done)
#TODO overlapping features (done)
#TODO away from the reference encoding (done, but currently hard coded)

NUCLEOTIDE_AMB = 'N'
AMINO_ACID_AMB = 'X'

missing_char = NUCLEOTIDE_AMB

def get_amino_acids(var2geno, pos_dict, variants, pseudoalignment):
    with open(variants, 'r') as f:
        header = f.readline().strip().split('\t')
        key2col = dict([(key, i) for key, i in zip(header, range(len(header)))])
        for l in f:
            fields = l.split('\t')
            impact = fields[key2col['impact']]
            PAO1 = fields[key2col['PAO1']]
            ref = fields[key2col['ref']]
            alt = fields[key2col['alt']]
            mutation_type = fields[key2col['mutation_type']]
            if len(ref) == len(alt) and ((impact == "MODERATE" and mutation_type != "non-coding") or impact == "HIGH"):
                us_codon = fields[key2col['upstream_codon']]
                ds_codon = fields[key2col['downstream_codon']]
                codon_pos = None 
                if us_codon[0] != ds_codon[0]:
                    codon_pos = 0
                elif us_codon[1] != ds_codon[1]:
                    codon_pos = 1
                else:
                    codon_pos = 2
                pos = int(fields[key2col['pos']])
                strand = fields[key2col['Strand']]
                if fields[key2col['Strand']]  == '+':
                    pos = pos - codon_pos
                else:
                    pos = pos + codon_pos
                ref = fields[key2col['upstream_aa']] 
                alt = fields[key2col['downstream_aa']]
                pos_dict[(int(pos), PAO1)] = (ref, alt, codon_pos, strand) 
                var2geno[(int(pos), PAO1)] = []
    #iterate through pseudoalignemnt
    j = 0
    sample2id = {}
    isols = []
    with open(pseudoalignment) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sample2id[record.id] = j
            isols.append(record.id)
            print(record.id)
            for var, info in pos_dict.items():
                pos, feature_id = var
                if info[3] == "+": 
                    codon = record.seq[pos - 1: pos + 2]
                    aa = str(codon.translate())

                else:
                    codon = record.seq[pos - 3 : pos]
                    reverse_codon = codon.reverse_complement()
                    aa = str(reverse_codon.translate())
                var2geno[(pos, feature_id)].append(aa)
            j += 1
            #DEBUG
            if j == 5:
                pass 
    return sample2id, var2geno, pos_dict, isols 
    

def extract(variants, pseudoalignment, is_one_hot, group_impute, is_coding, sample2variant_out):
    """main function"""
    #determine missing char
    if is_coding:
        global missing_char
        missing_char = AMINO_ACID_AMB
    pos_dict = {}
    var2geno = {}
    #read in variant annotation
    if not is_coding:
        isols = []
        with open(variants, 'r') as f:
            f.readline()
            for l in f:
                fields = l.split('\t')
                pos, ref, alt = fields[2:5]
                #ignore indels
                if len(ref) == len(alt):
                    pos_dict[(int(pos), "")] = (ref, alt) 
                    var2geno[(int(pos), "")] = []
        #iterate through pseudoalignemnt
        j = 0
        sample2id = {}
        with open(pseudoalignment) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                print(record.id, j)
                sample2id[record.id] = j
                j += 1
                isols.append(record.id)
                #0 vs 1-indexing
                for pos, _ in pos_dict:
                    var2geno[(pos, _)].append(record.seq[pos - 1])
                if j == 5:
                    pass 
                     
    else:
        sample2id, var2geno, pos_dict, isols = get_amino_acids(var2geno, pos_dict, variants, pseudoalignment)
    print(var2geno.keys(), pos_dict)

    #read in group data if given
    if group_impute:
        groups = []
        sample2group = {}
        group2samples = {}
        with open(group_impute) as f:
            f.readline()
            sample2group['PAO1'] = "reference"
            for line in f:
                sample, group = line.strip().split()
                if sample in sample2id:
                    sample2group[sample] = group
        for sample, group in sample2group.items():
            if group in group2samples:
                group2samples[group].append(sample)
            else:
                group2samples[group] = [sample]

    id2sample = dict([(i, sample) for sample, i in sample2id.items()])
    if group_impute:
        groups = [sample2group[id2sample[i]] for i in range(len(sample2id))]
    #write output variant sample file
    with open(sample2variant_out, 'w') as out:
        out.write("pos\t{}\n".format("\t".join(isols),"\n"))
        for pos, feature_id in pos_dict:
            uq_chars = set() 
            for char in var2geno[(pos, feature_id)]:
                uq_chars.add(char)
            if len(uq_chars) == 2 and missing_char not in uq_chars or len(uq_chars) > 2:
                if group_impute:
                    var2geno[(pos, feature_id)] = impute_var(var2geno[(pos, feature_id)], groups, sample2id, group2samples)
                if is_one_hot:
                    encoding = one_hot(var2geno[(pos, feature_id)], pos, uq_chars, is_coding, feature_id) 
                    out.write(encoding)
                else:
                    out.write("{}{}{}\t{}\n".format(pos, "_" if is_coding else "", feature_id, "\t".join(var2geno[(pos, feature_id)])))

def impute_var(var, groups, sample2id, group2samples):
    """do groupwise or normal imputation"""
    var_imputed = []
    for char, i in zip(var, range(len(var))):
        if char == missing_char:
            #try groupwise imputation
            #get allele for other members in this group
            group = groups[i]
            group_alleles = [var[sample2id[sample]] for sample in group2samples[groups[i]]]
            maj_allele = get_majority_allele(group_alleles)
            #check if we succesful and do imputation across all samples
            if maj_allele == missing_char:
                maj_allele = get_majority_allele(var)
            else:
                var, group2samples, maj_allele
            var_imputed.append(maj_allele)
        else:
            var_imputed.append(char)
    return var_imputed
            
    #try group impute


def get_majority_allele(var):
    """get the majority allele"""
    #get majority allele
    allele2count = {}
    for char in var:
        if char not in allele2count:
            allele2count[char] = 1
        else:
            allele2count[char] += 1
    all_sorted = sorted((value, key) for (key,value) in allele2count.items())
    if all_sorted[-1][1] == missing_char:
        if len(all_sorted) > 1: 
            return all_sorted[-2][1] 
        else:
            return all_sorted[-1][1]
    else:
        return all_sorted[-1][1]
    

def one_hot(chars, pos, uq_chars, is_coding, feature_id):
    """do one hot encoding"""
    out_encoding = []
    i = 0
    uq_chars = list(uq_chars)
    #filter missing char and reference char
    #this will make sure the encoding is always away from the reference
    #PAO1 needs to be at index 0 for this to work
    uq_chars = [i for i in uq_chars if i != missing_char and i != chars[0]]
    for target in uq_chars:
        cur_encoding = []
        for char in chars:
            if char == target:
                cur_encoding.append("1")
            else:
                cur_encoding.append("0")
        out_encoding.append("{}{}{}_{}\t{}\n".format(pos, "_" if is_coding else "", feature_id, target, "\t".join(cur_encoding)))
    out_encoding = "".join(out_encoding) 
    return out_encoding


if __name__ == "__main__":
    parse()
