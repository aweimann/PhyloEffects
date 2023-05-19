class VariantEffect:
    def __init__(self, position, upstream_allele, downstream_allele, node):
        self.node = node
        self.position = position
        self.upstream_aa = ""
        self.downstream_aa = ""
        self.reference_aa = ""
        self.upstream_codon = ""
        self.downstream_codon = ""
        self.reference_codon = ""
        self.impact = ""
        self.aa_change = ""
        self.multi_codon_substitution = ""
        self.upstream_allele = upstream_allele
        self.downstream_allele = downstream_allele
        self.pseudogene = ""

    def to_string(self):
        return '\t'.join([str(i) for i in [self.node, self.position, self.upstream_allele, self.downstream_allele,
                                                    self.mutation_type, self.upstream_aa, self.downstream_aa,
                                                    self.reference_aa, self.upstream_codon, self.downstream_codon,
                                                    self.reference_codon, self.impact, self.aa_change,
                                                    self.multi_codon_substitution, self.locus_tag, self.pseudogene]])

    def __key(self):
        return self.position, self.upstream_allele, self.downstream_allele, self.locus_tag

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, VariantEffect):
            return self.__key() == other.__key()
        return NotImplemented