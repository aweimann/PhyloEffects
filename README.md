# Introduction
PhyloEffects is a piece of software to reconstruct variant effect annotations in their phylogenetic context.

## Installation
```conda install -c bioconda gffutils==0.12 treetime==0.11.2 cyvcf2==0.30.28  pyranges==0.0.129```   

## Example Usage
To run the example data navigate to the repository and run    
```python phyloeffects.py -a example_data/100_aln.fasta   -r example_data/Pseudomonas_aeruginosa_PAO1_v2.fa -t example_data/100_rescaled.nwk   -o example_data/phyloeffects -c  example_data/100_pos_mapping.txt   -g example_data/Pseudomonas_aeruginosa_PAO1_107_converted.gff```
