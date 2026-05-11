#!/usr/bin/env python3
"""
Filter indels from VCF files and aggregate across bacterial genomes.
Groups identical variants and tracks which samples have each variant.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def get_sample_name(vcf_file):
    """Extract sample name from VCF header (last column of #CHROM line)."""
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                if len(fields) > 9:
                    return fields[9]  # Sample name is after FORMAT field
    return None


def parse_vcf_indels(vcf_file, sample_name):
    """
    Parse VCF file and extract indels with quality information.

    Yields tuples of (chrom, pos, ref, alt, qual, ao, dp, type, ref_len, alt_len)
    """
    with open(vcf_file) as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            info = fields[7]

            # Extract TYPE from INFO field
            variant_type = None
            ao = None
            dp = None

            for item in info.split(';'):
                if item.startswith('TYPE='):
                    variant_type = item.split('=')[1]
                elif item.startswith('AO='):
                    ao = item.split('=')[1]
                elif item.startswith('DP='):
                    dp = item.split('=')[1]

            # Filter for indels only
            if variant_type in ('ins', 'del'):
                yield {
                    'chrom': chrom,
                    'pos': int(pos),
                    'ref': ref,
                    'alt': alt,
                    'qual': float(qual),
                    'ao': int(ao) if ao else None,
                    'dp': int(dp) if dp else None,
                    'type': variant_type,
                    'ref_len': len(ref),
                    'alt_len': len(alt),
                    'sample': sample_name,
                }


def aggregate_indels(vcf_files, output_file):
    """
    Aggregate indels from multiple VCF files, grouping identical variants.

    Outputs two files:
    1. Aggregated TSV with quality metrics and sample lists
    2. Binary matrix for parsimony analysis

    Args:
        vcf_files: List of VCF file paths
        output_file: Output TSV file path (aggregated indels)
    """
    # Dictionary with key (chrom, pos, ref, alt) -> variant info + set of samples
    variants = defaultdict(lambda: {
        'type': None,
        'ref_len': None,
        'alt_len': None,
        'samples': set(),
        'quals': [],
    })

    all_samples = set()

    for vcf_file in vcf_files:
        sample_name = get_sample_name(vcf_file)
        if not sample_name:
            print(f"Warning: Could not extract sample name from {vcf_file}", file=sys.stderr)
            sample_name = Path(vcf_file).stem

        print(f"Processing {vcf_file} (sample: {sample_name})...", file=sys.stderr)
        try:
            for indel in parse_vcf_indels(vcf_file, sample_name):
                key = (indel['chrom'], indel['pos'], indel['ref'], indel['alt'])
                variants[key]['type'] = indel['type']
                variants[key]['ref_len'] = indel['ref_len']
                variants[key]['alt_len'] = indel['alt_len']
                variants[key]['samples'].add(indel['sample'])
                variants[key]['quals'].append(indel['qual'])
                all_samples.add(indel['sample'])
        except Exception as e:
            print(f"Error processing {vcf_file}: {e}", file=sys.stderr)
            continue

    # Write aggregated TSV output
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join([
            'chrom', 'pos', 'ref', 'alt', 'type',
            'ref_len', 'alt_len', 'num_samples', 'samples',
            'min_qual', 'max_qual', 'mean_qual'
        ]) + '\n')

        # Write variants sorted by position
        for (chrom, pos, ref, alt), info in sorted(variants.items(), key=lambda x: (x[0][0], x[0][1])):
            quals = info['quals']
            samples_str = ','.join(sorted(info['samples']))
            mean_qual = sum(quals) / len(quals) if quals else 0

            f.write('\t'.join([
                chrom,
                str(pos),
                ref,
                alt,
                info['type'],
                str(info['ref_len']),
                str(info['alt_len']),
                str(len(info['samples'])),
                samples_str,
                f"{min(quals):.2f}" if quals else 'NA',
                f"{max(quals):.2f}" if quals else 'NA',
                f"{mean_qual:.2f}" if quals else 'NA',
            ]) + '\n')

    # Write binary matrix for parsimony input
    binary_output = output_file.replace('.tsv', '_binary_matrix.tsv')

    # Sort samples and variants for reproducibility
    sorted_samples = sorted(all_samples)
    sorted_variants = sorted(variants.items(), key=lambda x: (x[0][0], x[0][1]))

    # Write binary matrix
    with open(binary_output, 'w') as f:
        # Write header with sample names
        f.write('indel_id\t' + '\t'.join(sorted_samples) + '\n')

        # Write each indel as a row
        for (chrom, pos, ref, alt), info in sorted_variants:
            indel_id = f"{chrom}:{pos}:{ref}:{alt}"
            row = [indel_id]

            for sample in sorted_samples:
                row.append('1' if sample in info['samples'] else '0')

            f.write('\t'.join(row) + '\n')

    total_variants = len(variants)
    total_observations = sum(len(info['samples']) for info in variants.values())
    print(f"Found {total_variants} unique indels across {len(all_samples)} samples", file=sys.stderr)
    print(f"Total observations: {total_observations}", file=sys.stderr)
    print(f"Aggregated indels: {output_file}", file=sys.stderr)
    print(f"Binary matrix: {binary_output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Filter indels from VCF files and aggregate quality scores'
    )
    parser.add_argument(
        'vcf_files',
        nargs='+',
        help='VCF file(s) to process'
    )
    parser.add_argument(
        '-o', '--output',
        default='indels_aggregated.tsv',
        help='Output TSV file (default: indels_aggregated.tsv)'
    )

    args = parser.parse_args()

    # Verify input files exist
    for vcf_file in args.vcf_files:
        if not Path(vcf_file).exists():
            print(f"Error: VCF file not found: {vcf_file}", file=sys.stderr)
            sys.exit(1)

    aggregate_indels(args.vcf_files, args.output)


if __name__ == '__main__':
    main()
