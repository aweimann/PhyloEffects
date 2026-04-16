#!/usr/bin/env python3
"""Port of the legacy Bash workflow in pipeline.txt.

This script mirrors the active steps from the Bash pipeline while keeping
argument compatibility with the original positional interface.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

from parse_gubbins import extract as parse_gubbins_extract
from phyloeffects import main as run_phyloeffects
from rescale_tree import prune_tree
from summarise_embl_output import extract as summarise_embl_extract


def run_command(cmd: list[str]) -> None:
    """Run a command and fail fast on non-zero exit status."""
    subprocess.run(cmd, check=True)


def write_position_mapping(vcf_path: Path, output_path: Path) -> None:
    """Create the index-to-genome-position mapping used by PhyloEffects.

    Replicates:
    cat <(echo $'index\tposition') <(cut -f2 vcf | tail -n+5 | awk '{print NR"\t"$0}')
    """
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with vcf_path.open("r", newline="") as handle_in, output_path.open("w", newline="") as handle_out:
        writer = csv.writer(handle_out, delimiter="\t")
        writer.writerow(["index", "position"])

        idx = 1
        for line in handle_in:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            writer.writerow([idx, fields[1]])
            idx += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the PhyloEffects workflow (Python port of pipeline.txt). "
            "Positional arguments intentionally match the legacy Bash script."
        )
    )
    parser.add_argument("cluster", help="Name of cluster")
    parser.add_argument("snp_sites", help="Path prefix for snp-sites outputs")
    parser.add_argument("gubbins", help="Path prefix for Gubbins outputs")
    parser.add_argument(
        "outgroup",
        help="Outgroup for tree pruning",
    )
    parser.add_argument("parsimony", help="Parsimony output path (currently unused)")
    parser.add_argument("pseudoalignment", help="Pseudoalignment input path")
    parser.add_argument("phyloeffects", help="Output path for phylogenetic effect outputs")
    parser.add_argument("gff", help="Genome annotation GFF input")
    parser.add_argument("genome", help="Reference genome FASTA input")
    parser.add_argument(
        "--run-snp-sites",
        action="store_true",
        help="Also run snp-sites VCF and FASTA generation (commented out in legacy Bash script)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cluster = args.cluster
    # gubbins output dir
    gubbins_dir = Path(args.gubbins)

    # snp-sites output dir
    snp_sites_dir = Path(args.snp_sites)
    phyloeffects_dir = Path(args.phyloeffects)
    # check if these directories exist, if not create them
    for directory in [snp_sites_dir, phyloeffects_dir]:
        directory.mkdir(parents=True, exist_ok=True)


    vcf_file = snp_sites_dir / f"{cluster}.vcf"
    aln_file = snp_sites_dir / f"{cluster}.fasta"
    pos_map_file = snp_sites_dir / f"{cluster}_pos_mapping.txt"

    if args.run_snp_sites:
        print("snp-sites")
        run_command(["snp-sites", "-v", "-o", str(vcf_file), str(args.pseudoalignment)])
        run_command(["snp-sites", "-o", str(aln_file), str(args.pseudoalignment)])

    write_position_mapping(vcf_file, pos_map_file)

    print("parse gubbins output")
    parse_gubbins_extract(
        str(gubbins_dir / f"{cluster}.recombination_predictions.gff"),
        str(snp_sites_dir / f"{cluster}.recombination_prediction.txt"),
    )
    summarise_embl_extract(
        str(gubbins_dir / f"{cluster}.branch_base_reconstruction.embl"),
        str(snp_sites_dir / f"{cluster}.recombination_pos.txt"),
    )

    print("rescale tree")
    phyloeffects_dir.mkdir(parents=True, exist_ok=True)
    prune_tree(
        input_tree=str(gubbins_dir / f"{cluster}.node_labelled.final_tree.tre"),
        output_tree=str(phyloeffects_dir / f"{cluster}_rescaled.nwk"),
        alignment=str(aln_file),
        alignment_out=str(phyloeffects_dir / f"{cluster}_aln.fasta"),
        outgroup=args.outgroup,
        midpoint=False,
        rescale=True,
        relabel=False,
    )

    print("PhyloEffects")
    (phyloeffects_dir / cluster).mkdir(parents=True, exist_ok=True)
    run_phyloeffects(
        [
            "-a",
            str(phyloeffects_dir / f"{cluster}_aln.fasta"),
            "-t",
            str(phyloeffects_dir / f"{cluster}_rescaled.nwk"),
            "-r",
            str(args.genome),
            "-o",
            str(phyloeffects_dir / cluster),
            "-c",
            str(pos_map_file),
            "-g",
            str(args.gff),
        ]
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())