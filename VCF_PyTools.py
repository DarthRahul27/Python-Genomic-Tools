#!/usr/bin/env python3
## Python script for handling and manipulating VCF files ##

import os
import sys

# --- Allow local imports whether run as script or module ---
sys.path.append(os.path.dirname(__file__))

import argparse
from modules import vcf, fasta, gff


def main():
    script_name = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(
        prog=script_name,
        description="Tools to parse and manipulate VCF files"
    )

    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file [REQUIRED]")
    parser.add_argument("--filter", choices=("y", "n"), default="n",
                        help="Option to filter vcf file and output into stdout")
    parser.add_argument("--depth_filter", type=int, default=20,
                        help="If --filter=y, minimum depth filter (default: 20)")
    parser.add_argument("--remove_ambiguous", choices=("y", "n"), default="n",
                        help="If --filter=y, remove any sites that are ambiguous")
    parser.add_argument("--stats", choices=("y", "n"), default="n",
                        help="Option to print basic stats from VCF file")
    parser.add_argument("--annotate", choices=("y", "n"), default="n",
                        help="Option to annotate variants in VCF file [requires --gff3 and --fasta]")
    parser.add_argument("--convert_to_fasta", choices=("y", "n"), default="n",
                        help="Convert VCF file to FASTA per-sample consensus [requires --fasta reference]")
    parser.add_argument("--gff3", help="GFF3 file (required if --annotate=y)")
    parser.add_argument("--fasta", help="FASTA file (required if --annotate=y or --convert_to_fasta=y)")

    args = parser.parse_args()

    # ----------------------------------------------------------------
    # Validation
    # ----------------------------------------------------------------
    if args.annotate == "y" and (not args.gff3 or not args.fasta):
        parser.error("--annotate=y requires both --gff3 and --fasta files")

    if args.convert_to_fasta == "y" and not args.fasta:
        parser.error("--convert_to_fasta=y requires --fasta reference")

    # ----------------------------------------------------------------
    # Load reference FASTA if needed
    # ----------------------------------------------------------------
    fasta_sequences = None
    if args.fasta:
        fasta_sequences = fasta.parse_fasta(args.fasta)

    # ----------------------------------------------------------------
    # Option 1: Filter VCF
    # ----------------------------------------------------------------
    if args.filter == "y":
        print(f"\n[INFO] Filtering {args.vcf} with depth >= {args.depth_filter}, "
              f"remove ambiguous = {args.remove_ambiguous}")
        vcf.filter_vcf(args.vcf, args.depth_filter, args.remove_ambiguous == "y")
        print("[INFO] Filtering complete. Filtered VCF printed to stdout.")

    # ----------------------------------------------------------------
    # Option 2: Annotation
    # ----------------------------------------------------------------
    elif args.annotate == "y":
        print("\n[INFO] Annotating VCF file using GFF3 and FASTA reference...")
        gff_records = gff.parse_gff(args.gff3)
        with open(args.vcf, "r") as fh:
            vcf_struct = {}
            for line in fh:
                if line.startswith("##"):
                    continue
                vcf_info = vcf.parse_vcf_line(line, vcf_struct)
                if vcf_info.get("header") == "Y":
                    continue
                vcf.annotate_vcf_file(vcf_info, fasta_sequences, gff_records)
        print("[INFO] Annotation complete. Summary written to summary_annotation.tab")

    # ----------------------------------------------------------------
    # Option 3: Convert to FASTA
    # ----------------------------------------------------------------
    elif args.convert_to_fasta == "y":
        print("\n[INFO] Converting VCF to FASTA sequences...")
        vcf.vcf_file_to_fasta(args.vcf, fasta_sequences)
        print("[INFO] FASTA conversion complete.")

    # ----------------------------------------------------------------
    # Option 4: VCF Statistics
    # ----------------------------------------------------------------
    elif args.stats == "y":
        print("\n[INFO] Calculating VCF summary statistics...")
        vcf.vcf_stats(args.vcf)
        print("[INFO] VCF statistics summary complete.")

    else:
        print("\n[INFO] No operation selected. Use --annotate=y, --filter=y, or --convert_to_fasta=y.")


if __name__ == "__main__":
    main()

