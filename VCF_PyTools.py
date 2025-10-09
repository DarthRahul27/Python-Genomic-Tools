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
        description="Tools to parse, filter, annotate, and merge VCF files"
    )

    parser.add_argument("-v", "--vcf", help="Input VCF file (required for most operations)")
    parser.add_argument("--filter", choices=("y", "n"), default="n",
                        help="Filter VCF file and print output to stdout")
    parser.add_argument("--depth_filter", type=int, default=20,
                        help="If --filter=y, minimum depth filter (default: 20)")
    parser.add_argument("--remove_ambiguous", choices=("y", "n"), default="n",
                        help="If --filter=y, remove any sites that are ambiguous")
    parser.add_argument("--stats", choices=("y", "n"), default="n",
                        help="Calculate and print summary statistics from VCF file")
    parser.add_argument("--annotate", choices=("y", "n"), default="n",
                        help="Annotate variants using GFF3 and FASTA [requires --gff3 and --fasta]")
    parser.add_argument("--convert_to_fasta", choices=("y", "n"), default="n",
                        help="Convert VCF file to FASTA per-sample consensus [requires --fasta reference]")
    parser.add_argument("--merge", choices=("y", "n"), default="n",
                        help="Merge multiple VCFs listed in a name-location file into one multi-sample VCF")
    parser.add_argument("--name_location_file",
                        help="Tab-delimited file with sample names and VCF paths (for use with --merge=y)")
    parser.add_argument("--output_merged_vcf", default="merged.vcf",
                        help="Output filename for merged VCF (default: merged.vcf)")
    parser.add_argument("--gff3", help="GFF3 file (required if --annotate=y)")
    parser.add_argument("--fasta", help="FASTA file (required if --annotate=y or --convert_to_fasta=y)")

    args = parser.parse_args()

    # ----------------------------------------------------------------
    # Validate input combinations
    # ----------------------------------------------------------------
    if args.annotate == "y" and (not args.gff3 or not args.fasta):
        parser.error("--annotate=y requires both --gff3 and --fasta files")

    if args.convert_to_fasta == "y" and not args.fasta:
        parser.error("--convert_to_fasta=y requires --fasta reference")

    if args.merge == "y" and not args.name_location_file:
        parser.error("--merge=y requires --name_location_file with sample\tvcf paths")

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
        if not args.vcf:
            parser.error("--vcf required when using --filter=y")
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
        if not args.vcf:
            parser.error("--vcf required when using --convert_to_fasta=y")
        print("\n[INFO] Converting VCF to FASTA sequences...")
        vcf.vcf_file_to_fasta(args.vcf, fasta_sequences)
        print("[INFO] FASTA conversion complete.")

    # ----------------------------------------------------------------
    # Option 4: VCF Statistics
    # ----------------------------------------------------------------
    elif args.stats == "y":
        if not args.vcf:
            parser.error("--vcf required when using --stats=y")
        print("\n[INFO] Calculating VCF summary statistics...")
        vcf.vcf_stats(args.vcf)
        print("[INFO] VCF statistics summary complete.")

    # ----------------------------------------------------------------
    # Option 5: Merge multiple VCFs
    # ----------------------------------------------------------------
    elif args.merge == "y":
        print(f"\n[INFO] Merging VCFs listed in {args.name_location_file}...")
        vcf.vcf_name_location_merge(args.name_location_file, args.output_merged_vcf)
        print(f"[INFO] Merged VCF written to {args.output_merged_vcf}")

    else:
        print("\n[INFO] No operation selected. "
              "Use --annotate=y, --filter=y, --convert_to_fasta=y, --stats=y, or --merge=y.")


if __name__ == "__main__":
    main()

