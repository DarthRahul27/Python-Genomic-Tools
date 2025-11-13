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
        description="Tools to parse, filter, annotate, merge VCFs, and compute LOH/window stats"
    )

    parser.add_argument("-v", "--vcf", help="Input VCF file (required for most operations)")
    parser.add_argument("--filter", choices=("y", "n"), default="n",
                        help="Filter VCF file and print output to stdout")
    parser.add_argument("--depth_filter", type=int, default=20,
                        help="If --filter=y, minimum depth filter (default: 20)")
    parser.add_argument("--remove_ambiguous", choices=("y", "n"), default="n",
                        help="If --filter=y, remove ambiguous sites")
    parser.add_argument("--remove_ref_sites", choices=("y", "n"), default="n",
                        help="Remove reference-only sites from VCF")
    parser.add_argument("--stats", choices=("y", "n"), default="n",
                        help="Print summary VCF statistics")
    parser.add_argument("--annotate", choices=("y", "n"), default="n",
                        help="Annotate variants using GFF3 + FASTA (--gff3 & --fasta required)")
    parser.add_argument("--convert_to_fasta", choices=("y", "n"), default="n",
                        help="Convert VCF variants to per-sample FASTA consensus")
    parser.add_argument("--merge", choices=("y", "n"), default="n",
                        help="Merge multiple VCFs from name-location file")
    parser.add_argument("--name_location_file",
                        help="Tab-delimited file with sample name + VCF path")
    parser.add_argument("--output_merged_vcf", default="merged.vcf",
                        help="Output filename for merged VCF")
    parser.add_argument("--make_10kb_windows", choices=("y", "n"), default="n",
                        help="Generate 10 kb genomic windows with variant + LOH statistics")
    parser.add_argument("--gff3", help="GFF3 file (required if --annotate=y)")
    parser.add_argument("--fasta", help="FASTA file (required for annotate/convert/windows)")

    args = parser.parse_args()

    # ----------------------------------------------------------------
    # Validate inputs
    # ----------------------------------------------------------------
    if args.annotate == "y" and (not args.gff3 or not args.fasta):
        parser.error("--annotate=y requires both --gff3 and --fasta")

    if args.convert_to_fasta == "y" and not args.fasta:
        parser.error("--convert_to_fasta=y requires --fasta reference")

    if args.merge == "y" and not args.name_location_file:
        parser.error("--merge=y requires --name_location_file")

    if not args.vcf and any(
        opt == "y" for opt in [
            args.filter, args.stats, args.annotate,
            args.convert_to_fasta, args.remove_ref_sites,
            args.make_10kb_windows
        ]
    ):
        parser.error("--vcf is required for this operation")

    # ----------------------------------------------------------------
    # Load FASTA if needed
    # ----------------------------------------------------------------
    fasta_sequences = None
    if args.fasta:
        fasta_sequences = fasta.parse_fasta(args.fasta)

    # ----------------------------------------------------------------
    # Option 1: Filter VCF
    # ----------------------------------------------------------------
    if args.filter == "y":
        print(f"[INFO] Filtering {args.vcf}...")
        vcf.filter_vcf(args.vcf, args.depth_filter, args.remove_ambiguous == "y")
        return

    # ----------------------------------------------------------------
    # Option 2: Remove reference-only sites
    # ----------------------------------------------------------------
    if args.remove_ref_sites == "y":
        print(f"[INFO] Removing reference-only sites from {args.vcf}...")
        vcf.remove_ref_sites(args.vcf)
        return

    # ----------------------------------------------------------------
    # Option 3: Annotation
    # ----------------------------------------------------------------
    if args.annotate == "y":
        print("[INFO] Annotating variants...")
        gff_records = gff.parse_gff(args.gff3)
        vcf_struct = {}
        with open(args.vcf) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                vcf_info = vcf.parse_vcf_line(line, vcf_struct)
                if vcf_info.get("header") == "Y":
                    continue
                vcf.annotate_vcf_file(vcf_info, fasta_sequences, gff_records)
        print("[INFO] Annotation complete")
        return

    # ----------------------------------------------------------------
    # Option 4: Convert to FASTA
    # ----------------------------------------------------------------
    if args.convert_to_fasta == "y":
        print("[INFO] Converting VCF → FASTA...")
        vcf.vcf_file_to_fasta(args.vcf, fasta_sequences)
        return

    # ----------------------------------------------------------------
    # Option 5: Summary statistics
    # ----------------------------------------------------------------
    if args.stats == "y":
        print("[INFO] Computing VCF statistics...")
        vcf.vcf_stats(args.vcf)
        return

    # ----------------------------------------------------------------
    # Option 6: Merge VCFs
    # ----------------------------------------------------------------
    if args.merge == "y":
        print("[INFO] Merging VCFs...")
        vcf.vcf_name_location_merge(args.name_location_file, args.output_merged_vcf)
        return

    # ----------------------------------------------------------------
    # Option 7: 10 kb WINDOWS + LOH
    # ----------------------------------------------------------------
    if args.make_10kb_windows == "y":
        if not args.fasta:
            parser.error("--make_10kb_windows=y requires --fasta")

        print("\n[INFO] Generating 10 kb windowed variant + LOH statistics...")

        stats = vcf.windowed_loh_stats(
            vcf_file=args.vcf,
            fasta_file=args.fasta,
            window_size=10_000,
            het_gap=100,
        )

        vcf.write_windowed_loh_tsv(stats, output_prefix="windows10kb")
        print("[INFO] Window files written: windows10kb.SAMPLE.tsv")
        return

    # ----------------------------------------------------------------
    # Default: Nothing chosen
    # ----------------------------------------------------------------
    print("[INFO] No operation selected.")


if __name__ == "__main__":
    main()

