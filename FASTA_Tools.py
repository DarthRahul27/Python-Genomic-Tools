#!/usr/bin/env python
### Python script that allows us to parse FASTA files, get stats and other related functions ###


import argparse
import sys
from modules import fasta

def main():
	
	script_name = sys.argv[0].split("/")[-1]
	parser = argparse.ArgumentParser(
		prog=script_name,
		description="FASTA Pasrer and Statitics Tools"
		)
	parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
	parser.add_argument("-n", "--names", help="Tab delimted files with old chromosome names tab new chromosome names")
	parser.add_argument("-o", "--output_file_name", help="Output file if -n with renamed chromosome")
	parser.add_argument(
		"-s", "--stats",
		choices=("y", "n"),
		default="y",
		help="Show FASTA stats in table [y/n] default y"
		)

	args = parser.parse_args()

	records = fasta.parse_fasta(args.fasta)

	if args.names:
		mapping = {}
		with open(args.names) as fh:
			for line in fh:
				old, new = line.strip().split("\t")
				mapping[old] = new
		records = {mapping.get(h, h): seq for h, seq in records.items()}

	if args.output_file_name:
		with open(args.output_file_name, "w") as out:
			for h, seq in records.items():
				out.write(f">{h}\n")

				for i in range(0, len(seq), 60):
					out.write(seq[i:i+60] + "\n")
		print(f"FASTA written to {args.output}")

	if args.stats == "y":
		fasta.print_fasta_stats(records)


if __name__ == "__main__":
	main()