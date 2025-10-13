# r.anand@exeter.ac.uk
# Module for VCF file manipulation and extracting features

try:
    # Works when used inside a package
    from .fasta import translate_sequences, CODON_TABLE
except ImportError:
    # Works when run as a standalone script
    from modules.fasta import translate_sequences, CODON_TABLE


# IUPAC codes for heterozygous bases
IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
}


# --------------------------------------------------------------------
# header + genotype parsing helpers
# --------------------------------------------------------------------

def get_vcf_header(vcf_line, vcf_struct):
    """Parse the VCF header line (#CHROM...) and extract sample names."""
    bits = vcf_line.strip().split("\t")
    vcf_struct["next"] = 1
    vcf_struct["header"] = "Y"

    if vcf_line.startswith("#CHROM"):
        vcf_struct["isolate_names"] = {}
        for i in range(9, len(bits)):
            vcf_struct["isolate_names"][i - 9] = bits[i]
    return vcf_struct


def get_vcf_struct_and_determine_base_type(vcf_struct, gt_value):
    """Determine base1, base2, and base_type from genotype and REF/ALT alleles."""
    ref = vcf_struct.get("ref", "N")
    alt = vcf_struct.get("alt", "N")

    if ref == "N" or alt == "N" or gt_value in (".", "./.", ".|."):
        return "N", None, "ambiguous"

    bases = [ref] + alt.split(",")

    # Normalize GT (split phasing/slash)
    gt_parts = gt_value.replace("|", "/").split("/")

    # Single allele (haploid VCF)
    if len(gt_parts) == 1:
        allele_idx = gt_parts[0]
        if allele_idx == "0":
            return ref, None, "reference"
        else:
            allele_idx = int(allele_idx)
            if allele_idx < len(bases):
                allele = bases[allele_idx]
                if len(allele) == len(ref):
                    if len(ref) == 1:
                        return allele, None, "snp"
                    return allele, None, "snp_multi"
                elif len(allele) > len(ref):
                    return allele, None, "insertion"
                elif len(allele) < len(ref):
                    return allele, None, "deletion"
            return "N", None, "ambiguous"

    # Diploid or polyploid
    a1, a2 = gt_parts[0], gt_parts[1]

    # Homozygous reference
    if a1 == a2 == "0":
        return ref, None, "reference"

    # Homozygous alternate (1/1, 2/2, etc.)
    if a1 == a2 and a1 != "0":
        allele_idx = int(a1)
        if allele_idx < len(bases):
            allele = bases[allele_idx]
            if len(allele) == len(ref):
                if len(ref) == 1:
                    return allele, None, "snp"
                return allele, None, "snp_multi"
            elif len(allele) > len(ref):
                return allele, None, "insertion"
            elif len(allele) < len(ref):
                return allele, None, "deletion"
        return "N", None, "ambiguous"

    # Heterozygous (0/1, 1/0, 1/2, etc.)
    if a1 != a2:
        try:
            base1 = bases[int(a1)]
        except (IndexError, ValueError):
            base1 = "N"
        try:
            base2 = bases[int(a2)]
        except (IndexError, ValueError):
            base2 = "N"

        # Determine indel vs SNP heterozygotes
        if len(base1) != len(base2):
            if len(base1) > len(base2):
                base_type = "het_insertion"
            elif len(base1) < len(base2):
                base_type = "het_deletion"
            else:
                base_type = "heterozygous"
        else:
            base_type = "heterozygous"

        return base1, base2, base_type

    return "N", None, "ambiguous"


def parse_vcf_line(vcf_line, vcf_struct=None):
    """Parse a single VCF line into a dict of fields and sample info."""
    bits = vcf_line.strip().split("\t")
    vcf_info = {} if vcf_struct is None else vcf_struct

    if vcf_line.startswith("#"):
        return get_vcf_header(vcf_line, vcf_info)

    vcf_info["next"] = 0
    vcf_info["header"] = "N"

    if len(bits) < 9:
        print(f"bad vcf line with <9 columns: {vcf_line}")
        vcf_info["next"] = 1
        return vcf_info

    if len(bits) > 9:
        for i in range(9, len(bits)):
            vcf_info[f"sample_info_{i-9}"] = bits[i]

    vcf_info.update({
        "supercontig": bits[0],
        "position": bits[1],
        "id": bits[2],
        "ref": bits[3],
        "alt": bits[4],
        "cons_qual": bits[5],
        "filter": bits[6],
        "info": bits[7],
        "format": bits[8],
    })

    format_parts = vcf_info["format"].split(":")
    vcf_info["format_parts_available"] = {part: True for part in format_parts}
    vcf_info["number_of_samples"] = len(bits) - 9

    # per-sample data
    for i in range(9, len(bits)):
        isolate_number = i - 9
        sample_data = bits[i]
        sample_parts = sample_data.split(":")
        if len(sample_parts) < len(format_parts):
            raise ValueError(f"{vcf_line} sample info has fewer parts than FORMAT")
        sample_parts = sample_parts[:len(format_parts)]
        sample_dict = dict(zip(format_parts, sample_parts))
        for key, val in sample_dict.items():
            vcf_info[f"{key}{isolate_number}"] = val
        gt_value = sample_dict.get("GT", ".")
        base1, base2, base_type = get_vcf_struct_and_determine_base_type(vcf_info, gt_value)
        vcf_info[f"{isolate_number}_base1"] = base1
        vcf_info[f"{isolate_number}_base2"] = base2
        vcf_info[f"base_type{isolate_number}"] = base_type
        if base_type == "heterozygous" and base1 and base2:
            code = IUPAC_CODES.get(frozenset([base1, base2]), "N")
            vcf_info[f"amb_char{isolate_number}"] = code
    return vcf_info

# --------------------------------------------------------------------
# filter logic
# --------------------------------------------------------------------
def filter_vcf(vcf_file, min_depth=20, remove_ambiguous=False):
    """
    Filter a VCF file by minimum read depth and optionally remove ambiguous sites.
    Writes filtered lines to stdout.
    """
    with open(vcf_file, "r") as fh:
        for line in fh:
            if line.startswith("##") or line.startswith("#CHROM"):
                print(line.strip())
                continue

            bits = line.strip().split("\t")
            if len(bits) < 10:
                continue

            # parse FORMAT column
            fmt = bits[8].split(":")
            if "DP" not in fmt:
                print(line.strip())  # keep if no depth info
                continue

            dp_index = fmt.index("DP")

            # check all sample depths
            sample_ok = True
            for sample_col in bits[9:]:
                sample_parts = sample_col.split(":")
                if len(sample_parts) > dp_index:
                    try:
                        dp_val = int(sample_parts[dp_index])
                        if dp_val < min_depth:
                            sample_ok = False
                            break
                    except ValueError:
                        sample_ok = False
                        break

            if not sample_ok:
                continue

            if remove_ambiguous and ("N" in bits[3] or "N" in bits[4]):
                continue

            print(line.strip())
# --------------------------------------------------------------------
# remove ref sites logic 
# --------------------------------------------------------------------
def remove_ref_sites(vcf_file):
    """
    Remove all non-snp ref sites from a vcf file
    Output to stdout
    """
    vcf_struct = {}
    with open(vcf_file, "r") as fh:
        for line in fh:
            if line.startswith("##") or line.startswith("#CHROM"):
                print(line.strip())
                continue
            vcf_info = parse_vcf_line(line, vcf_struct)
            if vcf_info.get("header") == "Y":
                continue
            num_samples = vcf_info["number_of_samples"]
            variant_found = False

            for i in range(num_samples):
                base_type = vcf_info.get("base_type{i}", "unknown")

                if base_type not in ("reference", "ambiguous", "unknown"):
                    variant_found = True
                    break
            if variant_found:
                print(line.strip())

# --------------------------------------------------------------------
# stats logic 
# --------------------------------------------------------------------
def vcf_stats(vcf_file):
    """
    Print detailed summary statistics for a VCF file.
    Includes separate counts for heterozygous vs homozygous variant types.
    """
    total_variants = 0
    snps = 0
    het_snps = 0
    insertions = 0
    het_insertions = 0
    deletions = 0
    het_deletions = 0
    ambiguous = 0
    unknown = 0

    with open(vcf_file, "r") as fh:
        vcf_struct = {}
        for line in fh:
            if line.startswith("#"):
                continue

            vcf_info = parse_vcf_line(line, vcf_struct)
            if vcf_info.get("header") == "Y":
                continue

            total_variants += 1

            for i in range(vcf_info["number_of_samples"]):
                base_type = vcf_info.get(f"base_type{i}", "unknown")

                if base_type.startswith("snp"):
                    snps += 1
                elif base_type == "heterozygous":
                    het_snps += 1
                elif base_type == "insertion":
                    insertions += 1
                elif base_type == "het_insertion":
                    het_insertions += 1
                elif base_type == "deletion":
                    deletions += 1
                elif base_type == "het_deletion":
                    het_deletions += 1
                elif base_type == "ambiguous":
                    ambiguous += 1
                else:
                    unknown += 1

    # --- Print summary ---
    print("\nVCF Summary Statistics")
    print("----------------------")
    print(f"Total variant sites:        {total_variants}")
    print()
    print(f"Homozygous SNPs:            {snps}")
    print(f"Heterozygous SNPs:          {het_snps}")
    print(f"Homozygous insertions:      {insertions}")
    print(f"Heterozygous insertions:    {het_insertions}")
    print(f"Homozygous deletions:       {deletions}")
    print(f"Heterozygous deletions:     {het_deletions}")
    print()
    print(f"Ambiguous sites:            {ambiguous}")
    print(f"Unknown / other:            {unknown}")
    print("----------------------")
    print(f"Total variants classified:  {snps + het_snps + insertions + het_insertions + deletions + het_deletions + ambiguous + unknown}")

# --------------------------------------------------------------------
# annotation logic
# --------------------------------------------------------------------

def annotate_vcf_file(vcf_info, fasta_sequences, gff_records, summary_filename="summary_annotation.tab"):
    """
    Annotate variants from VCF using FASTA and GFF info.
    Prints annotated variants to stdout and appends CDS/gene variants to summary_annotation.tab.
    """
    chrom = vcf_info["supercontig"]
    pos = int(vcf_info["position"])
    ref = vcf_info["ref"]
    alt = vcf_info["alt"]

    id_to_record = {rec["attributes"].get("ID"): rec for rec in gff_records if "ID" in rec["attributes"]}
    per_sample_annotations = {}

    for i in range(vcf_info["number_of_samples"]):
        base_type = vcf_info.get(f"base_type{i}", "unknown")
        annotation = {
            "feature": "intergenic",
            "gene_id": None,
            "gene_alias": None,
            "gene_name": None,
            "effect": "N/A",
            "base_type": base_type,
        }

        for rec in gff_records:
            if rec["seqid"] != chrom or not (rec["start"] <= pos <= rec["end"]):
                continue

            annotation["feature"] = rec["type"]
            annotation["gene_id"] = rec["attributes"].get("ID")
            annotation["gene_alias"] = rec["attributes"].get("Alias")
            annotation["gene_name"] = rec["attributes"].get("Name")

            # climb to gene
            parent_id = rec["attributes"].get("Parent")
            while parent_id and parent_id in id_to_record:
                parent_rec = id_to_record[parent_id]
                if parent_rec["type"] == "gene":
                    annotation["gene_id"] = parent_rec["attributes"].get("ID")
                    annotation["gene_alias"] = parent_rec["attributes"].get("Alias")
                    annotation["gene_name"] = parent_rec["attributes"].get("Name")
                    break
                parent_id = parent_rec["attributes"].get("Parent")

            # CDS logic (uppercase)
            if rec["type"] == "CDS":
                seq = fasta_sequences[chrom]
                codon_index = (pos - rec["start"]) // 3
                codon_start = rec["start"] + codon_index * 3
                ref_codon = seq[codon_start - 1 : codon_start + 2]
                codon_pos = pos - codon_start

                if base_type.startswith("snp"):
                    if len(ref) == 1 and len(alt) == 1:
                        alt_codon = ref_codon[:codon_pos] + alt + ref_codon[codon_pos + 1 :]
                        ref_aa = translate_sequences(ref_codon)
                        alt_aa = translate_sequences(alt_codon)
                        if ref_aa != "*" and alt_aa == "*":
                            annotation["effect"] = "stop_gained"
                        elif ref_aa == "*" and alt_aa != "*":
                            annotation["effect"] = "stop_lost"
                        elif ref_aa == alt_aa:
                            annotation["effect"] = "synonymous"
                        else:
                            annotation["effect"] = "nonsynonymous"

                elif base_type == "heterozygous":
                    base1 = vcf_info.get(f"{i}_base1")
                    base2 = vcf_info.get(f"{i}_base2")
                    ref_aa = translate_sequences(ref_codon)
                    codon1 = ref_codon[:codon_pos] + base1 + ref_codon[codon_pos + 1 :]
                    codon2 = ref_codon[:codon_pos] + base2 + ref_codon[codon_pos + 1 :]
                    aa1 = translate_sequences(codon1)
                    aa2 = translate_sequences(codon2)
                    if "*" in (aa1, aa2) and "*" not in (ref_aa,):
                        annotation["effect"] = "stop_gained_het"
                    elif ref_aa == "*" and ("*" not in (aa1, aa2)):
                        annotation["effect"] = "stop_lost_het"
                    elif aa1 == ref_aa and aa2 == ref_aa:
                        annotation["effect"] = "synonymous_het"
                    elif aa1 != ref_aa and aa2 != ref_aa:
                        annotation["effect"] = "nonsynonymous_het"
                    else:
                        annotation["effect"] = "mixed_het"

                elif base_type in ("insertion", "deletion", "het_insertion", "het_deletion"):
                    indel_len = abs(len(alt) - len(ref))
                    if indel_len % 3 == 0:
                        annotation["effect"] = "inframe_indel"
                    else:
                        window = seq[pos - 1 : pos + 30]
                        shifted_seq = window[len(alt) - len(ref) :]
                        aa_seq = translate_sequences(shifted_seq)
                        annotation["effect"] = "frameshift_stop_gained" if "*" in aa_seq[:10] else "frameshift_indel"

            elif rec["type"] in ("five_prime_utr", "3utr", "three_prime_utr"):
                annotation["effect"] = "utr_variant"
            elif rec["type"] == "intron":
                annotation["effect"] = "intron_variant"
            elif rec["type"] == "exon":
                annotation["effect"] = "exonic_non_coding"
            elif rec["type"] == "gene":
                annotation["effect"] = "genic_non_coding"
            else:
                annotation["effect"] = f"noncoding_{rec['type']}"
            break

        per_sample_annotations[i] = annotation

    vcf_info["annotations"] = per_sample_annotations

    # print results
    for i, ann in per_sample_annotations.items():
        print(
            f"{chrom}\t{pos}\t{ref}->{alt}\t"
            f"sample{i}\t{ann['feature']}\t"
            f"gene_id={ann['gene_id']}\tgene_alias={ann['gene_alias']}\tgene_name={ann['gene_name']}\t"
            f"{ann['effect']}\t{ann['base_type']}"
        )

    # write summary
    with open(summary_filename, "a") as out:
        for i, ann in per_sample_annotations.items():
            if ann["feature"] in ("CDS", "gene"):
                out.write(
                    f"sample{i}\t{ann['gene_id'] or 'NA'}\t"
                    f"{ann['gene_name'] or 'NA'}\t{ann['gene_alias'] or 'NA'}\t"
                    f"{ann['effect']}\n"
                )
    return vcf_info


def vcf_file_to_fasta(vcf_file, fasta_sequences):
    """
    Python function to convert VCF file to fasta sequence

    """

    sample_sequences = {}
    sample_names = []

    vcf_struct = {}

    with open (vcf_file, "r") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            vcf_info = parse_vcf_line(vcf_file, vcf_struct)

            if "isolate_names" in vcf_info and not sample_names:
                sample_names = [vcf_info["isolate_names"][i] for i in sorted(vcf_info["isolate_names"].keys())]
                for name in sample_names:
                    sample_sequences[name] = fasta.sequences.copy()

            if vcf_info.get("header") == "Y":
                continue

            chrom = vcf_info["supercontig"]
            pos = int(vcf_info["position"])
            ref = vcf_info["ref"]
            alt = vcf_info["alt"]

            for i, sample_name in enumerate(sample_names):
                seq_dict = sample_sequences[sample_name]
                if chrom not in seq_dict:
                    continue

                seq_list = list(seq_dict)

                gt = vcf_info.get(f"GT{i}", 0)
                base_type = vcf_info.get(f"base_type{i}", "reference")

                if base_type == "reference" or gt == "0":
                    continue

                if base_type.startswith("snp"):
                    seq_list[pos - 1] = alt[0]

                elif base_type == "heterozygous":
                    base1 = vcf_info.get(f"{i}_base1")
                    base2 = vcf_info.get(f"{i}_base2")
                    code = IUPAC_CODES.gets(frozenset([base1, base2]), "N")
                    seq_list[pos - 1] = code

                elif base_type in ("insertion", "het_insertion"):
                    seq_list[pos - 1 : pos - 1 + len(ref)] = list(alt)

                elif base_type in ("deletion", "het_deletion"):
                    seq_list[pos - 1 : pos - 1 + len(ref)] = []

                seq_dict[chrom] = "".join(seq_list)
                sample_sequences[sample_name] = seq_dict

    for sample_name, seqs in sample_sequences.items():
        outname = "{sample_name}.fasta"
        with open(outname, "w") as out:
            for chrom, seq in seq.items():
                out.write(f">{chrom}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")
        print(f"Wrote {outname}")
    print("All sample FASTA sequences generated successfully")


def vcf_name_location_merge(name_location_file, output_file="merged.vcf"):
    """
    Merge multiple single-sample VCF files into one multi-sample VCF.

    Input file format (tab-separated):
        Sample1   path/to/Sample1.vcf
        Sample2   path/to/Sample2.vcf
        Sample3   path/to/Sample3.vcf

    Output:
        merged.vcf  – combined VCF with all samples as separate columns.

    Notes:
        - Assumes all VCFs have the same reference and coordinate system.
        - Lines are matched by CHROM and POS.
        - INFO and FILTER fields are taken from the first VCF if duplicated.
    """

    import collections

    sample_names = []
    vcf_files = []

    # --- Read sample name + VCF mapping ---
    with open(name_location_file, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            name, path = line.strip().split("\t")
            sample_names.append(name)
            vcf_files.append(path)

    if not vcf_files:
        print("[ERROR] No VCFs provided for merging.")
        return

    print(f"[INFO] Merging {len(vcf_files)} VCF files into {output_file}")

    # --- Read headers and variant lines ---
    merged_header = []
    variant_dict = collections.OrderedDict()  # key: (CHROM, POS, REF, ALT)

    for idx, vcf_path in enumerate(vcf_files):
        with open(vcf_path, "r") as fh:
            for line in fh:
                if line.startswith("##"):
                    if idx == 0:  # keep meta-info only once
                        merged_header.append(line.strip())
                    continue

                if line.startswith("#CHROM"):
                    if idx == 0:
                        header_parts = line.strip().split("\t")
                        merged_header_line = "\t".join(header_parts[:9] + sample_names)
                        merged_header.append(merged_header_line)
                    continue

                bits = line.strip().split("\t")
                if len(bits) < 10:
                    continue

                chrom, pos, _id, ref, alt, qual, flt, info, fmt = bits[:9]
                sample_data = bits[9]

                key = (chrom, pos, ref, alt)
                if key not in variant_dict:
                    # Initialize with placeholders for all samples
                    variant_dict[key] = {
                        "chrom": chrom,
                        "pos": pos,
                        "id": _id,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual,
                        "filter": flt,
                        "info": info,
                        "format": fmt,
                        "samples": ["./."] * len(vcf_files),
                    }

                variant_dict[key]["samples"][idx] = sample_data

    # --- Write merged VCF ---
    with open(output_file, "w") as out:
        for hline in merged_header:
            out.write(hline + "\n")
        for key, record in variant_dict.items():
            out.write(
                "\t".join([
                    record["chrom"],
                    record["pos"],
                    record["id"],
                    record["ref"],
                    record["alt"],
                    record["qual"],
                    record["filter"],
                    record["info"],
                    record["format"],
                ] + record["samples"]) + "\n"
            )

    print(f"[INFO] Merged VCF written to: {output_file}")
