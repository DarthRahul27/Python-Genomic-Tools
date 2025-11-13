# r.anand@exeter.ac.uk
# Module for VCF file manipulation and extracting features

try:
    # When inside package
    from .fasta import parse_fasta, translate_sequences, CODON_TABLE
except ImportError:
    # When standalone
    from modules.fasta import parse_fasta, translate_sequences, CODON_TABLE


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
# HEADER & GENOTYPE PARSING HELPERS
# --------------------------------------------------------------------

def get_vcf_header(vcf_line, vcf_struct):
    """Parse #CHROM line and extract sample names."""
    bits = vcf_line.strip().split("\t")
    vcf_struct["next"] = 1
    vcf_struct["header"] = "Y"

    if vcf_line.startswith("#CHROM"):
        vcf_struct["isolate_names"] = {}
        for i in range(9, len(bits)):
            vcf_struct["isolate_names"][i - 9] = bits[i]
    return vcf_struct


def get_vcf_struct_and_determine_base_type(vcf_struct, gt_value):
    """Determine genotype → base type."""
    ref = vcf_struct.get("ref", "N")
    alt = vcf_struct.get("alt", "N")

    if ref == "N" or alt == "N" or gt_value in (".", "./.", ".|."):
        return "N", None, "ambiguous"

    bases = [ref] + alt.split(",")
    gt_parts = gt_value.replace("|", "/").split("/")

    # Haploid GT
    if len(gt_parts) == 1:
        a = gt_parts[0]
        if a == "0":
            return ref, None, "reference"
        try:
            idx = int(a)
            allele = bases[idx]
        except:
            return "N", None, "ambiguous"

        # SNPs & indels
        if len(allele) == len(ref):
            if len(ref) == 1:
                return allele, None, "snp"
            else:
                return allele, None, "snp_multi"
        elif len(allele) > len(ref):
            return allele, None, "insertion"
        elif len(allele) < len(ref):
            return allele, None, "deletion"

        return "N", None, "ambiguous"

    # Diploid
    a1, a2 = gt_parts[:2]

    # Homozygous reference
    if a1 == a2 == "0":
        return ref, None, "reference"

    # Homozygous alt
    if a1 == a2:
        try:
            idx = int(a1)
            allele = bases[idx]
        except:
            return "N", None, "ambiguous"

        if len(allele) == len(ref):
            if len(ref) == 1:
                return allele, None, "snp"
            else:
                return allele, None, "snp_multi"
        elif len(allele) > len(ref):
            return allele, None, "insertion"
        elif len(allele) < len(ref):
            return allele, None, "deletion"
        return "N", None, "ambiguous"

    # Heterozygous
    try:
        b1 = bases[int(a1)]
    except:
        b1 = "N"
    try:
        b2 = bases[int(a2)]
    except:
        b2 = "N"

    if len(b1) != len(b2):
        # heterozygous indel
        if len(b1) > len(b2):
            return b1, b2, "het_insertion"
        elif len(b1) < len(b2):
            return b1, b2, "het_deletion"
    return b1, b2, "heterozygous"


# --------------------------------------------------------------------
# VCF LINE PARSER
# --------------------------------------------------------------------

def parse_vcf_line(vcf_line, vcf_struct=None):
    bits = vcf_line.strip().split("\t")
    vcf_info = vcf_struct if vcf_struct is not None else {}

    if vcf_line.startswith("#"):
        return get_vcf_header(vcf_line, vcf_info)

    vcf_info["next"] = 0
    vcf_info["header"] = "N"

    if len(bits) < 9:
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

    fmt_parts = bits[8].split(":")
    vcf_info["format_parts_available"] = {p: True for p in fmt_parts}
    vcf_info["number_of_samples"] = len(bits) - 9

    # sample fields
    for i in range(9, len(bits)):
        isolate_number = i - 9
        sample_data = bits[i]
        sample_parts = sample_data.split(":")
        sample_parts = sample_parts[:len(fmt_parts)]
        sample_dict = dict(zip(fmt_parts, sample_parts))

        for key, val in sample_dict.items():
            vcf_info[f"{key}{isolate_number}"] = val

        gt = sample_dict.get("GT", ".")
        base1, base2, base_type = get_vcf_struct_and_determine_base_type(vcf_info, gt)

        vcf_info[f"{isolate_number}_base1"] = base1
        vcf_info[f"{isolate_number}_base2"] = base2
        vcf_info[f"base_type{isolate_number}"] = base_type

        if base_type == "heterozygous":
            code = IUPAC_CODES.get(frozenset([base1, base2]), "N")
            vcf_info[f"amb_char{isolate_number}"] = code

    return vcf_info


# --------------------------------------------------------------------
# BASIC FILTERING
# --------------------------------------------------------------------

def filter_vcf(vcf_file, min_depth=20, remove_ambiguous=False):
    vcf_struct = {}
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("##") or line.startswith("#CHROM"):
                print(line.strip())
                continue

            try:
                vcf_info = parse_vcf_line(line, vcf_struct)
            except Exception:
                continue

            fmt = vcf_info.get("format", "")
            parts = fmt.split(":")
            if "DP" not in parts:
                print(line.strip())
                continue
            dp_index = parts.index("DP")

            sample_ok = True
            for i in range(vcf_info["number_of_samples"]):
                dp = vcf_info.get(f"DP{i}", "0")
                try:
                    if int(dp) < min_depth:
                        sample_ok = False
                        break
                except:
                    sample_ok = False
                    break

            if not sample_ok:
                continue

            if remove_ambiguous and (
                "N" in vcf_info.get("ref", "") or "N" in vcf_info.get("alt", "")
            ):
                continue

            print(line.strip())


# --------------------------------------------------------------------
# REMOVE REFERENCE SITES
# --------------------------------------------------------------------

def remove_ref_sites(vcf_file):
    vcf_struct = {}
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("##") or line.startswith("#CHROM"):
                print(line.strip())
                continue
            vcf_info = parse_vcf_line(line, vcf_struct)
            if vcf_info.get("header") == "Y":
                continue
            n = vcf_info["number_of_samples"]
            keep = False
            for i in range(n):
                base_type = vcf_info.get(f"base_type{i}", "unknown")
                if base_type not in ("reference", "ambiguous", "unknown"):
                    keep = True
                    break
            if keep:
                print(line.strip())


# --------------------------------------------------------------------
# BASIC LOH HELPERS
# --------------------------------------------------------------------

def _merge_simple(pos_list, max_gap=100):
    if not pos_list:
        return []
    pos_list = sorted(pos_list)
    blocks = []
    start = prev = pos_list[0]

    for p in pos_list[1:]:
        if p - prev <= max_gap:
            prev = p
        else:
            blocks.append((start, prev))
            start = p
            prev = p

    blocks.append((start, prev))
    return blocks


def _invert_simple(blocks, chrom_length):
    if not blocks:
        return [(0, chrom_length - 1)]
    blocks = sorted(blocks)
    loh = []
    cur = 0
    for (s, e) in blocks:
        if s > cur:
            loh.append((cur, s - 1))
        cur = e + 1
    if cur <= chrom_length - 1:
        loh.append((cur, chrom_length - 1))
    return loh


def _overlap_interval(a, b):
    s1, e1 = a
    s2, e2 = b
    s = max(s1, s2)
    e = min(e1, e2)
    return max(0, e - s + 1)


def _sum_overlap_with_blocks(win, blocks):
    tot = 0
    for b in blocks:
        tot += _overlap_interval(win, b)
    return tot


def _make_windows_zero_based(seq_len, window_size):
    """Make 0-based, inclusive windows: (0,999), (1000,1999), ..."""
    windows = []
    start = 0
    while start < seq_len:
        end = min(start + window_size - 1, seq_len - 1)
        windows.append((start, end))
        start += window_size
    return windows


# --------------------------------------------------------------------
# FIXED-SIZE WINDOWED VARIANT & LOH STATS
# --------------------------------------------------------------------

def windowed_loh_stats(
    vcf_file,
    fasta_file,
    window_size=10_000,
    het_gap=100,
    min_loh_size=0,
    include_het_indels=True,
):
    """
    Produce per-window variant tallies + LOH blocks + LOH bp.
    Windows are 0-based inclusive.
    """

    fasta_seqs = parse_fasta(fasta_file)
    chrom_lengths = {c: len(s) for c, s in fasta_seqs.items()}

    chrom_windows = {
        chrom: _make_windows_zero_based(chrom_lengths[chrom], window_size)
        for chrom in fasta_seqs
    }

    per_sample_windows = {}
    het_positions = {}

    vcf_struct = {}
    sample_names = None

    # ------------------------------------------------------------
    # FIRST PASS: fill variant stats + heterozygous positions
    # ------------------------------------------------------------
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                vcf_struct = get_vcf_header(line, vcf_struct)
                if "isolate_names" in vcf_struct:
                    sample_names = [
                        vcf_struct["isolate_names"][i]
                        for i in sorted(vcf_struct["isolate_names"])
                    ]
                continue

            vcf_info = parse_vcf_line(line, vcf_struct)
            if vcf_info.get("header") == "Y":
                continue

            chrom = vcf_info["supercontig"]
            if chrom not in chrom_lengths:
                continue

            n = vcf_info["number_of_samples"]
            if sample_names is None:
                sample_names = [f"sample{i}" for i in range(n)]

            # lazy init per-sample per-chrom windows
            for si in range(n):
                sname = sample_names[si]
                per_sample_windows.setdefault(sname, {})
                het_positions.setdefault(sname, {})
                if chrom not in per_sample_windows[sname]:
                    per_sample_windows[sname][chrom] = []
                    for idx, (ws, we) in enumerate(chrom_windows[chrom], start=1):
                        per_sample_windows[sname][chrom].append({
                            "contig": chrom,
                            "win_start": ws,
                            "win_end": we,
                            "window_id": idx,
                            "het_snps": 0,
                            "hom_snps": 0,
                            "het_indels": 0,
                            "hom_indels": 0,
                            "total_variants": 0,
                            "het_sites": 0,
                            "loh_blocks": 0,
                            "loh_bp": 0,
                        })
                if chrom not in het_positions[sname]:
                    het_positions[sname][chrom] = []

            pos = int(vcf_info["position"])
            pos0 = pos - 1
            widx = pos0 // window_size

            for si in range(n):
                sname = sample_names[si]
                base_type = vcf_info.get(f"base_type{si}", "unknown")

                if base_type in ("reference", "ambiguous", "unknown"):
                    continue

                windows = per_sample_windows[sname][chrom]
                if widx >= len(windows):
                    continue

                row = windows[widx]

                if base_type.startswith("snp") and base_type != "heterozygous":
                    row["hom_snps"] += 1
                    row["total_variants"] += 1

                elif base_type == "heterozygous":
                    row["het_snps"] += 1
                    row["total_variants"] += 1
                    row["het_sites"] += 1
                    het_positions[sname][chrom].append(pos0)

                elif base_type in ("insertion", "deletion"):
                    row["hom_indels"] += 1
                    row["total_variants"] += 1

                elif base_type in ("het_insertion", "het_deletion"):
                    row["het_indels"] += 1
                    row["total_variants"] += 1
                    if include_het_indels:
                        row["het_sites"] += 1
                        het_positions[sname][chrom].append(pos0)

    # ------------------------------------------------------------
    # SECOND PASS: merge het → het blocks → LOH blocks → annotate windows
    # ------------------------------------------------------------
    het_blocks = {}
    loh_blocks = {}

    for sname in per_sample_windows:
        het_blocks[sname] = {}
        loh_blocks[sname] = {}

        for chrom in per_sample_windows[sname]:
            seq_len = chrom_lengths[chrom]
            pos_list = het_positions[sname][chrom]

            hblocks = _merge_simple(pos_list, max_gap=het_gap)
            het_blocks[sname][chrom] = hblocks

            loh = _invert_simple(hblocks, seq_len)
            if min_loh_size > 0:
                loh = [(s, e) for (s, e) in loh if (e - s + 1) >= min_loh_size]
            loh_blocks[sname][chrom] = loh

            for row in per_sample_windows[sname][chrom]:
                win = (row["win_start"], row["win_end"])

                row["loh_bp"] = _sum_overlap_with_blocks(win, loh)

                count_blocks = 0
                for b in loh:
                    if _overlap_interval(win, b) > 0:
                        count_blocks += 1
                row["loh_blocks"] = count_blocks

    return {
        "windows": per_sample_windows,
        "het_blocks": het_blocks,
        "loh_blocks": loh_blocks,
        "params": {
            "window_size": window_size,
            "het_gap": het_gap,
            "min_loh_size": min_loh_size,
            "include_het_indels": include_het_indels,
        }
    }


# --------------------------------------------------------------------
# Output Tab WRITER
# --------------------------------------------------------------------

def write_windowed_loh_tsv(window_stats, output_prefix="loh_windows"):
    windows_by_sample = window_stats["windows"]

    for sample, chrom_dict in windows_by_sample.items():
        outname = f"{output_prefix}.{sample}.tab"
        with open(outname, "w") as out:
            out.write(
                "Contig\tWindow_Start\tWindow_Stop\tRunning_Position\t"
                "Het_SNPs\tHom_SNPs\tHet_Indels\tHom_Indels\t"
                "Total_Variants\tHet_Sites\tLOH_Blocks\tLOH_bp\n"
            )

            for chrom, win_list in chrom_dict.items():
                for row in win_list:
                    out.write(
                        f"{row['contig']}\t"
                        f"{row['win_start']}\t{row['win_end']}\t{row['window_id']}\t"
                        f"{row['het_snps']}\t{row['hom_snps']}\t"
                        f"{row['het_indels']}\t{row['hom_indels']}\t"
                        f"{row['total_variants']}\t{row['het_sites']}\t"
                        f"{row['loh_blocks']}\t{row['loh_bp']}\n"
                    )
        print(f"[INFO] wrote {outname}")

