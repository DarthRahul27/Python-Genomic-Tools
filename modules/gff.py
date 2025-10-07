# r.anand@exeter.ac.uk
# Module for GFF file manipulation and extracting GFF features

# Global list to collect all stats runs
gff_stats = []


def parse_gff(filepath):
    """
    Parse a GFF3 file and return all records (genes, mRNAs, CDS, exons, etc.).

    :param filepath: str, path to the GFF3 file
    :return: list of dicts with parsed GFF fields
    """
    print(f"parse_gff {filepath} ...")

    gff_records = []

    with open(filepath, "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            bits = line.strip().split("\t")
            if len(bits) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attributes = bits

            # Parse attributes into a dictionary
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value

            record = {
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": int(start),
                "end": int(end),
                "score": None if score == "." else float(score),
                "strand": strand,
                "phase": None if phase == "." else int(phase),
                "attributes": attr_dict,
            }

            gff_records.append(record)

    print(f"parse_gff {filepath} complete! Found {len(gff_records)} records.")
    return gff_records


def filter_gff_by_feature(gff_records, feature):
    """
    Filter GFF records by feature type.
    
    :param gff_records: list of dicts from parse_gff()
    :param feature: str, feature type (e.g., "gene", "mRNA", "CDS", "exon")
    :return: filtered list of records
    """
    return [rec for rec in gff_records if rec["type"].lower() == feature.lower()]


def get_gff_stats(gff_records, feature_name="all"):
    """
    Compute summary statistics from GFF records.
    Appends results to global gff_stats list.
    
    :param gff_records: list of dicts (possibly filtered by feature type)
    :param feature_name: str, label for stats (e.g., "gene", "CDS")
    """
    global gff_stats

    if not gff_records:
        result = {"feature": feature_name, "count": 0}
        gff_stats.append(result)
        return gff_stats

    count = len(gff_records)

    # Chrom distribution
    chrom_distribution = {}
    for rec in gff_records:
        chrom = rec["seqid"]
        chrom_distribution[chrom] = chrom_distribution.get(chrom, 0) + 1

    # Strand distribution
    strand_distribution = {}
    for rec in gff_records:
        strand = rec["strand"]
        if strand in ("+", "-"):
            strand_distribution[strand] = strand_distribution.get(strand, 0) + 1

    # Length stats
    lengths = [rec["end"] - rec["start"] + 1 for rec in gff_records]
    lengths.sort()
    length_min = lengths[0]
    length_max = lengths[-1]
    length_mean = sum(lengths) / len(lengths)

    n = len(lengths)
    if n % 2 == 1:
        length_median = lengths[n // 2]
    else:
        length_median = (lengths[n // 2 - 1] + lengths[n // 2]) / 2

    # Score stats
    scores = [rec["score"] for rec in gff_records if rec["score"] is not None]
    score_min = score_max = score_mean = None
    if scores:
        scores.sort()
        score_min = scores[0]
        score_max = scores[-1]
        score_mean = sum(scores) / len(scores)

    result = {
        "feature": feature_name,
        "count": count,
        "chrom_distribution": chrom_distribution,
        "strand_distribution": strand_distribution,
        "length_min": length_min,
        "length_max": length_max,
        "length_mean": length_mean,
        "length_median": length_median,
        "score_min": score_min,
        "score_max": score_max,
        "score_mean": score_mean,
    }

    gff_stats.append(result)
    return gff_stats
