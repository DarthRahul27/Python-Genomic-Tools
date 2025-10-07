# r.anand@exeter.ac.uk
# Module for fasta file manipulation and getting fasta stats

def parse_fasta(filepath):
    print(f"parse_fasta {filepath} ...")
    sequences = {}
    header = None
    seq_chunks = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            sequences[header] = "".join(seq_chunks)
    print(f"parse_fasta {filepath} complete!")
    return sequences


def gc_content(seq):
    gc = sum(1 for base in seq if base in "GCgc")
    return (gc / len(seq)) * 100 if seq else 0.0


def count_ambiguous(seq):
    return sum(1 for base in seq if base in "Nn")


def fasta_stats(records):
    per_seq = {}
    lengths, ambigs = [], []

    for h, seq in records.items():
        length = len(seq)
        gc = gc_content(seq)
        n_count = count_ambiguous(seq)

        lengths.append(length)
        ambigs.append(n_count)
        per_seq[h] = {"length": length, "gc": gc, "ambiguous": n_count}

    if not lengths:
        overall = {"count": 0, "total_length": 0, "avg_len": 0, "max_len": 0, "ambiguous": 0}
    else:
        total_len = sum(lengths)
        total_ambig = sum(ambigs)
        overall = {
            "count": len(lengths),
            "total_len": total_len,
            "avg_len": total_len / len(lengths),
            "max_len": max(lengths),
            "ambiguous": total_ambig,
        }

    return {"per_sequence": per_seq, "overall": overall}


def print_fasta_stats(records):
    stats = fasta_stats(records)
    print("Chromosome\tLength\tGC%\tAmbiguous")
    for h, s in stats["per_sequence"].items():
        print(f"{h}\t{s['length']}\t{round(s['gc'],2)}\t{s['ambiguous']}")
    print("\nOverall stats")
    overall = stats["overall"]
    print(f"Sequences\t{overall['count']}")
    print(f"Total length\t{overall['total_len']}")
    print(f"Average length\t{round(overall['avg_len'],2)}")
    print(f"Max length\t{overall['max_len']}")
    print(f"Total ambiguous\t{overall['ambiguous']}")


CODON_TABLE = {
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F',
    'TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y',
    'TAA':'*','TAG':'*','TGA':'*',
    'TGC':'C','TGT':'C',
    'TGG':'W',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H',
    'CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'ATA':'I','ATC':'I','ATT':'I',
    'ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N',
    'AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S',
    'AGA':'R','AGG':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D',
    'GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G'
}

def translate_sequences(seq, frame=0):
    """Translate nucleotide sequence into amino acids using standard codon table."""
    seq = seq.upper().replace("U", "T")
    aa_seq = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "X")
        aa_seq.append(aa)
    return "".join(aa_seq)



