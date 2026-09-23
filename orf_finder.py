"""Find candidate bacterial genes (ORFs) with regular expressions, translate them,
and compare each call with the CDS features of a GenBank record."""
import re

GENBANK_FILE = "data/J01636_lac_operon.gb"
MIN_CODONS = 100
STARTS = "ATG|GTG|TTG"   # try "ATG" to see which lac genes you lose

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16*i + 4*j + k]
               for i, a in enumerate(BASES)
               for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}

# start, then MIN_CODONS or more non-stop codons, then a stop; the lookahead allows overlaps
ORF_RE = re.compile(r"(?=((?:%s)(?:(?!TAA|TAG|TGA)[ACGT]{3}){%d,}(?:TAA|TAG|TGA)))" % (STARTS, MIN_CODONS))
# Shine-Dalgarno core, a 4-12 base spacer, then the start codon at the very end
RBS_RE = re.compile(r"(?:GGAG|GAGG|AGGA)[ACGT]{4,12}(?:%s)$" % STARTS)
CDS_RE = re.compile(r"^ {5}CDS +(complement\()?<?(\d+)\.\.>?(\d+)\)?(.*?)(?=^ {5}\S|^ORIGIN)", re.M | re.S)


def read_genbank(path):
    """Return the sequence and a dict of CDS features {(strand, start, end): gene}."""
    with open(path) as fh:
        text = fh.read()
    origin = text.split("\nORIGIN", 1)[1].split("\n", 1)[1]   # skip the ORIGIN line itself
    seq = re.sub(r"[^ACGT]", "", origin.split("//")[0].upper())
    cds = {}
    for m in CDS_RE.finditer(text):
        gene = re.search(r'/gene="([^"]+)"', m.group(4))
        strand = "-" if m.group(1) else "+"
        cds[(strand, int(m.group(2)), int(m.group(3)))] = gene.group(1) if gene else "CDS"
    return seq, cds


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(orf):
    codons = re.findall(r"[ACGT]{3}", orf)
    return "M" + "".join(CODON_TABLE.get(c, "X") for c in codons[1:])   # any start codon reads as fMet


def has_rbs(s, start):
    return bool(RBS_RE.search(s[max(0, start - 30): start + 3]))


def find_orfs(seq):
    """One ORF per stop codon: the most upstream start with an RBS, else the most upstream start."""
    length = len(seq)
    orfs = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        starts_by_stop = {}
        for m in ORF_RE.finditer(s):
            starts_by_stop.setdefault(m.start(1) + len(m.group(1)), []).append(m.start(1))
        for end, starts in starts_by_stop.items():
            with_rbs = [p for p in starts if has_rbs(s, p)]
            start = with_rbs[0] if with_rbs else starts[0]
            orf = s[start:end]
            if strand == "+":
                a, b = start + 1, end                    # 1-based, like GenBank
            else:
                a, b = length - end + 1, length - start
            orfs.append((a, b, strand, orf, bool(with_rbs)))
    return sorted(orfs)


def compare(orf_call, cds):
    a, b, strand = orf_call[:3]
    if (strand, a, b) in cds:
        return "= " + cds[(strand, a, b)]
    for (cs, ca, cb), gene in cds.items():
        same_stop = cb == b if strand == "+" else ca == a
        if cs == strand and same_stop:
            return "~ %s (annotated %d..%d)" % (gene, ca, cb)
    return "not annotated"


if __name__ == "__main__":
    seq, cds = read_genbank(GENBANK_FILE)
    print("%d bp, %d annotated CDS" % (len(seq), len(cds)))
    for orf_call in find_orfs(seq):
        a, b, strand, orf, rbs = orf_call
        protein = translate(orf).rstrip("*")
        print("%5d..%-5d %s  %4d aa  RBS:%-3s  %-32s %s..." % (
            a, b, strand, len(protein), "yes" if rbs else "no", compare(orf_call, cds), protein[:20]))
