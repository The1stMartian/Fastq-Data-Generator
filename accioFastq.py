#!/usr/bin/env python3
"""
Generates a tiny FASTQ dataset (3 control + 3 test) from a template sequence. Fastq options include:
  - Paired-end (PE) or single-end (SE) reads
  - Optional in-read barcodes:
        * UMI only (prefix)
        * Cell barcode + UMI (prefix CB then UMI)
        * No barcodes (set lengths to 0)
  - Realistic Phred qualities that decline near BOTH 5' and 3' ends

R1 (PE):  [CB][UMI][biological sequence]
R2 (PE):  [biological sequence]                 (optionally also prefixed if --barcode-in both)
SE:       [CB][UMI][biological sequence]

# Flag Details:
    --fasta yourSequencechrM.fa 		# Nucleotide sequence input
    --outdir fastq 						# output path
    --layout PE 						# PE or SE)
    --reads 5000 						# Number of reads
    --read-len 75 						# read length
    --cb-len 0 							# 0 = no cell bc, or X nts
    --umi-len 8 						# 0=no UMI, or X nts
    --up-region 3307:4262 				# Up regulated region (nts)
    --down-region 9207:9990 			# Down regulated region (nts)
    --up-name MT-ND1 					# Manually entered region name
    --down-name MT-CO3 					# Manually entered region name
    --up-mult 4 						# One gene is upregulated by X multiplier vs control
    --down-mult 0.05 					# One gene is down regulated by x multiplier vs control
    --seed 42 							# random generator
    --verbose							# Script output progress


# Usage example:
python makeFqData.py \
    --fasta chrM.fa \
    --outdir fastq \
    --layout PE \
    --reads 5000 \
    --read-len 75 \
    --umi-len 10 \
    --cb-len 0 \
    --up-region 3307:4262 \
    --down-region 9207:9990 \
    --up-name MT-ND1 \
    --down-name MT-CO3 \
    --up-mult 4 \
    --down-mult 0.05 \
    --seed 42 \
    --verbose
"""

import argparse
import gzip
import os
import random
import sys
from pathlib import Path
from typing import Tuple, List

DNA = "ACGT"

# ----------------- I/O helpers -----------------

def read_fasta_single(fa_path: str) -> str:
    if not os.path.exists(fa_path):
        raise FileNotFoundError(f"FASTA not found: {fa_path}")
    seq = []
    with open(fa_path, 'r') as f:
        for line in f:
            if not line:
                continue
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
    s = "".join(seq)
    if not s:
        raise ValueError(f"No sequence parsed from FASTA: {fa_path}")
    return s

def write_fastq_gz(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as gz:
        for name, seq, qual in records:
            gz.write(f"@{name}\n{seq}\n+\n{qual}\n")

# ----------------- Sequence helpers -----------------

def revcomp(s: str) -> str:
    tab = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tab)[::-1]

def mutate(seq: str, mut_rate: float, rng: random.Random) -> str:
    if mut_rate <= 0:
        return seq
    out = []
    for c in seq:
        if c not in "ACGT":
            out.append(c)
            continue
        if rng.random() < mut_rate:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
        else:
            out.append(c)
    return "".join(out)

# Circular genome fragment

def sample_fragment(ref: str, frag_len: int, start: int) -> str:
    L = len(ref)
    idxs = [(start + i) % L for i in range(frag_len)]
    return "".join(ref[i] for i in idxs)

# ----------------- Regions & sampling -----------------

def parse_region(arg: str) -> Tuple[int, int]:
    # "start:end" -> (int start, int end), 0-based half-open
    s, e = arg.split(":")
    s, e = int(s), int(e)
    if not (0 <= s < e):
        raise ValueError(f"Bad region spec '{arg}' (need 0-based half-open start:end with start<end)")
    return s, e

def region_len(region: Tuple[int, int]) -> int:
    s, e = region
    return e - s

# circular interval overlap: does [a, a+alen) intersect [b, b+blen) on circle of length L?

def circ_overlap(a: int, alen: int, b: int, blen: int, L: int) -> bool:
    def segs(start, length):
        end = start + length
        if end <= L:
            return [(start, end)]
        else:
            return [(start, L), (0, end % L)]
    for A0, A1 in segs(a, alen):
        for B0, B1 in segs(b, blen):
            if A0 < B1 and B0 < A1:
                return True
    return False

# choose bucket only (bg/up/down) by weights

def pick_bucket(w_bg: float, w_up: float, w_down: float, rng: random.Random) -> str:
    buckets = ["bg", "up", "down"]
    weights = [max(0.0, w_bg), max(0.0, w_up), max(0.0, w_down)]
    return rng.choices(buckets, weights=weights, k=1)[0]

# choose a start so that an interval of 'span' is fully contained in 'region'

def choose_start_in_region(region: Tuple[int, int], span: int, rng: random.Random) -> int:
    s, e = region
    Lr = e - s
    if span > Lr:
        # cannot fully contain; fall back to region start (will spill)
        return s
    return s + rng.randrange(0, Lr - span + 1)

# choose a start so an interval of 'span' on the circle avoids both regions

def choose_start_in_bg(genome_len: int, up: Tuple[int, int], dn: Tuple[int, int], span: int, rng: random.Random) -> int:
    # rejection sampling; fast because regions are small relative to genome
    for _ in range(10000):
        st = rng.randrange(0, genome_len)
        if not circ_overlap(st, span, up[0], region_len(up), genome_len) and \
           not circ_overlap(st, span, dn[0], region_len(dn), genome_len):
            return st
    # extremely unlikely; fallback
    return 0

# ----------------- Quality modeling -----------------

def clamp(v, lo, hi):
    return max(lo, min(hi, v))

def phred_char(q_int: float) -> str:
    q = clamp(int(round(q_int)), 3, 40)  # Phred+33, clamp 3..40
    return chr(q + 33)

def u_shaped_qual_profile(n: int, rng: random.Random, center_q: float = 33,
                          end_drop: float = 10, jitter_sd: float = 2.0,
                          low_q_prob: float = 0.04) -> List[float]:
    """U-shaped quality: higher mid-read, lower at both ends, with jitter and rare low-Q spikes."""
    if n <= 1:
        return [center_q]
    quals = []
    mid = (n - 1) / 2.0
    for i in range(n):
        c = abs((i - mid) / mid)  # 0 at center, 1 at ends
        q = center_q - end_drop * (c ** 1.2) + rng.gauss(0, jitter_sd)
        if rng.random() < low_q_prob:
            q = rng.uniform(3, 8)
        quals.append(clamp(q, 3, 40))
    return quals

def quals_to_string(quals: List[float]) -> str:
    return "".join(phred_char(q) for q in quals)

# ----------------- Main -----------------

def main():
    ap = argparse.ArgumentParser(description="Generate small chrM FASTQ set with optional in-read UMI/CB; SE or PE; realistic qualities.")
    ap.add_argument("--fasta", required=True, help="FASTA with chrM sequence")
    ap.add_argument("--outdir", default="fastq", help="Output directory")
    ap.add_argument("--layout", choices=["PE", "SE"], default="PE", help="Paired-end (PE) or single-end (SE)")
    ap.add_argument("--reads", type=int, default=2000, help="PE: read PAIRS per sample; SE: reads per sample")
    ap.add_argument("--read-len", type=int, default=75, help="Biological read length (without barcodes)")

    # Barcode options (in-sequence only; NEVER in headers)
    ap.add_argument("--umi-len", type=int, default=0, help="UMI length (prefix). Set 0 to disable.")
    ap.add_argument("--cb-len", type=int, default=0, help="Cell barcode length (prefix before UMI). Set 0 to disable.")
    ap.add_argument("--barcode-in", choices=["R1", "both"], default="R1",
                    help="For PE only, place CB/UMI in R1 (typical) or in both reads")

    # Fragment / error model
    ap.add_argument("--ins-mean", type=int, default=200)
    ap.add_argument("--ins-sd", type=float, default=20)
    ap.add_argument("--mut-rate-ctrl", type=float, default=0.001)
    ap.add_argument("--mut-rate-test", type=float, default=0.003)

    # Regions controlling up/down genes
    ap.add_argument("--up-region", required=True, help="0-based half-open region for 'up' gene, e.g. 3307:4262")
    ap.add_argument("--down-region", required=True, help="0-based half-open region for 'down' gene, e.g. 9207:9990")
    ap.add_argument("--up-name", default="gene_up")
    ap.add_argument("--down-name", default="gene_down")
    ap.add_argument("--weight-bg", type=float, default=98.0)
    ap.add_argument("--weight-up", type=float, default=1.0)
    ap.add_argument("--weight-down", type=float, default=1.0)
    ap.add_argument("--up-mult", type=float, default=4.0, help="Multiply up weight in TEST")
    ap.add_argument("--down-mult", type=float, default=0.25, help="Multiply down weight in TEST")

    # Quality profile knobs (optional)
    ap.add_argument("--center-q", type=float, default=33.0, help="Typical mid-read quality")
    ap.add_argument("--end-drop", type=float, default=10.0, help="Drop from center to ends")
    ap.add_argument("--jitter-sd", type=float, default=2.0, help="Gaussian jitter SD for per-base Q")
    ap.add_argument("--low-q-prob", type=float, default=0.04, help="Chance of very low-Q spike per base")

    # Misc
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--verbose", action="store_true", help="Print progress messages")

    args = ap.parse_args()

    rng = random.Random(args.seed)

    # Read reference
    try:
        chrM = read_fasta_single(args.fasta)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    L = len(chrM)

    try:
        region_up = parse_region(args.up_region)
        region_down = parse_region(args.down_region)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    # sanity: warn if region shorter than read or typical insert (enrichment will degrade)
    if region_len(region_up) < args.read_len:
        print(f"[WARN] Up-region ({region_up}) shorter than read_len={args.read_len}; reads may spill outside region.", file=sys.stderr)
    if region_len(region_down) < args.read_len:
        print(f"[WARN] Down-region ({region_down}) shorter than read_len={args.read_len}; reads may spill outside region.", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "regions.bed", "w") as bed:
        bed.write(f"chrM\t{region_up[0]}\t{region_up[1]}\t{args.up_name}\n")
        bed.write(f"chrM\t{region_down[0]}\t{region_down[1]}\t{args.down_name}\n")

    # Samples: (name, is_test, mut_rate)
    samples = (
        [("ctrl1", False, args.mut_rate_ctrl), ("ctrl2", False, args.mut_rate_ctrl), ("ctrl3", False, args.mut_rate_ctrl)] +
        [("test1", True,  args.mut_rate_test), ("test2", True,  args.mut_rate_test), ("test3", True,  args.mut_rate_test)]
    )

    def make_barcode() -> str:
        cb = "".join(rng.choice(DNA) for _ in range(args.cb_len)) if args.cb_len > 0 else ""
        umi = "".join(rng.choice(DNA) for _ in range(args.umi_len)) if args.umi_len > 0 else ""
        return cb + umi

    if args.verbose:
        print(f"[INFO] layout={args.layout} reads={args.reads} read_len={args.read_len} CB={args.cb_len} UMI={args.umi_len}")
        print(f"[INFO] Up={args.up_name} {args.up_region}  Down={args.down_name} {args.down_region}")

    for sample, is_test, mut_rate in samples:
        region_counts = {"bg": 0, "up": 0, "down": 0}
        w_bg = args.weight_bg
        w_up = args.weight_up * (args.up_mult if is_test else 1.0)
        w_dn = args.weight_down * (args.down_mult if is_test else 1.0)

        if args.layout == "SE":
            out_se = outdir / f"{sample}.fastq.gz"
            with gzip.open(out_se, "wt") as gz:
                for i in range(1, args.reads + 1):
                    span = args.read_len
                    bucket = pick_bucket(w_bg, w_up, w_dn, rng)

                    if bucket == "up":
                        start = choose_start_in_region(region_up, span, rng)
                    elif bucket == "down":
                        start = choose_start_in_region(region_down, span, rng)
                    else:
                        start = choose_start_in_bg(L, region_up, region_down, span, rng)

                    region_counts[bucket] += 1
                    frag = sample_fragment(chrM, span, start)

                    bio = mutate(frag, mut_rate, rng)
                    barcode = make_barcode()
                    seq = barcode + bio

                    quals = u_shaped_qual_profile(len(seq), rng, args.center_q, args.end_drop, args.jitter_sd, args.low_q_prob)
                    qual_str = quals_to_string(quals)

                    qname = f"{sample}:read{i}"
                    gz.write(f"@{qname}\n{seq}\n+\n{qual_str}\n")

                    if args.verbose and (i % max(1, args.reads // 10) == 0):
                        print(f"[SE] {sample} {i}/{args.reads}", flush=True)

            if args.verbose:
                print(f"[DONE] {sample} -> {out_se}")

        else:
            # Paired-end: write R1 & R2 in single pass with identical QNAMEs
            r1_path = outdir / f"{sample}_R1.fastq.gz"
            r2_path = outdir / f"{sample}_R2.fastq.gz"
            with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
                for i in range(1, args.reads + 1):
                    ins = max(2 * args.read_len + 10, int(rng.gauss(args.ins_mean, args.ins_sd)))
                    bucket = pick_bucket(w_bg, w_up, w_dn, rng)

                    if bucket == "up":
                        ins_eff = min(ins, region_len(region_up))
                        start = choose_start_in_region(region_up, ins_eff, rng)
                    elif bucket == "down":
                        ins_eff = min(ins, region_len(region_down))
                        start = choose_start_in_region(region_down, ins_eff, rng)
                    else:
                        ins_eff = ins
                        start = choose_start_in_bg(L, region_up, region_down, ins_eff, rng)

                    region_counts[bucket] += 1
                    frag = sample_fragment(chrM, ins_eff, start)

                    r1_bio = mutate(frag[:args.read_len], mut_rate, rng)
                    r2_bio = mutate(revcomp(frag[-args.read_len:]), mut_rate, rng)

                    barcode = make_barcode()
                    r1_seq = barcode + r1_bio
                    if args.barcode_in == "both" and (args.cb_len + args.umi_len) > 0:
                        r2_seq = barcode + r2_bio
                    else:
                        r2_seq = r2_bio

                    r1_q = u_shaped_qual_profile(len(r1_seq), rng, args.center_q, args.end_drop, args.jitter_sd, args.low_q_prob)
                    r2_q = u_shaped_qual_profile(len(r2_seq), rng, args.center_q, args.end_drop, args.jitter_sd, args.low_q_prob)
                    r1_qual = quals_to_string(r1_q)
                    r2_qual = quals_to_string(r2_q)

                    qname = f"{sample}:read{i}"
                    r1.write(f"@{qname}\n{r1_seq}\n+\n{r1_qual}\n")
                    r2.write(f"@{qname}\n{r2_seq}\n+\n{r2_qual}\n")

                    if args.verbose and (i % max(1, args.reads // 10) == 0):
                        print(f"[PE] {sample} {i}/{args.reads}", flush=True)

            if args.verbose:
                print(f"[DONE] {sample} -> {r1_path}, {r2_path}")

        # Per-sample region counts (based on bucket choices)
        rc_tsv = outdir / f"{sample}.region_counts.tsv"
        with open(rc_tsv, "w") as t:
            t.write("bucket\treads\n")
            t.write(f"bg\t{region_counts['bg']}\n")
            t.write(f"{args.up_name}\t{region_counts['up']}\n")
            t.write(f"{args.down_name}\t{region_counts['down']}\n")

    if args.verbose:
        print(f"Regions BED: {outdir/'regions.bed'}")
    print("Done.")

if __name__ == "__main__":
    main()
