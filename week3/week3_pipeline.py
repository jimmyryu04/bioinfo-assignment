#!/usr/bin/env python3
"""
Week 3 — Mission 3 (Bioinformatics 1, SNU 2026 Spring)
CIMS (Crosslinking-induced mutation sites) analysis via per-position
Shannon entropy on CLIP-seq pileups.

This script ports week3.ipynb to a single runnable file and finishes the
assignment at the end of the notebook:

    1. Count base calls at each position.
    2. Compute Shannon entropy per position.
    3. Emit a 4-column bedGraph for the UCSC Genome Browser (mm39).
    4. Repeat for Mirlet7g, Mirlet7f-1, and Mirlet7d.

Step 5 onward in the notebook (UCSC upload / PDF screenshot) is a manual
browser action — the bedGraph files produced here are what gets uploaded.

Intermediate files are cached in WORK_DIR so re-runs only redo what is
missing.
"""

import math
import os
import re
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration

DATA_DIR = Path(
    "/qbio/ryu/project/etc/bioinformatics/assignments/week1/work/binfo1-datapack1"
)
WORK_DIR = Path("/qbio/ryu/project/etc/bioinformatics/assignments/week3/work")
WORK_DIR.mkdir(exist_ok=True, parents=True)

SAMTOOLS = "/blaze/ryu/conda/envs/lab/bin/samtools"
CLIP_BAM = DATA_DIR / "CLIP-35L33G.bam"
GTF = DATA_DIR / "gencode.gtf"

# Three miRNA loci on mm39 (from gencode.gtf gene entries).
GENES = [
    {"name": "Mirlet7g",   "chrom": "chr9",  "start": 106056039, "end": 106056126},
    {"name": "Mirlet7d",   "chrom": "chr13", "start": 48689488,  "end": 48689590},
    {"name": "Mirlet7f-1", "chrom": "chr13", "start": 48691305,  "end": 48691393},
]


# ---------------------------------------------------------------------------
# Helpers


def sh(cmd: str) -> None:
    print(f"$ {cmd}", flush=True, file=sys.stderr)
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def cached(path: Path) -> bool:
    if path.exists() and path.stat().st_size > 0:
        print(f"[skip] {path}", file=sys.stderr)
        return True
    return False


# ---------------------------------------------------------------------------
# Stage 1 — subset CLIP BAM to the gene region


def subset_bam(gene: dict, out_bam: Path) -> None:
    if cached(out_bam):
        return
    region = f"{gene['chrom']}:{gene['start']}-{gene['end']}"
    sh(f"{SAMTOOLS} view -b -o {out_bam} {CLIP_BAM} {region}")


# ---------------------------------------------------------------------------
# Stage 2 — mpileup over the region BAM


def mpileup(in_bam: Path, out_pileup: Path) -> None:
    if cached(out_pileup):
        return
    sh(f"{SAMTOOLS} mpileup {in_bam} > {out_pileup}")


# ---------------------------------------------------------------------------
# Stage 3 — clip pileup to the gene's exact coordinates
#
# `samtools mpileup` emits every position covered by any overlapping read,
# so the raw pileup runs well past the 88 bp gene. The notebook uses awk to
# trim it down; we do the same here in pure Python.


def clip_pileup_to_gene(in_pileup: Path, gene: dict, out_pileup: Path) -> None:
    if cached(out_pileup):
        return
    lo, hi = gene["start"], gene["end"]
    with open(in_pileup) as f, open(out_pileup, "w") as out:
        for line in f:
            parts = line.split("\t", 3)
            if len(parts) < 2:
                continue
            pos = int(parts[1])
            if lo <= pos <= hi:
                out.write(line)


# ---------------------------------------------------------------------------
# Stage 4 — parse a pileup "matches" string into base counts
#
# The pileup matches column is a packed format. To get clean per-base counts
# we have to skip past structural characters that don't represent a base
# call at THIS position:
#   ^X    start of a read; X is the mapping quality (skip both chars)
#   $     end of a read (1 char)
#   *, #  deletion placeholder
#   <, >  reference skip (e.g. spliced read)
#   +N[bases]  insertion right after this position (skip the run)
#   -N[bases]  deletion right after this position (skip the run)
# What remains is one character per aligned read:
#   .  or  ,   match to the reference (forward / reverse strand)
#   A C G T (or lowercase)   mismatch base call
#   N  or  n   ambiguous
#
# Since the BAM was mpileup'd without a -f reference, the ref column is "N"
# and matches are still represented as ./, — that's fine for entropy, where
# we only need a category for each distinguishable outcome.


_INDEL_RE = re.compile(r"[+-](\d+)")


def parse_pileup_bases(s: str) -> dict:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "match": 0, "N": 0}
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "^":
            i += 2  # skip the mapping-quality char too
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            m = _INDEL_RE.match(s, i)
            if m is None:
                i += 1
                continue
            indel_len = int(m.group(1))
            i = m.end() + indel_len
            continue
        if c in "*<>#":
            i += 1
            continue
        if c in ".,":
            counts["match"] += 1
        else:
            up = c.upper()
            if up in "ACGT":
                counts[up] += 1
            elif up == "N":
                counts["N"] += 1
            # anything else is silently ignored
        i += 1
    return counts


# ---------------------------------------------------------------------------
# Stage 5 — Shannon entropy from counts


def shannon_entropy(counts: dict) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


# ---------------------------------------------------------------------------
# Stage 6 — write a 4-column bedGraph
#
# bedGraph fields: chrom, start (0-based), end (exclusive), value
# Pileup positions are 1-based; convert with start = pos-1, end = pos.


def write_bedgraph(in_pileup: Path, gene: dict, out_bg: Path) -> None:
    track_name = f"ShannonH_{gene['name']}"
    description = f"Per-base Shannon entropy from CLIP-35L33G pileup at {gene['name']}"
    with open(in_pileup) as f, open(out_bg, "w") as out:
        out.write(
            f'track type=bedGraph name="{track_name}" '
            f'description="{description}" visibility=full color=0,0,200\n'
        )
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            basereads = parts[4]
            counts = parse_pileup_bases(basereads)
            h = shannon_entropy(counts)
            out.write(f"{chrom}\t{pos - 1}\t{pos}\t{h:.6f}\n")


# ---------------------------------------------------------------------------
# Driver


def run_gene(gene: dict) -> Path:
    name = gene["name"]
    print(f"\n=== {name} ({gene['chrom']}:{gene['start']}-{gene['end']}) ===")

    region_bam = WORK_DIR / f"CLIP-{name}.bam"
    pileup_all = WORK_DIR / f"CLIP-{name}.pileup"
    pileup_gene = WORK_DIR / f"CLIP-{name}-gene.pileup"
    bedgraph = WORK_DIR / f"CLIP-{name}-entropy.bedgraph"

    print("[1] Subset CLIP BAM to gene region")
    subset_bam(gene, region_bam)

    print("[2] Run samtools mpileup")
    mpileup(region_bam, pileup_all)

    print("[3] Clip pileup to gene coordinates")
    clip_pileup_to_gene(pileup_all, gene, pileup_gene)

    print("[4] Parse bases, compute Shannon entropy, write bedGraph")
    write_bedgraph(pileup_gene, gene, bedgraph)

    n_pos = sum(1 for _ in open(bedgraph)) - 1  # minus header line
    print(f"     wrote {bedgraph} ({n_pos} positions)")
    return bedgraph


def main() -> None:
    os.chdir(WORK_DIR)
    print(f"Working in {WORK_DIR}", file=sys.stderr)

    outputs = [run_gene(g) for g in GENES]

    print("\n=== Done ===")
    print("BedGraph files (upload via UCSC 'add custom tracks' on mm39):")
    for p in outputs:
        print(f"  {p}")


if __name__ == "__main__":
    main()
