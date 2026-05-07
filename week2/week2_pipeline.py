#!/usr/bin/env python3
"""
Week 2 — Mission 2 (Bioinformatics 1, SNU 2026 Spring)
Ribosome footprint density near start/stop codons (Figure S5A).

This script ports the pipeline from week2.ipynb to a single runnable file and
extends it to draw Figure S5A. Intermediate files are cached in WORK_DIR so
re-runs only redo missing steps.

Stages:
  1. Filter GTF for start_codon / stop_codon entries on the + strand with
     transcript_support_level "1".
  2. Filter GTF for + strand exons.
  3. Use bedtools intersect to keep exons that contain a start (or stop) codon
     of the same transcript, emitting a small BED with the codon position.
  4. Filter each RPF BAM to + strand reads >= 25 nt.
  5. Compute 5'-end coverage (bedtools genomecov -bg -5).
  6. Intersect coverage with the codon-containing-exon BED.
  7. Aggregate counts by position relative to the codon and plot Figure S5A.
"""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Configuration

DATA_DIR = Path(
    "/qbio/ryu/project/etc/bioinformatics/assignments/week1/work/binfo1-datapack1"
)
WORK_DIR = Path("/qbio/ryu/project/etc/bioinformatics/assignments/week2/work")
WORK_DIR.mkdir(exist_ok=True, parents=True)

BEDTOOLS = "/blaze/ryu/conda/envs/lab/bin/bedtools"
SAMTOOLS = "/blaze/ryu/conda/envs/lab/bin/samtools"

GTF = DATA_DIR / "gencode.gtf"
SAMPLES = ["RPF-siLuc", "RPF-siLin28a"]


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
# Stage 1 — codon GTF (TSL=1, + strand)


def extract_codon_gtf(codon_type: str, out_path: Path) -> None:
    if cached(out_path):
        return
    sh(
        rf"""grep $'\t{codon_type}\t.*\t+\t.*transcript_support_level "1"' {GTF} """
        rf"""| sed -e $'s/\t[^\t]*transcript_id "\\([^"]*\\)".*$/\t\\1/g' """
        rf"""> {out_path}"""
    )


# ---------------------------------------------------------------------------
# Stage 2 — + strand exon GTF


def extract_plus_exon_gtf(out_path: Path) -> None:
    if cached(out_path):
        return
    sh(
        rf"""grep $'\texon\t.*\t+\t' {GTF} """
        rf"""| sed -e $'s/\t[^\t]*transcript_id "\\([^"]*\\)".*$/\t\\1/g' """
        rf"""> {out_path}"""
    )


# ---------------------------------------------------------------------------
# Stage 3 — exons containing the codon (BED with codon coord in col 5)


def make_codon_exon_bed(codon_gtf: Path, exon_gtf: Path, out_bed: Path) -> None:
    if cached(out_bed):
        return
    # Output BED columns:
    #   chr  exon_start(0-based)  exon_end  transcript_id  codon_start(0-based)  strand
    sh(
        rf"""{BEDTOOLS} intersect -a {codon_gtf} -b {exon_gtf} -wa -wb """
        rf"""| awk -F'\t' -v OFS='\t' '$9 == $18 {{ print $10, $13-1, $14, $18, $4-1, $16; }}' """
        rf"""| sort -k1,1 -k2,3n -k4,4 > {out_bed}"""
    )


# ---------------------------------------------------------------------------
# Stage 4 — filter BAM (+ strand, length >= 25 nt)


def filter_bam_plus_long(sample: str, out_bam: Path) -> None:
    if cached(out_bam):
        return
    bam = DATA_DIR / f"{sample}.bam"
    # bioawk -c sam '{ if (length($seq) >= 25) print $0; }'
    # Replaced with plain awk on SAM column 10 (the SEQ field).
    sh(
        rf"""({SAMTOOLS} view -H {bam}; """
        rf""" {SAMTOOLS} view -F20 {bam} """
        rf""" | awk -F'\t' '{{ if (length($10) >= 25) print $0; }}') """
        rf"""| {SAMTOOLS} view -b -o {out_bam}"""
    )


# ---------------------------------------------------------------------------
# Stage 5 — 5'-end bedgraph


def fivep_bedgraph(filtered_bam: Path, out_bg: Path) -> None:
    if cached(out_bg):
        return
    sh(rf"""{BEDTOOLS} genomecov -ibam {filtered_bam} -bg -5 > {out_bg}""")


# ---------------------------------------------------------------------------
# Stage 6 — intersect 5' coverage with codon-containing exons


def intersect_fivep_with_codon(fivep_bg: Path, codon_exon_bed: Path, out_path: Path) -> None:
    if cached(out_path):
        return
    sh(
        rf"""{BEDTOOLS} intersect -a {fivep_bg} -b {codon_exon_bed} """
        rf"""-wa -wb -nonamecheck > {out_path}"""
    )


# ---------------------------------------------------------------------------
# Stage 7 — aggregate counts by relative position


def aggregate_relative_counts(intersect_path: Path, x_min: int, x_max: int):
    """
    For each row in `<bedgraph entry> | <codon-containing exon entry>`,
    add bg_count to every position from bg_start to bg_end-1, indexed by
    its offset to the codon's first base.
    """
    counts = defaultdict(int)
    with open(intersect_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # bg: chr(0) start(1) end(2) count(3)
            # exon: chr(4) start(5) end(6) tid(7) codon_pos(8) strand(9)
            bg_start = int(parts[1])
            bg_end = int(parts[2])
            bg_count = int(parts[3])
            codon_pos = int(parts[8])
            for pos in range(bg_start, bg_end):
                rel = pos - codon_pos
                if x_min <= rel <= x_max:
                    counts[rel] += bg_count
    return counts


# ---------------------------------------------------------------------------
# Stage 8 — plot Figure S5A


def plot_figure_s5a(start_data, stop_data, out_png: Path) -> None:
    samples_order = ["RPF-siLuc", "RPF-siLin28a"]
    labels = {"RPF-siLuc": "siLuc", "RPF-siLin28a": "siLin28a"}

    fig, axes = plt.subplots(
        2, 2,
        figsize=(13, 5.2),
        gridspec_kw={"width_ratios": [10, 7], "wspace": 0.20, "hspace": 0.45},
    )

    for row, sample in enumerate(samples_order):
        # --- start codon panel: -50 .. +50 ---
        ax = axes[row, 0]
        x = np.arange(-50, 51)
        y = np.array([start_data[sample].get(int(p), 0) for p in x]) / 1000.0
        ax.bar(x, y, width=0.7, color="black", linewidth=0)
        ax.axvline(0, color="red", linewidth=1.2)
        ax.set_xlim(-51, 51)
        ax.set_ylim(0, 130)
        ax.set_xticks(np.arange(-50, 51, 10))
        ax.set_yticks([0, 40, 80, 120])
        ax.tick_params(axis="both", direction="out", length=4, labelsize=9)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        if row == 0:
            ax.set_title("start codon", fontsize=11)
        if row == 1:
            ax.set_xlabel(
                "Relative position to start codon of 5'-end of reads",
                fontsize=11,
            )
        ax.set_ylabel(
            f"{labels[sample]}\n\nRaw read count\n(x1000)",
            fontsize=10, rotation=0, ha="right", va="center", labelpad=15,
            fontweight="bold",
        )

        # --- stop codon panel: -50 .. +20 ---
        ax = axes[row, 1]
        x = np.arange(-50, 21)
        y = np.array([stop_data[sample].get(int(p), 0) for p in x]) / 1000.0
        ax.bar(x, y, width=0.7, color="black", linewidth=0)
        ax.axvline(0, color="red", linewidth=1.2)
        ax.set_xlim(-51, 21)
        ax.set_ylim(0, 130)
        ax.set_xticks(np.arange(-50, 21, 10))
        ax.set_yticks([0, 40, 80, 120])
        ax.tick_params(axis="both", direction="out", length=4, labelsize=9)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        if row == 0:
            ax.set_title("stop codon", fontsize=11)

    fig.suptitle(
        "A   Ribosome footprint density near start and stop codons",
        x=0.02, y=0.99, ha="left", fontweight="bold", fontsize=12,
    )
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out_png}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Driver


def main() -> None:
    os.chdir(WORK_DIR)
    print(f"Working in {WORK_DIR}", file=sys.stderr)

    start_gtf = WORK_DIR / "gencode-start.gtf"
    stop_gtf = WORK_DIR / "gencode-stop.gtf"
    plus_exon_gtf = WORK_DIR / "gencode-plusexon.gtf"

    print("[1] Extract start codons (TSL=1, + strand)")
    extract_codon_gtf("start_codon", start_gtf)
    print("[1] Extract stop codons (TSL=1, + strand)")
    extract_codon_gtf("stop_codon", stop_gtf)
    print("[2] Extract + strand exons")
    extract_plus_exon_gtf(plus_exon_gtf)

    start_exon_bed = WORK_DIR / "gencode-exons-containing-startcodon.bed"
    stop_exon_bed = WORK_DIR / "gencode-exons-containing-stopcodon.bed"
    print("[3] Build BED of exons containing a start codon")
    make_codon_exon_bed(start_gtf, plus_exon_gtf, start_exon_bed)
    print("[3] Build BED of exons containing a stop codon")
    make_codon_exon_bed(stop_gtf, plus_exon_gtf, stop_exon_bed)

    start_data = {}
    stop_data = {}
    for sample in SAMPLES:
        print(f"\n=== {sample} ===")
        filt_bam = WORK_DIR / f"filtered-{sample}.bam"
        print(f"[4] Filter BAM (+ strand, length >= 25 nt)")
        filter_bam_plus_long(sample, filt_bam)

        fivep_bg = WORK_DIR / f"fivepcounts-{sample}.bed"
        print(f"[5] 5'-end coverage")
        fivep_bedgraph(filt_bam, fivep_bg)

        start_int = WORK_DIR / f"fivepcounts-startcodon-{sample}.txt"
        stop_int = WORK_DIR / f"fivepcounts-stopcodon-{sample}.txt"
        print(f"[6] Intersect with start-codon exons")
        intersect_fivep_with_codon(fivep_bg, start_exon_bed, start_int)
        print(f"[6] Intersect with stop-codon exons")
        intersect_fivep_with_codon(fivep_bg, stop_exon_bed, stop_int)

        print(f"[7] Aggregate by relative position")
        start_data[sample] = aggregate_relative_counts(start_int, -50, 50)
        stop_data[sample] = aggregate_relative_counts(stop_int, -50, 20)

    print("\n[8] Plot Figure S5A")
    plot_figure_s5a(start_data, stop_data, WORK_DIR / "figs5a_reproduced.png")


if __name__ == "__main__":
    main()
