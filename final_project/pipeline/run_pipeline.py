#!/usr/bin/env python
"""Final project MVP analysis pipeline.

This script implements the amount/position analysis described in
data_analysis_plan.md. It avoids optional dependencies (pysam, statsmodels)
and uses samtools/bedtools plus pandas/numpy/scipy/matplotlib.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ryu/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# ggplot aesthetic (gray panel + white grid, ggplot2-style) with Inter typeface.
plt.style.use("ggplot")
sns.set_context("paper")  # font scaling only; does not override the ggplot style

# Single source of truth for font sizes used across every figure.
FS_BASE = 11   # default text
FS_TITLE = 12  # axes titles
FS_LABEL = 11  # x/y axis labels
FS_TICK = 10   # tick labels
FS_ANNOT = 9   # in-plot text boxes, bar labels, legends, coefficient labels
FS_PANEL = 13  # bold A/B/C/D panel letters

mpl.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        # --- typeface: Inter everywhere, including mathtext ($\Delta$, $R^2$, italics) ---
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Inter",
        "mathtext.it": "Inter:italic",
        "mathtext.bf": "Inter:bold",
        "mathtext.sf": "Inter",
        "mathtext.cal": "Inter:italic",
        "font.size": FS_BASE,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 8.0,
        "axes.labelsize": FS_LABEL,
        "axes.labelpad": 5.0,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_ANNOT,
        "legend.frameon": False,
        # --- all in-plot text black (override ggplot's gray) ---
        "text.color": "black",
        "axes.titlecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.labelcolor": "black",
        "ytick.labelcolor": "black",
        "legend.labelcolor": "black",
    }
)


SAMPLE_COLS = [
    "CLIP-35L33G.bam",
    "RNA-control.bam",
    "RNA-siLin28a.bam",
    "RNA-siLuc.bam",
    "RPF-siLin28a.bam",
    "RPF-siLuc.bam",
]

REGIONS = ["5utr", "cds", "3utr"]
REGION_LABELS = {"5utr": "5' UTR", "cds": "CDS", "3utr": "3' UTR"}
REGION_COLORS = {"5utr": "#4C72B0", "cds": "#DD8452", "3utr": "#55A868"}
MODEL_COLORS = {"amount": "#6B8BA4", "position": "#B56B5B", "baseline": "#9A9A9A"}
ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class Paths:
    project: Path
    pipeline: Path
    data_dir: Path
    read_counts: Path
    gtf: Path
    clip_bam: Path
    rna_control_bam: Path
    results: Path
    work: Path
    tables: Path
    figures: Path
    logs: Path
    samtools: str
    bedtools: str


def log(message: str) -> None:
    print(f"[pipeline] {message}", flush=True)


def fmt_p(p_value: float) -> str:
    if p_value is None or not math.isfinite(float(p_value)):
        return "NA"
    p_value = float(p_value)
    if p_value == 0.0:
        return "<1e-300"
    if p_value < 1e-3:
        return f"{p_value:.1e}"
    return f"{p_value:.3g}"


def p_equals(p_value: float) -> str:
    text = fmt_p(p_value)
    return f"p{text}" if text.startswith("<") else f"p={text}"


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FS_PANEL,
        fontweight="bold",
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    """Save both PNG and PDF versions for report drafting."""
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))


def cleanup_stale_figures(paths: Paths) -> None:
    """Remove figures from earlier design iterations that are no longer produced."""
    stale_stems = [
        "fig1_workflow",
        "fig2_baseline_total_clip_vs_delta_rd",
        "fig3_region_composition_summary",
        "fig4_nested_model_delta_r2",
        "fig5_binding_dose_response",
    ]
    for stem in stale_stems:
        for ext in (".png", ".pdf"):
            stale = paths.figures / f"{stem}{ext}"
            if stale.exists():
                stale.unlink()


def parse_args() -> argparse.Namespace:
    default_project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, default=default_project)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument(
        "--samtools",
        default="/blaze/ryu/conda/envs/lab/bin/samtools",
        help="Path to samtools.",
    )
    parser.add_argument(
        "--bedtools",
        default="/blaze/ryu/conda/envs/lab/bin/bedtools",
        help="Path to bedtools.",
    )
    parser.add_argument(
        "--min-response-count",
        type=int,
        default=30,
        help="Minimum RNA/RPF raw count for response-eligible genes.",
    )
    parser.add_argument(
        "--min-baseline-count",
        type=int,
        default=30,
        help="Minimum raw count in all six samples for the baseline figure.",
    )
    parser.add_argument(
        "--min-clip-points",
        type=int,
        default=10,
        help="Minimum exon-assigned CLIP 5-prime points for composition models.",
    )
    parser.add_argument(
        "--min-rna-region",
        type=int,
        default=5,
        help="Minimum RNA-control reads in each 5UTR/CDS/3UTR region.",
    )
    parser.add_argument(
        "--pseudocount-cpm",
        type=float,
        default=0.1,
        help="CPM pseudocount used in log-ratios.",
    )
    return parser.parse_args()


def init_paths(args: argparse.Namespace) -> Paths:
    project = args.project_dir.resolve()
    data_dir = project / "week1" / "work" / "binfo1-datapack1"
    results = (args.results_dir or (project / "results")).resolve()
    paths = Paths(
        project=project,
        pipeline=project / "pipeline",
        data_dir=data_dir,
        read_counts=project / "week1" / "work" / "read-counts.txt",
        gtf=data_dir / "gencode.gtf",
        clip_bam=data_dir / "CLIP-35L33G.bam",
        rna_control_bam=data_dir / "RNA-control.bam",
        results=results,
        work=results / "work",
        tables=results / "tables",
        figures=results / "figures",
        logs=results / "logs",
        samtools=args.samtools,
        bedtools=args.bedtools,
    )
    for d in [paths.results, paths.work, paths.tables, paths.figures, paths.logs]:
        d.mkdir(parents=True, exist_ok=True)
    return paths


def run_command(cmd: list[str], log_path: Path | None = None) -> None:
    log("running: " + " ".join(map(str, cmd)))
    if log_path is None:
        subprocess.run(cmd, check=True)
        return
    with log_path.open("w") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, text=True)


def capture_int(cmd: list[str]) -> int:
    out = subprocess.check_output(cmd, text=True).strip()
    return int(out)


def strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def parse_attrs(attr_text: str) -> dict[str, object]:
    attrs: dict[str, object] = {}
    tags: list[str] = []
    for key, value in ATTR_RE.findall(attr_text):
        if key == "tag":
            tags.append(value)
        elif key not in attrs:
            attrs[key] = value
    if tags:
        attrs["tag"] = tags
    return attrs


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    items = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not items:
        return []
    merged = [items[0]]
    for start, end in items[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def load_read_counts(paths: Paths) -> pd.DataFrame:
    df = pd.read_csv(paths.read_counts, sep="\t", comment="#")
    df["gene_id"] = df["Geneid"]
    df["gene_id_base"] = df["Geneid"].map(strip_version)
    return df


def compute_gene_response(paths: Paths, args: argparse.Namespace) -> pd.DataFrame:
    log("Step 1: computing gene-level response from read-counts.txt")
    cnts = load_read_counts(paths)
    lib_sizes = cnts[SAMPLE_COLS].sum()
    cpm = cnts[SAMPLE_COLS].div(lib_sizes, axis=1) * 1e6

    pc = args.pseudocount_cpm
    out = cnts[["gene_id", "gene_id_base", "Length"] + SAMPLE_COLS].copy()
    out = out.rename(columns={"Length": "featurecounts_gene_length"})
    for col in SAMPLE_COLS:
        out[col.replace(".bam", "_raw")] = cnts[col].astype(float)
        out[col.replace(".bam", "_cpm")] = cpm[col]

    out["total_CLIP_enrichment_featureCounts"] = np.log2(
        (out["CLIP-35L33G_cpm"] + pc) / (out["RNA-control_cpm"] + pc)
    )
    te_kd = (out["RPF-siLin28a_cpm"] + pc) / (out["RNA-siLin28a_cpm"] + pc)
    te_ctrl = (out["RPF-siLuc_cpm"] + pc) / (out["RNA-siLuc_cpm"] + pc)
    out["delta_rd"] = np.log2(te_kd / te_ctrl)
    out["rna_expression_cpm"] = out[
        ["RNA-control_cpm", "RNA-siLuc_cpm", "RNA-siLin28a_cpm"]
    ].mean(axis=1)
    out["log_rna_expression"] = np.log2(out["rna_expression_cpm"] + pc)
    out["log_featurecounts_gene_length"] = np.log2(
        out["featurecounts_gene_length"].astype(float) + 1.0
    )

    response_cols = [
        "RNA-control.bam",
        "RNA-siLuc.bam",
        "RNA-siLin28a.bam",
        "RPF-siLuc.bam",
        "RPF-siLin28a.bam",
    ]
    out["response_eligible"] = (
        cnts[response_cols].astype(float) >= args.min_response_count
    ).all(axis=1)
    out["baseline_eligible"] = (
        cnts[SAMPLE_COLS].astype(float) >= args.min_baseline_count
    ).all(axis=1)

    lib_sizes.rename("library_size").to_csv(paths.tables / "library_sizes.tsv", sep="\t")
    out.to_csv(paths.tables / "gene_response.tsv", sep="\t", index=False)
    return out


def build_representative_regions(paths: Paths) -> pd.DataFrame:
    reps_path = paths.work / "representative_transcripts.tsv"
    regions_path = paths.work / "regions.all.bed"
    if reps_path.exists() and regions_path.exists():
        log("Step 2-3: reusing representative transcript regions")
        return pd.read_csv(reps_path, sep="\t")

    log("Step 2-3: parsing GTF and building representative transcript regions")
    transcripts: dict[str, dict[str, object]] = {}

    with paths.gtf.open() as fh:
        for line_no, line in enumerate(fh, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attr_text = fields
            if feature not in {"transcript", "exon", "CDS", "UTR"}:
                continue
            attrs = parse_attrs(attr_text)
            tid = attrs.get("transcript_id")
            gid = attrs.get("gene_id")
            if not tid or not gid:
                continue
            tid = str(tid)
            gid = str(gid)
            record = transcripts.setdefault(
                tid,
                {
                    "gene_id": gid,
                    "gene_id_base": strip_version(gid),
                    "gene_name": attrs.get("gene_name", strip_version(gid)),
                    "transcript_id": tid,
                    "chrom": chrom,
                    "strand": strand,
                    "gene_type": attrs.get("gene_type", ""),
                    "transcript_type": attrs.get("transcript_type", ""),
                    "transcript_support_level": attrs.get(
                        "transcript_support_level", "NA"
                    ),
                    "tags": set(),
                    "exon": [],
                    "CDS": [],
                    "UTR": [],
                },
            )
            for key in [
                "gene_type",
                "transcript_type",
                "transcript_support_level",
                "gene_name",
            ]:
                if attrs.get(key):
                    record[key] = attrs[key]
            if attrs.get("tag"):
                record["tags"].update(attrs["tag"])  # type: ignore[index, union-attr]
            if feature in {"exon", "CDS", "UTR"}:
                start = int(start_s) - 1
                end = int(end_s)
                record[feature].append((start, end))  # type: ignore[index, union-attr]
            if line_no % 2_000_000 == 0:
                log(f"  parsed {line_no:,} GTF lines")

    candidates = []
    for tx in transcripts.values():
        if tx["gene_type"] != "protein_coding":
            continue
        if tx["transcript_type"] != "protein_coding":
            continue
        cds_intervals = merge_intervals(tx["CDS"])  # type: ignore[arg-type]
        exon_intervals = merge_intervals(tx["exon"])  # type: ignore[arg-type]
        utr_intervals = merge_intervals(tx["UTR"])  # type: ignore[arg-type]
        cds_len = interval_length(cds_intervals)
        tx_len = interval_length(exon_intervals)
        utr_len = interval_length(utr_intervals)
        if cds_len <= 0 or tx_len <= 0:
            continue
        tags = set(tx["tags"])  # type: ignore[arg-type]
        tsl = str(tx["transcript_support_level"])
        score_tuple = (
            1 if tsl == "1" else 0,
            1 if any(str(t).startswith("appris_principal") for t in tags) else 0,
            1 if "CCDS" in tags else 0,
            1 if "basic" in tags else 0,
            cds_len,
            tx_len,
        )
        tx["cds_length"] = cds_len
        tx["transcript_length"] = tx_len
        tx["utr_length"] = utr_len
        tx["selection_score"] = score_tuple
        candidates.append(tx)

    by_gene: dict[str, dict[str, object]] = {}
    for tx in candidates:
        gid = str(tx["gene_id_base"])
        if gid not in by_gene or tx["selection_score"] > by_gene[gid]["selection_score"]:
            by_gene[gid] = tx

    reps = pd.DataFrame(
        [
            {
                "gene_id": tx["gene_id"],
                "gene_id_base": tx["gene_id_base"],
                "gene_name": tx["gene_name"],
                "transcript_id": tx["transcript_id"],
                "chrom": tx["chrom"],
                "strand": tx["strand"],
                "transcript_length": tx["transcript_length"],
                "cds_length": tx["cds_length"],
                "utr_length": tx["utr_length"],
                "selection_score": "|".join(map(str, tx["selection_score"])),
            }
            for tx in by_gene.values()
        ]
    ).sort_values(["chrom", "gene_id_base"])
    reps.to_csv(paths.work / "representative_transcripts.tsv", sep="\t", index=False)

    region_records = []
    for tx in by_gene.values():
        cds_intervals = merge_intervals(tx["CDS"])  # type: ignore[arg-type]
        utr_intervals = merge_intervals(tx["UTR"])  # type: ignore[arg-type]
        if not cds_intervals:
            continue
        cds_min = min(start for start, _ in cds_intervals)
        cds_max = max(end for _, end in cds_intervals)
        base = {
            "chrom": tx["chrom"],
            "gene_id_base": tx["gene_id_base"],
            "score": 0,
            "strand": tx["strand"],
            "transcript_id": tx["transcript_id"],
            "gene_name": tx["gene_name"],
        }
        for start, end in cds_intervals:
            region_records.append(
                {**base, "start": start, "end": end, "region": "cds"}
            )
        for start, end in utr_intervals:
            region = None
            if tx["strand"] == "+":
                if end <= cds_min:
                    region = "5utr"
                elif start >= cds_max:
                    region = "3utr"
            else:
                if start >= cds_max:
                    region = "5utr"
                elif end <= cds_min:
                    region = "3utr"
            if region is None:
                continue
            region_records.append({**base, "start": start, "end": end, "region": region})

    regions = pd.DataFrame(region_records)
    regions["region_length"] = regions["end"] - regions["start"]
    regions = regions.sort_values(["chrom", "start", "end", "strand", "gene_id_base"])

    bed_cols = [
        "chrom",
        "start",
        "end",
        "gene_id_base",
        "score",
        "strand",
        "transcript_id",
        "region",
        "region_length",
        "gene_name",
    ]
    regions[bed_cols].to_csv(
        paths.work / "regions.all.bed", sep="\t", header=False, index=False
    )
    for region in REGIONS:
        regions.loc[regions["region"] == region, bed_cols].to_csv(
            paths.work / f"regions.{region}.bed", sep="\t", header=False, index=False
        )
    return reps


def find_gene_locus(gtf: Path, gene_name: str) -> dict[str, object] | None:
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("gene_name") == gene_name:
                return {
                    "chrom": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene_id": attrs.get("gene_id"),
                    "gene_name": gene_name,
                }
    return None


def check_strandedness(paths: Paths) -> pd.DataFrame:
    log("Checking strandedness direction with Actb locus")
    locus = find_gene_locus(paths.gtf, "Actb")
    if locus is None:
        out = pd.DataFrame([{"gene_name": "Actb", "status": "not_found"}])
        out.to_csv(paths.tables / "strandedness_check.tsv", sep="\t", index=False)
        return out

    region = f"{locus['chrom']}:{locus['start']}-{locus['end']}"
    cmd = [paths.samtools, "view", "-F", "2308", str(paths.clip_bam), region]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    fwd = 0
    rev = 0
    total = 0
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.split("\t", 3)
        if len(fields) < 2:
            continue
        flag = int(fields[1])
        total += 1
        if flag & 16:
            rev += 1
        else:
            fwd += 1
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("samtools view failed during strandedness check")

    same = rev if locus["strand"] == "-" else fwd
    opposite = fwd if locus["strand"] == "-" else rev
    same_fraction = same / total if total else float("nan")
    out = pd.DataFrame(
        [
            {
                **locus,
                "query_region": region,
                "forward_reads": fwd,
                "reverse_reads": rev,
                "same_strand_reads": same,
                "opposite_strand_reads": opposite,
                "same_strand_fraction": same_fraction,
                "inferred_featureCounts_strand": "-s 1" if same_fraction >= 0.8 else "check",
            }
        ]
    )
    out.to_csv(paths.tables / "strandedness_check.tsv", sep="\t", index=False)
    return out


def cigar_reference_length(cigar: str) -> int:
    length = 0
    for size_s, op in CIGAR_RE.findall(cigar):
        if op in {"M", "D", "N", "=", "X"}:
            length += int(size_s)
    return length


def extract_clip_5end_points(paths: Paths) -> dict[str, int]:
    out_bed = paths.work / "clip_5end_points.bed"
    stats_path = paths.tables / "clip_point_extraction_stats.tsv"
    if out_bed.exists() and stats_path.exists():
        log("Step 4: reusing existing CLIP 5-prime point BED")
        return pd.read_csv(stats_path, sep="\t").iloc[0].to_dict()

    log("Step 4: extracting NH=1 CLIP 5-prime-end points")
    cmd = [paths.samtools, "view", "-@", "4", "-F", "2308", str(paths.clip_bam)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
    assert proc.stdout is not None

    total = 0
    nh1 = 0
    skipped_no_nh = 0
    skipped_bad_cigar = 0
    with out_bed.open("w") as out:
        for line in proc.stdout:
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            nh = None
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    nh = int(tag.rsplit(":", 1)[-1])
                    break
            if nh is None:
                skipped_no_nh += 1
                continue
            if nh != 1:
                continue
            chrom = fields[2]
            pos0 = int(fields[3]) - 1
            cigar = fields[5]
            ref_len = cigar_reference_length(cigar)
            if ref_len <= 0:
                skipped_bad_cigar += 1
                continue
            flag = int(fields[1])
            reverse = bool(flag & 16)
            if reverse:
                start = pos0 + ref_len - 1
                strand = "-"
            else:
                start = pos0
                strand = "+"
            end = start + 1
            out.write(f"{chrom}\t{start}\t{end}\t{fields[0]}\t0\t{strand}\n")
            nh1 += 1
            if total % 5_000_000 == 0:
                log(f"  processed {total:,} CLIP alignments; kept {nh1:,} NH=1 points")
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("samtools view failed during CLIP point extraction")

    stats_df = pd.DataFrame(
        [
            {
                "sam_alignments_seen": total,
                "nh1_points_written": nh1,
                "skipped_no_nh": skipped_no_nh,
                "skipped_bad_cigar": skipped_bad_cigar,
            }
        ]
    )
    stats_df.to_csv(stats_path, sep="\t", index=False)
    return stats_df.iloc[0].to_dict()


def bedtools_count_regions(paths: Paths) -> None:
    log("Step 4: counting CLIP points and RNA-control reads by region")
    regions_bed = paths.work / "regions.all.bed"
    clip_points = paths.work / "clip_5end_points.bed"
    clip_out = paths.work / "regions.clip_point_counts.bed"
    rna_out = paths.work / "regions.rna_control_counts.bed"

    if not clip_out.exists():
        with clip_out.open("w") as out, (paths.logs / "bedtools_clip_counts.err").open(
            "w"
        ) as err:
            subprocess.run(
                [
                    paths.bedtools,
                    "coverage",
                    "-counts",
                    "-s",
                    "-a",
                    str(regions_bed),
                    "-b",
                    str(clip_points),
                ],
                check=True,
                stdout=out,
                stderr=err,
                text=True,
            )

    if not rna_out.exists():
        with rna_out.open("w") as out, (paths.logs / "bedtools_rna_counts.err").open(
            "w"
        ) as err:
            subprocess.run(
                [
                    paths.bedtools,
                    "coverage",
                    "-counts",
                    "-s",
                    "-split",
                    "-a",
                    str(regions_bed),
                    "-b",
                    str(paths.rna_control_bam),
                ],
                check=True,
                stdout=out,
                stderr=err,
                text=True,
            )


def aggregate_region_counts(paths: Paths) -> pd.DataFrame:
    log("Step 5: aggregating region-level counts")
    cols = [
        "chrom",
        "start",
        "end",
        "gene_id_base",
        "score",
        "strand",
        "transcript_id",
        "region",
        "region_length",
        "gene_name",
        "count",
    ]
    clip = pd.read_csv(paths.work / "regions.clip_point_counts.bed", sep="\t", names=cols)
    rna = pd.read_csv(paths.work / "regions.rna_control_counts.bed", sep="\t", names=cols)

    key_cols = ["gene_id_base", "gene_name", "transcript_id", "region"]
    clip_g = (
        clip.groupby(key_cols, as_index=False)
        .agg(clip_points=("count", "sum"), region_length=("region_length", "sum"))
    )
    rna_g = (
        rna.groupby(key_cols, as_index=False)
        .agg(rna_control_reads=("count", "sum"), region_length_rna=("region_length", "sum"))
    )
    long = clip_g.merge(rna_g, on=key_cols, how="outer")
    long["clip_points"] = long["clip_points"].fillna(0).astype(int)
    long["rna_control_reads"] = long["rna_control_reads"].fillna(0).astype(int)
    long["region_length"] = long["region_length"].fillna(long["region_length_rna"]).fillna(0)
    long = long.drop(columns=["region_length_rna"])
    long.to_csv(paths.work / "region_counts.long.tsv", sep="\t", index=False)

    wide_parts = []
    for value in ["clip_points", "rna_control_reads", "region_length"]:
        pivot = long.pivot_table(
            index=["gene_id_base", "gene_name", "transcript_id"],
            columns="region",
            values=value,
            aggfunc="sum",
            fill_value=0,
        )
        pivot.columns = [f"{value}_{c}" for c in pivot.columns]
        wide_parts.append(pivot)
    wide = pd.concat(wide_parts, axis=1).reset_index()
    for region in REGIONS:
        for value in ["clip_points", "rna_control_reads", "region_length"]:
            col = f"{value}_{region}"
            if col not in wide.columns:
                wide[col] = 0
    wide.to_csv(paths.work / "region_counts.tsv", sep="\t", index=False)
    return wide


def make_analysis_table(
    paths: Paths,
    gene_response: pd.DataFrame,
    region_counts: pd.DataFrame,
    args: argparse.Namespace,
    clip_stats: dict[str, int],
) -> pd.DataFrame:
    log("Step 5: computing adjusted composition and model-ready table")
    reps_path = paths.work / "representative_transcripts.tsv"
    if reps_path.exists():
        reps = pd.read_csv(reps_path, sep="\t")
        rep_cols = [
            "gene_id_base",
            "chrom",
            "strand",
            "transcript_length",
            "cds_length",
            "utr_length",
        ]
        region_counts = region_counts.merge(
            reps[[c for c in rep_cols if c in reps.columns]],
            on="gene_id_base",
            how="left",
        )
    table = gene_response.merge(region_counts, on="gene_id_base", how="inner")

    clip_lib = int(clip_stats["nh1_points_written"])
    rna_lib = capture_int(
        [paths.samtools, "view", "-@", "4", "-c", "-F", "2308", str(paths.rna_control_bam)]
    )
    pc = args.pseudocount_cpm

    for region in REGIONS:
        table[f"clip_points_{region}_cpm"] = (
            table[f"clip_points_{region}"].astype(float) / clip_lib * 1e6
        )
        table[f"rna_control_reads_{region}_cpm"] = (
            table[f"rna_control_reads_{region}"].astype(float) / rna_lib * 1e6
        )
        table[f"E_{region}"] = (
            table[f"clip_points_{region}_cpm"] + pc
        ) / (table[f"rna_control_reads_{region}_cpm"] + pc)

    table["clip_points_total"] = table[[f"clip_points_{r}" for r in REGIONS]].sum(axis=1)
    table["rna_control_total"] = table[[f"rna_control_reads_{r}" for r in REGIONS]].sum(
        axis=1
    )
    table["clip_points_total_cpm"] = table["clip_points_total"].astype(float) / clip_lib * 1e6
    table["rna_control_total_cpm"] = table["rna_control_total"].astype(float) / rna_lib * 1e6
    table["total_CLIP_amount_point"] = np.log2(
        (table["clip_points_total_cpm"] + pc)
        / (table["rna_control_total_cpm"] + pc)
    )

    e_sum = table[[f"E_{r}" for r in REGIONS]].sum(axis=1)
    for region in REGIONS:
        table[f"f_{region}"] = table[f"E_{region}"] / e_sum
        table[f"raw_f_{region}"] = table[f"clip_points_{region}"] / table[
            "clip_points_total"
        ].replace(0, np.nan)

    table["has_all_regions"] = (
        table[[f"region_length_{r}" for r in REGIONS]].astype(float) > 0
    ).all(axis=1)
    table["composition_eligible"] = (
        table["response_eligible"]
        & table["has_all_regions"]
        & (table["clip_points_total"] >= args.min_clip_points)
        & (
            table[[f"rna_control_reads_{r}" for r in REGIONS]].astype(float)
            >= args.min_rna_region
        ).all(axis=1)
        & np.isfinite(table["delta_rd"])
        & np.isfinite(table["total_CLIP_amount_point"])
        & np.isfinite(table["f_cds"])
        & np.isfinite(table["f_3utr"])
    )

    if "transcript_length" not in table.columns:
        table["transcript_length"] = table[[f"region_length_{r}" for r in REGIONS]].sum(
            axis=1
        )
    table["log_transcript_length"] = np.log2(table["transcript_length"].astype(float) + 1.0)
    table.to_csv(paths.work / "analysis_table.tsv", sep="\t", index=False)

    eligibility = pd.DataFrame(
        [
            {"stage": "gene_response_rows", "n": len(gene_response)},
            {"stage": "response_eligible", "n": int(gene_response["response_eligible"].sum())},
            {"stage": "representative_region_matched", "n": len(table)},
            {"stage": "has_all_regions", "n": int(table["has_all_regions"].sum())},
            {
                "stage": f"clip_points_total_ge_{args.min_clip_points}",
                "n": int((table["clip_points_total"] >= args.min_clip_points).sum()),
            },
            {
                "stage": f"rna_each_region_ge_{args.min_rna_region}",
                "n": int(
                    (
                        table[[f"rna_control_reads_{r}" for r in REGIONS]].astype(float)
                        >= args.min_rna_region
                    )
                    .all(axis=1)
                    .sum()
                ),
            },
            {"stage": "composition_eligible_model_set", "n": int(table["composition_eligible"].sum())},
        ]
    )
    eligibility.to_csv(paths.tables / "eligibility_summary.tsv", sep="\t", index=False)

    lib = pd.DataFrame(
        [
            {"library": "CLIP_NH1_5end_points", "size": clip_lib},
            {"library": "RNA_control_primary_mapped_reads", "size": rna_lib},
        ]
    )
    lib.to_csv(paths.tables / "point_signal_library_sizes.tsv", sep="\t", index=False)
    return table


def standardize_columns(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    vals = df[cols].astype(float).to_numpy()
    means = np.nanmean(vals, axis=0)
    sds = np.nanstd(vals, axis=0)
    sds[sds == 0] = 1.0
    return (vals - means) / sds


def ols_fit(y: np.ndarray, x: np.ndarray) -> dict[str, object]:
    x = np.column_stack([np.ones(len(y)), x])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    yhat = x @ beta
    resid = y - yhat
    sse = float(np.sum(resid**2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else float("nan")
    return {"beta": beta, "yhat": yhat, "sse": sse, "r2": r2, "n": len(y), "p": x.shape[1]}


def nested_f_test(reduced: dict[str, object], full: dict[str, object]) -> dict[str, float]:
    sse_r = float(reduced["sse"])
    sse_f = float(full["sse"])
    p_r = int(reduced["p"])
    p_f = int(full["p"])
    n = int(full["n"])
    df_num = p_f - p_r
    df_den = n - p_f
    if df_num <= 0 or df_den <= 0 or sse_f <= 0:
        return {"F": float("nan"), "p": float("nan"), "df_num": df_num, "df_den": df_den}
    numerator = max(sse_r - sse_f, 0.0) / df_num
    denominator = sse_f / df_den
    f_stat = numerator / denominator if denominator > 0 else float("nan")
    p_value = float(stats.f.sf(f_stat, df_num, df_den)) if math.isfinite(f_stat) else float("nan")
    return {"F": f_stat, "p": p_value, "df_num": df_num, "df_den": df_den}


def fit_models(paths: Paths, table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    log("Step 6: fitting nested OLS models on identical composition-eligible gene set")
    model_df = table.loc[table["composition_eligible"]].copy()
    model_cols = [
        "delta_rd",
        "log_rna_expression",
        "log_transcript_length",
        "total_CLIP_amount_point",
        "f_cds",
        "f_3utr",
    ]
    model_df = model_df.replace([np.inf, -np.inf], np.nan).dropna(subset=model_cols)
    n_model = len(model_df)
    if n_model < 20:
        raise RuntimeError(f"Too few composition-eligible genes for modeling: n={n_model}")

    y = model_df["delta_rd"].astype(float).to_numpy()
    y_sd = float(np.nanstd(y))
    x0_cols = ["log_rna_expression", "log_transcript_length"]
    x1_cols = x0_cols + ["total_CLIP_amount_point"]
    x2_cols = x1_cols + ["f_cds", "f_3utr"]

    fit0 = ols_fit(y, standardize_columns(model_df, x0_cols))
    fit1 = ols_fit(y, standardize_columns(model_df, x1_cols))
    fit2 = ols_fit(y, standardize_columns(model_df, x2_cols))

    assert fit0["n"] == fit1["n"] == fit2["n"]
    f01 = nested_f_test(fit0, fit1)
    f12 = nested_f_test(fit1, fit2)

    model_df["pred_M1"] = fit1["yhat"]
    model_df["pred_M2"] = fit2["yhat"]
    model_df.to_csv(paths.tables / "model_gene_set.tsv", sep="\t", index=False)

    summary_rows = [
        {
            "model": "M0",
            "predictors": "log_rna_expression + log_transcript_length",
            "n": fit0["n"],
            "p_parameters": fit0["p"],
            "r2": fit0["r2"],
            "sse": fit0["sse"],
        },
        {
            "model": "M1",
            "predictors": "M0 + total_CLIP_amount_point",
            "n": fit1["n"],
            "p_parameters": fit1["p"],
            "r2": fit1["r2"],
            "sse": fit1["sse"],
            "delta_r2_vs_previous": float(fit1["r2"]) - float(fit0["r2"]),
            "nested_F_vs_previous": f01["F"],
            "nested_p_vs_previous": f01["p"],
            "df_num_vs_previous": f01["df_num"],
            "df_den_vs_previous": f01["df_den"],
        },
        {
            "model": "M2",
            "predictors": "M1 + f_cds + f_3utr",
            "n": fit2["n"],
            "p_parameters": fit2["p"],
            "r2": fit2["r2"],
            "sse": fit2["sse"],
            "delta_r2_vs_previous": float(fit2["r2"]) - float(fit1["r2"]),
            "nested_F_vs_previous": f12["F"],
            "nested_p_vs_previous": f12["p"],
            "df_num_vs_previous": f12["df_num"],
            "df_den_vs_previous": f12["df_den"],
        },
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(paths.tables / "model_summary.tsv", sep="\t", index=False)

    coef_rows = []
    for model_name, fit, cols in [
        ("M0", fit0, x0_cols),
        ("M1", fit1, x1_cols),
        ("M2", fit2, x2_cols),
    ]:
        terms = ["intercept"] + cols
        for term, beta in zip(terms, fit["beta"]):
            raw_effect = float(beta)
            standardized_beta = (
                raw_effect / y_sd
                if term != "intercept" and y_sd > 0 and math.isfinite(y_sd)
                else np.nan
            )
            coef_rows.append(
                {
                    "model": model_name,
                    "term": term,
                    "effect_per_1sd_predictor_delta_rd": raw_effect,
                    "standardized_beta": standardized_beta,
                    "n": fit["n"],
                }
            )
    pd.DataFrame(coef_rows).to_csv(paths.tables / "model_coefficients.tsv", sep="\t", index=False)

    assoc_rows = []
    for col in [
        "total_CLIP_enrichment_featureCounts",
        "total_CLIP_amount_point",
        "clip_points_total",
        "f_5utr",
        "f_cds",
        "f_3utr",
    ]:
        tmp = model_df[["delta_rd", col]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(tmp) >= 3:
            pear = stats.pearsonr(tmp[col], tmp["delta_rd"])
            spear = stats.spearmanr(tmp[col], tmp["delta_rd"])
            assoc_rows.append(
                {
                    "variable": col,
                    "n": len(tmp),
                    "pearson_r": pear.statistic,
                    "pearson_p": pear.pvalue,
                    "spearman_rho": spear.statistic,
                    "spearman_p": spear.pvalue,
                }
            )
    assoc = pd.DataFrame(assoc_rows)
    assoc.to_csv(paths.tables / "association_summary.tsv", sep="\t", index=False)
    return summary, model_df


def plot_baseline(paths: Paths, gene_response: pd.DataFrame) -> dict[str, float]:
    df = gene_response.loc[gene_response["baseline_eligible"]].replace(
        [np.inf, -np.inf], np.nan
    )
    df = df.dropna(subset=["total_CLIP_enrichment_featureCounts", "delta_rd"])
    pear = stats.pearsonr(df["total_CLIP_enrichment_featureCounts"], df["delta_rd"])
    spear = stats.spearmanr(df["total_CLIP_enrichment_featureCounts"], df["delta_rd"])

    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    ax.scatter(
        df["total_CLIP_enrichment_featureCounts"],
        df["delta_rd"],
        s=6,
        alpha=0.18,
        color="#333333",
        edgecolors="none",
        rasterized=True,
    )
    sns.regplot(
        data=df,
        x="total_CLIP_enrichment_featureCounts",
        y="delta_rd",
        scatter=False,
        color="#C0392B",
        ax=ax,
        line_kws={"lw": 1.6},
    )
    ax.set_xlabel("LIN28A CLIP enrichment ($\\log_2$)")
    ax.set_ylabel("$\\Delta$ ribosome density ($\\log_2$)")
    ax.text(
        0.03,
        0.97,
        f"n={len(df):,}\nPearson r={pear.statistic:.3f}, {p_equals(pear.pvalue)}\nSpearman rho={spear.statistic:.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=FS_ANNOT,
        bbox=dict(facecolor="white", edgecolor="#dddddd", alpha=0.85),
    )
    fig.tight_layout()
    save_figure(fig, paths.figures / "fig1_baseline_total_clip_vs_delta_rd.png")
    plt.close(fig)
    return {
        "baseline_n": len(df),
        "baseline_pearson_r": pear.statistic,
        "baseline_pearson_p": pear.pvalue,
        "baseline_spearman_rho": spear.statistic,
        "baseline_spearman_p": spear.pvalue,
    }


def plot_composition(paths: Paths, table: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5))
    length_cols = [f"region_length_{r}" for r in REGIONS]
    length_long = table[["gene_id_base"] + length_cols].melt(
        id_vars="gene_id_base", var_name="region", value_name="length"
    )
    length_long["region"] = length_long["region"].str.replace("region_length_", "", regex=False)
    length_long["region_label"] = length_long["region"].map(REGION_LABELS)
    length_long = length_long.loc[length_long["length"] > 0]
    length_long["log10_length"] = np.log10(length_long["length"])
    sns.boxplot(
        data=length_long,
        x="region_label",
        y="log10_length",
        order=[REGION_LABELS[r] for r in REGIONS],
        hue="region_label",
        palette={REGION_LABELS[r]: REGION_COLORS[r] for r in REGIONS},
        legend=False,
        ax=axes[0, 0],
        fliersize=1,
    )
    add_panel_label(axes[0, 0], "A")
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("$\\log_{10}$ length")

    comp = table.loc[table["composition_eligible"], ["gene_id_base"] + [f"f_{r}" for r in REGIONS]]
    comp_long = comp.melt(id_vars="gene_id_base", var_name="region", value_name="fraction")
    comp_long["region"] = comp_long["region"].str.replace("f_", "", regex=False)
    comp_long["region_label"] = comp_long["region"].map(REGION_LABELS)
    sns.violinplot(
        data=comp_long,
        x="region_label",
        y="fraction",
        order=[REGION_LABELS[r] for r in REGIONS],
        hue="region_label",
        palette={REGION_LABELS[r]: REGION_COLORS[r] for r in REGIONS},
        legend=False,
        inner="quartile",
        cut=0,
        ax=axes[0, 1],
    )
    add_panel_label(axes[0, 1], "B")
    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel("fraction of adjusted enrichment")

    raw_adj_rows = []
    for region in REGIONS:
        raw_adj_rows.append(
            {
                "region": region,
                "region_label": REGION_LABELS[region],
                "type": "Raw points",
                "mean": table.loc[table["composition_eligible"], f"raw_f_{region}"].mean(),
            }
        )
        raw_adj_rows.append(
            {
                "region": region,
                "region_label": REGION_LABELS[region],
                "type": "RNA-adjusted",
                "mean": table.loc[table["composition_eligible"], f"f_{region}"].mean(),
            }
        )
    raw_adj = pd.DataFrame(raw_adj_rows)
    sns.barplot(
        data=raw_adj,
        x="type",
        y="mean",
        hue="region_label",
        hue_order=[REGION_LABELS[r] for r in REGIONS],
        palette={REGION_LABELS[r]: REGION_COLORS[r] for r in REGIONS},
        ax=axes[1, 0],
    )
    add_panel_label(axes[1, 0], "C")
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("mean fraction")
    axes[1, 0].legend(
        title="",
        fontsize=FS_ANNOT,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=3,
        borderaxespad=0.0,
        columnspacing=0.9,
        handlelength=1.2,
        handletextpad=0.4,
    )

    eligibility = pd.read_csv(paths.tables / "eligibility_summary.tsv", sep="\t")
    def stage_label(stage: str) -> str:
        labels = {
            "gene_response_rows": "All count-table genes",
            "response_eligible": "Response-eligible genes",
            "representative_region_matched": "Representative transcript matched",
            "has_all_regions": "5'UTR/CDS/3'UTR annotated",
            "composition_eligible_model_set": "Final model set",
        }
        if stage in labels:
            return labels[stage]
        if stage.startswith("clip_points_total_ge_"):
            return f">={stage.rsplit('_', 1)[-1]} total CLIP points"
        if stage.startswith("rna_each_region_ge_"):
            return f"RNA-control >={stage.rsplit('_', 1)[-1]} in each region"
        return stage.replace("_", " ")

    eligibility["stage_label"] = eligibility["stage"].map(stage_label)
    sns.barplot(data=eligibility, y="stage_label", x="n", color="#B8B8B8", ax=axes[1, 1])
    add_panel_label(axes[1, 1], "D")
    axes[1, 1].set_xlabel("genes")
    axes[1, 1].set_ylabel("")
    axes[1, 1].set_xlim(0, eligibility["n"].max() * 1.15)
    for container in axes[1, 1].containers:
        axes[1, 1].bar_label(container, fmt="%.0f", fontsize=FS_ANNOT, padding=2)

    fig.tight_layout()
    save_figure(fig, paths.figures / "fig2_region_composition_summary.png")
    plt.close(fig)


def plot_models(paths: Paths, summary: pd.DataFrame, model_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(11.2, 3.6),
        gridspec_kw={"width_ratios": [1.15, 1.05, 1.2]},
    )

    model_palette = {
        "M0": MODEL_COLORS["baseline"],
        "M1": MODEL_COLORS["amount"],
        "M2": MODEL_COLORS["position"],
    }
    sns.barplot(data=summary, x="model", y="r2", hue="model", palette=model_palette, legend=False, ax=axes[0])
    add_panel_label(axes[0], "A")
    axes[0].set_ylabel("variance explained (R²)")
    axes[0].set_xlabel("")
    axes[0].set_ylim(0, max(summary["r2"]) * 1.18)
    for container in axes[0].containers:
        axes[0].bar_label(container, fmt="%.3f", fontsize=FS_ANNOT, padding=2)
    model_tick_labels = {
        "M0": "M0\ncovariates",
        "M1": "M1\n+ amount",
        "M2": "M2\n+ position",
    }
    axes[0].set_xticks(axes[0].get_xticks())
    axes[0].set_xticklabels(
        [model_tick_labels.get(t.get_text(), t.get_text()) for t in axes[0].get_xticklabels()]
    )

    delta = summary.loc[summary["model"].isin(["M1", "M2"])].copy()
    delta["term_added"] = ["amount", "position"]
    sns.barplot(
        data=delta,
        x="term_added",
        y="delta_r2_vs_previous",
        hue="term_added",
        palette={"amount": MODEL_COLORS["amount"], "position": MODEL_COLORS["position"]},
        legend=False,
        ax=axes[1],
    )
    add_panel_label(axes[1], "B")
    axes[1].set_ylabel("added variance explained (ΔR²)")
    axes[1].set_xlabel("")
    axes[1].set_ylim(0, max(delta["delta_r2_vs_previous"]) * 1.18)
    for container in axes[1].containers:
        axes[1].bar_label(container, fmt="%.4f", fontsize=FS_ANNOT, padding=2)
    added_tick_labels = {
        "amount": "binding\namount",
        "position": "region\ncomposition",
    }
    axes[1].set_xticks(axes[1].get_xticks())
    axes[1].set_xticklabels(
        [added_tick_labels.get(t.get_text(), t.get_text()) for t in axes[1].get_xticklabels()]
    )
    m2 = summary.loc[summary["model"] == "M2"].iloc[0]
    axes[1].text(
        0.98,
        0.96,
        f"M1 vs M2\nF={m2['nested_F_vs_previous']:.2f}\n{p_equals(m2['nested_p_vs_previous'])}",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=FS_ANNOT,
        bbox=dict(facecolor="white", edgecolor="#dddddd", alpha=0.85),
    )

    coefs = pd.read_csv(paths.tables / "model_coefficients.tsv", sep="\t")
    m2c = coefs.loc[coefs["model"] == "M2"].set_index("term")["standardized_beta"]
    coef_terms = [
        ("total_CLIP_amount_point", "total CLIP amount", MODEL_COLORS["amount"]),
        ("f_cds", "CDS fraction", REGION_COLORS["cds"]),
        ("f_3utr", "3' UTR fraction", REGION_COLORS["3utr"]),
    ]
    labels = [lab for _, lab, _ in coef_terms]
    betas = [float(m2c[key]) for key, _, _ in coef_terms]
    colors = [c for _, _, c in coef_terms]
    ypos = np.arange(len(labels))[::-1]
    axes[2].barh(ypos, betas, color=colors, edgecolor="none")
    axes[2].axvline(0, color="#555555", lw=0.9)
    axes[2].set_yticks(ypos)
    axes[2].set_yticklabels(labels)
    xmax = max(abs(b) for b in betas) * 1.45
    axes[2].set_xlim(-xmax, xmax)
    for y, b in zip(ypos, betas):
        axes[2].text(
            b + xmax * 0.03 * (1 if b >= 0 else -1),
            y,
            f"{b:+.3f}",
            va="center",
            ha="left" if b >= 0 else "right",
            fontsize=FS_ANNOT,
        )
    add_panel_label(axes[2], "C")
    axes[2].set_xlabel("standardized $\\beta$")

    fig.tight_layout()
    save_figure(fig, paths.figures / "fig3_nested_model.png")
    plt.close(fig)


def plot_dose(paths: Paths, table: pd.DataFrame) -> pd.DataFrame:
    df = table.loc[table["composition_eligible"]].copy()
    df = df.loc[df["clip_points_total"] > 0]
    df["dose_group"] = pd.qcut(
        df["clip_points_total"],
        q=[0, 0.25, 0.5, 0.75, 1.0],
        labels=["Q1", "Q2", "Q3", "Q4"],
        duplicates="drop",
    )
    dose_summary = (
        df.groupby("dose_group", observed=True)
        .agg(
            n=("delta_rd", "size"),
            median_delta_rd=("delta_rd", "median"),
            mean_delta_rd=("delta_rd", "mean"),
            median_clip_points=("clip_points_total", "median"),
        )
        .reset_index()
    )
    rho = stats.spearmanr(df["clip_points_total"], df["delta_rd"])
    dose_summary["overall_spearman_rho"] = rho.statistic
    dose_summary["overall_spearman_p"] = rho.pvalue
    dose_summary.to_csv(paths.tables / "dose_summary.tsv", sep="\t", index=False)

    order = [str(x) for x in dose_summary["dose_group"]]
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    sns.violinplot(
        data=df,
        x="dose_group",
        y="delta_rd",
        order=order,
        inner=None,
        cut=0,
        color=MODEL_COLORS["amount"],
        ax=ax,
    )
    sns.boxplot(
        data=df,
        x="dose_group",
        y="delta_rd",
        order=order,
        width=0.25,
        fliersize=0,
        color="white",
        ax=ax,
    )
    ax.plot(
        np.arange(len(dose_summary)),
        dose_summary["median_delta_rd"],
        color=MODEL_COLORS["position"],
        marker="D",
        markersize=4,
        linewidth=1.4,
        label="median",
    )
    ax.axhline(0, color="#888888", lw=0.8)
    ax.set_xlabel("LIN28A CLIP point-count quartile", fontsize=FS_TICK)
    ax.set_ylabel("$\\Delta$ ribosome density")
    ax.set_xticks(np.arange(len(dose_summary)))
    ax.set_xticklabels(
        [f"{row['dose_group']}\nn={int(row['n']):,}" for _, row in dose_summary.iterrows()]
    )
    ax.text(
        0.03,
        0.97,
        f"n={len(df):,}\nSpearman rho={rho.statistic:.3f}, {p_equals(rho.pvalue)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=FS_ANNOT,
        bbox=dict(facecolor="white", edgecolor="#dddddd", alpha=0.85),
    )
    ax.legend(loc="upper right", fontsize=FS_ANNOT, handlelength=1.2, borderpad=0.3)
    fig.tight_layout()
    save_figure(fig, paths.figures / "fig4_dose_response.png")
    plt.close(fig)
    return dose_summary


def write_run_manifest(paths: Paths, args: argparse.Namespace) -> None:
    manifest = pd.DataFrame(
        [
            {"key": "project_dir", "value": str(paths.project)},
            {"key": "results_dir", "value": str(paths.results)},
            {"key": "read_counts", "value": str(paths.read_counts)},
            {"key": "gtf", "value": str(paths.gtf)},
            {"key": "clip_bam", "value": str(paths.clip_bam)},
            {"key": "rna_control_bam", "value": str(paths.rna_control_bam)},
            {"key": "min_response_count", "value": args.min_response_count},
            {"key": "min_baseline_count", "value": args.min_baseline_count},
            {"key": "min_clip_points", "value": args.min_clip_points},
            {"key": "min_rna_region", "value": args.min_rna_region},
            {"key": "pseudocount_cpm", "value": args.pseudocount_cpm},
        ]
    )
    manifest.to_csv(paths.results / "run_manifest.tsv", sep="\t", index=False)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "(empty table)"
    view = df.copy()
    view = view.replace({np.nan: ""})
    headers = [str(c) for c in view.columns]
    rows = [[str(v) for v in row] for row in view.to_numpy()]
    widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows)) if rows else len(headers[i])
        for i in range(len(headers))
    ]
    header = "| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |"
    body = [
        "| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |"
        for row in rows
    ]
    return "\n".join([header, sep] + body)


def write_summary(
    paths: Paths,
    baseline_stats: dict[str, float],
    model_summary: pd.DataFrame,
    dose_summary: pd.DataFrame,
) -> None:
    eligibility = pd.read_csv(paths.tables / "eligibility_summary.tsv", sep="\t")
    stranded = pd.read_csv(paths.tables / "strandedness_check.tsv", sep="\t")
    clip_stats = pd.read_csv(paths.tables / "clip_point_extraction_stats.tsv", sep="\t")
    assoc = pd.read_csv(paths.tables / "association_summary.tsv", sep="\t")
    coeffs = pd.read_csv(paths.tables / "model_coefficients.tsv", sep="\t")

    m0 = model_summary.loc[model_summary["model"] == "M0"].iloc[0]
    m1 = model_summary.loc[model_summary["model"] == "M1"].iloc[0]
    m2 = model_summary.loc[model_summary["model"] == "M2"].iloc[0]
    position_sig = bool(m2["nested_p_vs_previous"] < 0.05)
    amount_sig = bool(m1["nested_p_vs_previous"] < 0.05)

    delta_r2_position = float(m2["delta_r2_vs_previous"])
    if position_sig and delta_r2_position < 0.02:
        position_text = (
            "position composition은 통계적으로 유의하지만 증분 설명력은 작다. 따라서 주효과는 total amount이고, 위치는 보조적인 신호로 해석하는 것이 안전하다."
        )
    elif position_sig:
        position_text = (
            "position composition이 total CLIP amount 이후에도 유의한 추가 설명력을 보였다."
        )
    else:
        position_text = (
            "position composition의 추가 설명력은 통계적으로 강하지 않았다. 이 경우 결론은 위치보다 총 결합량/binding burden 쪽에 더 가깝다."
        )

    if amount_sig:
        amount_text = "total CLIP amount는 covariate-only model보다 유의하게 DeltaRD를 더 설명했다."
    else:
        amount_text = "total CLIP amount의 증분 설명력도 강하지 않아, baseline/필터/point-count 정의를 함께 점검해야 한다."

    dose_rho = dose_summary["overall_spearman_rho"].iloc[0]
    dose_p = dose_summary["overall_spearman_p"].iloc[0]
    m2_coeff = coeffs.loc[coeffs["model"] == "M2"].set_index("term")[
        "standardized_beta"
    ]

    lines = [
        "# Pipeline result summary",
        "",
        "## 실행 산출물",
        "",
        "- `results/tables/`: gene-level response, region counts, model summaries",
        "- `results/work/`: BED/intermediate count files",
        "- `results/figures/`: final report figures",
        "",
        "## QC",
        "",
        f"- Actb strandedness check: same-strand fraction = {stranded.get('same_strand_fraction', pd.Series([np.nan])).iloc[0]:.4f}. 이 값이 높으면 sense-stranded counting(`-s 1`, `bedtools -s`)이 타당하다.",
        f"- CLIP alignments seen = {int(clip_stats['sam_alignments_seen'].iloc[0]):,}; NH=1 5-prime points written = {int(clip_stats['nh1_points_written'].iloc[0]):,}.",
        f"- Composition-eligible model set n = {int(m2['n']):,}. M0/M1/M2는 동일 n에서 적합했다.",
        "",
        "## Baseline 재현",
        "",
        f"- baseline featureCounts CLIP enrichment vs DeltaRD: n = {baseline_stats['baseline_n']:,}, Pearson r = {baseline_stats['baseline_pearson_r']:.4f}, p = {fmt_p(baseline_stats['baseline_pearson_p'])}; Spearman rho = {baseline_stats['baseline_spearman_rho']:.4f}.",
        "- 이 그림은 원 논문/W1의 transcript-level total binding signal이 knockdown 후 ribosome density change와 연결되는지 확인하는 QC 역할이다.",
        "",
        "## Nested model 결과",
        "",
        f"- M0 R2 = {m0['r2']:.4f}",
        f"- M1 R2 = {m1['r2']:.4f}; DeltaR2_amount = {m1['delta_r2_vs_previous']:.4f}; F-test p = {fmt_p(m1['nested_p_vs_previous'])}",
        f"- M2 R2 = {m2['r2']:.4f}; DeltaR2_position = {m2['delta_r2_vs_previous']:.4f}; F-test p = {fmt_p(m2['nested_p_vs_previous'])}",
        f"- M2 standardized beta: total_CLIP_amount_point = {m2_coeff.get('total_CLIP_amount_point', np.nan):.4f}, f_cds = {m2_coeff.get('f_cds', np.nan):.4f}, f_3utr = {m2_coeff.get('f_3utr', np.nan):.4f}.",
        "- 방향 해석: M2에서 `f_cds` 계수는 음수이므로, total amount를 통제한 뒤 CDS 쪽 adjusted composition이 높을수록 derepression이 강해진다는 방향은 아니다. 위치 효과는 작고 보조적이며, 특히 CDS-positive story로 과장하면 안 된다.",
        f"- 해석: {amount_text} {position_text}",
        "",
        "## Dose 분석",
        "",
        f"- total exon-assigned CLIP 5-prime point count vs DeltaRD Spearman rho = {dose_rho:.4f}, p = {fmt_p(dose_p)}.",
        "- dose quartile별 ribosome density 변화 분포는 `fig4_dose_response.png`와 `tables/dose_summary.tsv`에서 확인한다.",
        "",
        "## 변수별 상관 요약",
        "",
        markdown_table(assoc),
        "",
        "## 표준화 회귀계수",
        "",
        markdown_table(coeffs.loc[coeffs["model"] == "M2"]),
        "",
        "## Eligibility flow",
        "",
        markdown_table(eligibility),
        "",
        "## 결론 가이드",
        "",
        "- 유의미한 분석인가: 같은 gene set에서 nested comparison을 수행했고, amount와 RNA-control adjusted position composition을 같은 stranded point signal에서 유도했으므로 설계상 주제 질문에 직접 답한다.",
        "- 단, 조건당 생물학적 반복이 없으므로 gene별 differential translation 유의성은 주장하지 않는다. 결과는 across-gene association이다.",
        "- `DeltaR2_position`이 작거나 비유의적이면 negative result도 의미가 있다. 이는 LIN28A의 번역 억제 효과가 region-specific preference보다 전체 결합량에 더 잘 설명될 수 있음을 뜻한다.",
    ]
    (paths.results / "result_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    paths = init_paths(args)
    cleanup_stale_figures(paths)
    write_run_manifest(paths, args)

    for required in [paths.read_counts, paths.gtf, paths.clip_bam, paths.rna_control_bam]:
        if not required.exists():
            raise FileNotFoundError(required)

    gene_response = compute_gene_response(paths, args)
    build_representative_regions(paths)
    check_strandedness(paths)
    clip_stats = extract_clip_5end_points(paths)
    bedtools_count_regions(paths)
    region_counts = aggregate_region_counts(paths)
    table = make_analysis_table(paths, gene_response, region_counts, args, clip_stats)

    model_summary, model_df = fit_models(paths, table)
    # fig1 workflow schematic removed: the analysis flow is described as text in the report.
    baseline_stats = plot_baseline(paths, gene_response)
    plot_composition(paths, table)
    plot_models(paths, model_summary, model_df)
    dose_summary = plot_dose(paths, table)
    write_summary(paths, baseline_stats, model_summary, dose_summary)

    log(f"done. Results: {paths.results}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[pipeline:error] {exc}", file=sys.stderr)
        raise
