#!/usr/bin/env python3
"""Week 4: prepare gene-level LIN28A CLIP/RNA/RPF metrics.

This script starts from the featureCounts output generated in week1 and
creates a normalized gene-level metric table for the own analysis project.
The main question is whether LIN28A-bound genes show translation efficiency
derepression after Lin28a knockdown without a matching mRNA abundance change.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ASSIGNMENTS = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

COUNTS_PATH = ASSIGNMENTS / "week1" / "work" / "read-counts.txt"
GTF_PATH = ASSIGNMENTS / "week1" / "work" / "binfo1-datapack1" / "gencode.gtf"

SAMPLE_MAP = {
    "CLIP-35L33G.bam": "CLIP",
    "RNA-control.bam": "RNA_control",
    "RNA-siLin28a.bam": "RNA_siLin28a",
    "RNA-siLuc.bam": "RNA_siLuc",
    "RPF-siLin28a.bam": "RPF_siLin28a",
    "RPF-siLuc.bam": "RPF_siLuc",
}

RAW_SAMPLE_COLS = list(SAMPLE_MAP)
PSEUDOCOUNT_CPM = 0.1


def parse_gtf_attributes(attr_text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)";', attr_text))


def load_gene_annotation(gtf_path: Path) -> pd.DataFrame:
    rows: list[tuple[str | None, str | None, str | None]] = []
    with gtf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            rows.append(
                (
                    attrs.get("gene_id"),
                    attrs.get("gene_name"),
                    attrs.get("gene_type") or attrs.get("gene_biotype"),
                )
            )
    ann = pd.DataFrame(rows, columns=["Geneid", "gene_name", "gene_type"])
    return ann.drop_duplicates("Geneid")


def load_counts(counts_path: Path) -> pd.DataFrame:
    counts = pd.read_csv(counts_path, sep="\t", comment="#")
    missing = [col for col in RAW_SAMPLE_COLS if col not in counts.columns]
    if missing:
        raise ValueError(f"Missing expected count columns: {missing}")
    return counts


def add_cpm_and_metrics(counts: pd.DataFrame) -> pd.DataFrame:
    out = counts.copy()
    library_sizes = out[RAW_SAMPLE_COLS].sum(axis=0)

    for raw_col, clean_name in SAMPLE_MAP.items():
        out[f"CPM_{clean_name}"] = out[raw_col] / library_sizes[raw_col] * 1_000_000

    pc = PSEUDOCOUNT_CPM
    out["CLIP_enrichment"] = np.log2(
        (out["CPM_CLIP"] + pc) / (out["CPM_RNA_control"] + pc)
    )
    out["RNA_log2FC"] = np.log2(
        (out["CPM_RNA_siLin28a"] + pc) / (out["CPM_RNA_siLuc"] + pc)
    )
    out["RPF_log2FC"] = np.log2(
        (out["CPM_RPF_siLin28a"] + pc) / (out["CPM_RPF_siLuc"] + pc)
    )
    out["TE_log2FC"] = out["RPF_log2FC"] - out["RNA_log2FC"]

    out["library_size_normalized"] = True
    out["pseudocount_cpm"] = PSEUDOCOUNT_CPM
    return out


def add_filter_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    gene_name = out["gene_name"].fillna("")
    out["is_protein_coding"] = out["gene_type"].eq("protein_coding")
    out["is_non_histone"] = ~gene_name.str.startswith("Hist")
    out["pass_rna_count_filter"] = (
        (out["RNA-control.bam"] >= 30)
        & (out["RNA-siLuc.bam"] >= 30)
        & (out["RNA-siLin28a.bam"] >= 30)
    )
    out["pass_rpf_count_filter"] = out["RPF-siLuc.bam"] >= 80
    metric_cols = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    out["pass_finite_metric_filter"] = np.isfinite(out[metric_cols]).all(axis=1)
    out["pass_base_filter"] = (
        out["is_protein_coding"]
        & out["is_non_histone"]
        & out["pass_rna_count_filter"]
        & out["pass_rpf_count_filter"]
        & out["pass_finite_metric_filter"]
    )
    return out


def write_filter_summary(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    steps = [
        ("all_featureCounts_genes", pd.Series(True, index=df.index)),
        ("protein_coding", df["is_protein_coding"]),
        ("protein_coding_non_histone", df["is_protein_coding"] & df["is_non_histone"]),
        (
            "plus_RNA_count_filter",
            df["is_protein_coding"] & df["is_non_histone"] & df["pass_rna_count_filter"],
        ),
        ("plus_RPF_count_filter", df["pass_base_filter"]),
    ]
    rows = [{"filter_step": name, "n_genes": int(mask.sum())} for name, mask in steps]
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False)
    return summary


def write_metric_summary(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    metrics = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    filtered = df.loc[df["pass_base_filter"], metrics]
    rows = []
    for metric in metrics:
        values = filtered[metric].dropna()
        rows.append(
            {
                "metric": metric,
                "n": int(values.size),
                "mean": values.mean(),
                "median": values.median(),
                "std": values.std(),
                "min": values.min(),
                "q25": values.quantile(0.25),
                "q75": values.quantile(0.75),
                "max": values.max(),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return summary


def plot_metric_distributions(df: pd.DataFrame, out_path: Path) -> None:
    metrics = [
        ("CLIP_enrichment", "LIN28A CLIP enrichment"),
        ("RNA_log2FC", "RNA change after Lin28a KD"),
        ("RPF_log2FC", "RPF change after Lin28a KD"),
        ("TE_log2FC", "Translation efficiency change"),
    ]
    filtered = df.loc[df["pass_base_filter"]]

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.8))
    axes = axes.ravel()
    for ax, (metric, label) in zip(axes, metrics):
        values = filtered[metric].replace([np.inf, -np.inf], np.nan).dropna()
        lo, hi = values.quantile([0.01, 0.99])
        bins = np.linspace(lo, hi, 50)
        ax.hist(values.clip(lo, hi), bins=bins, color="#4C78A8", alpha=0.85)
        ax.axvline(values.median(), color="#D95F02", lw=1.5, label="median")
        ax.axvline(0, color="black", lw=0.8, alpha=0.5)
        ax.set_title(label, fontsize=10)
        ax.set_xlabel(metric)
        ax.set_ylabel("genes")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Week 4: normalized gene-level metrics after base filtering",
        x=0.01,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def write_report(filter_summary: pd.DataFrame, metric_summary: pd.DataFrame) -> None:
    n_total = int(filter_summary.loc[0, "n_genes"])
    n_filtered = int(filter_summary.loc[filter_summary["filter_step"] == "plus_RPF_count_filter", "n_genes"].iloc[0])
    te_median = metric_summary.loc[metric_summary["metric"] == "TE_log2FC", "median"].iloc[0]
    clip_median = metric_summary.loc[metric_summary["metric"] == "CLIP_enrichment", "median"].iloc[0]

    report = f"""# Week 4 Report: metric table preparation

## Goal

Prepare a reusable gene-level table for the own analysis project. The table
separates RNA abundance change from ribosome footprint change so that later
weeks can test translation efficiency derepression.

## Inputs

- Count table: `../week1/work/read-counts.txt`
- Annotation: `../week1/work/binfo1-datapack1/gencode.gtf`
- Pseudocount: `{PSEUDOCOUNT_CPM}` CPM

## Metrics

- `CLIP_enrichment = log2((CPM_CLIP + pc) / (CPM_RNA_control + pc))`
- `RNA_log2FC = log2((CPM_RNA_siLin28a + pc) / (CPM_RNA_siLuc + pc))`
- `RPF_log2FC = log2((CPM_RPF_siLin28a + pc) / (CPM_RPF_siLuc + pc))`
- `TE_log2FC = RPF_log2FC - RNA_log2FC`

## Base filter

The main analysis keeps protein-coding, non-histone genes with RNA counts at
least 30 in `RNA-control`, `RNA-siLuc`, and `RNA-siLin28a`, plus at least 80
RPF reads in `RPF-siLuc`.

Starting from {n_total:,} featureCounts rows, {n_filtered:,} genes passed the
base filter. The median CLIP enrichment among filtered genes was {clip_median:.3f},
and the median TE change was {te_median:.3f}.

## Outputs

- `results/gene_metrics.tsv`: all genes with raw counts, CPM values, metrics,
  and filter flags.
- `results/filter_summary.tsv`: number of genes retained at each filter step.
- `results/metric_summary.tsv`: summary statistics for the filtered genes.
- `figures/week4_metric_distributions.png`: distributions of the four core
  metrics.
"""
    (HERE / "report.md").write_text(report)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    counts = load_counts(COUNTS_PATH)
    annotation = load_gene_annotation(GTF_PATH)
    merged = counts.merge(annotation, on="Geneid", how="left")
    metrics = add_cpm_and_metrics(merged)
    metrics = add_filter_flags(metrics)

    metrics.to_csv(RESULTS / "gene_metrics.tsv", sep="\t", index=False, float_format="%.8g")
    filter_summary = write_filter_summary(metrics, RESULTS / "filter_summary.tsv")
    metric_summary = write_metric_summary(metrics, RESULTS / "metric_summary.tsv")
    plot_metric_distributions(metrics, FIGURES / "week4_metric_distributions.png")
    write_report(filter_summary, metric_summary)

    n_filtered = int(metrics["pass_base_filter"].sum())
    print(f"Week 4 complete: wrote metrics for {len(metrics):,} genes")
    print(f"Base-filtered genes: {n_filtered:,}")
    print(f"Output directory: {HERE}")


if __name__ == "__main__":
    main()
