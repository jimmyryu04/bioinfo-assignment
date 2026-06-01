#!/usr/bin/env python3
"""Week 6: select RNA-independent translational derepression candidates."""

from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ryu/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ASSIGNMENTS = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

GROUP_METRICS_PATH = ASSIGNMENTS / "week5" / "results" / "clip_group_gene_metrics.tsv"

CANDIDATE_RNA_ABS_MAX = 0.5
CANDIDATE_TE_MIN = 0.75
KNOWN_VALIDATION_GENES = ["Lamp1", "Epcam", "Cdh1"]


def load_group_metrics() -> pd.DataFrame:
    df = pd.read_csv(GROUP_METRICS_PATH, sep="\t")
    required = ["Geneid", "gene_name", "CLIP_enrichment", "RNA_log2FC", "TE_log2FC", "is_top20_clip"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns from week5 grouped metrics: {missing}")
    return df


def select_candidates(df: pd.DataFrame) -> pd.DataFrame:
    mask = (
        df["is_top20_clip"].astype(bool)
        & (df["RNA_log2FC"].abs() <= CANDIDATE_RNA_ABS_MAX)
        & (df["TE_log2FC"] >= CANDIDATE_TE_MIN)
    )
    candidates = df.loc[mask].copy()
    candidates = candidates.sort_values(
        ["TE_log2FC", "CLIP_enrichment"], ascending=[False, False]
    )
    candidates.insert(0, "candidate_rank", np.arange(1, len(candidates) + 1))
    candidates["candidate_rule"] = (
        f"CLIP top20%; abs(RNA_log2FC)<={CANDIDATE_RNA_ABS_MAX}; "
        f"TE_log2FC>={CANDIDATE_TE_MIN}"
    )
    candidates["symbol_keyword_category"] = candidates["gene_name"].fillna("").map(
        categorize_symbol
    )
    return candidates


def categorize_symbol(symbol: str) -> str:
    s = symbol.lower()
    if not s:
        return "other"
    if re.match(r"^tmem", s):
        return "transmembrane_symbol"
    if re.match(r"^(slc|kcn|clcn|cacn|atp)", s):
        return "transporter_or_channel_symbol"
    if re.match(r"^(lamp|laptm|sort|stx|vamp|rab|sec|cop|golg|golim|eea)", s):
        return "vesicle_secretory_symbol"
    if re.match(r"^(tmx|pdia|erp|hspa5|calr|canx|dnajc)", s):
        return "er_redox_chaperone_symbol"
    if re.match(r"^(elovl|acsl|scd|soat|sptlc|sgms|pla2|agpat)", s):
        return "lipid_membrane_metabolism_symbol"
    if re.match(r"^(cdh|itga|itgb|epcam|cldn|jam|icam|lamb|lamc|lama|lrp|notch|egfr)", s):
        return "cell_surface_adhesion_symbol"
    return "other"


def write_selection_steps(df: pd.DataFrame, candidates: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    top20 = df["is_top20_clip"].astype(bool)
    stable_rna = df["RNA_log2FC"].abs() <= CANDIDATE_RNA_ABS_MAX
    high_te = df["TE_log2FC"] >= CANDIDATE_TE_MIN
    rows = [
        {"selection_step": "base_filtered_genes_from_week5", "n_genes": int(len(df))},
        {"selection_step": "CLIP_top20", "n_genes": int(top20.sum())},
        {
            "selection_step": "CLIP_top20_and_abs_RNA_log2FC_le_0.5",
            "n_genes": int((top20 & stable_rna).sum()),
        },
        {
            "selection_step": "CLIP_top20_and_TE_log2FC_ge_0.75",
            "n_genes": int((top20 & high_te).sum()),
        },
        {"selection_step": "final_candidates", "n_genes": int(len(candidates))},
    ]
    table = pd.DataFrame(rows)
    table.to_csv(out_path, sep="\t", index=False)
    return table


def write_candidate_summary(candidates: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    metrics = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    rows = []
    for metric in metrics:
        values = candidates[metric].dropna()
        rows.append(
            {
                "metric": metric,
                "n_candidates": int(values.size),
                "mean": values.mean(),
                "median": values.median(),
                "min": values.min(),
                "q25": values.quantile(0.25),
                "q75": values.quantile(0.75),
                "max": values.max(),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return summary


def write_keyword_summary(candidates: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    counts = (
        candidates["symbol_keyword_category"]
        .value_counts()
        .rename_axis("symbol_keyword_category")
        .reset_index(name="n_candidates")
    )
    counts["fraction"] = counts["n_candidates"] / max(len(candidates), 1)
    counts.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return counts


def write_known_gene_table(df: pd.DataFrame, candidates: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    candidate_ids = set(candidates["Geneid"])
    rows = []
    for gene in KNOWN_VALIDATION_GENES:
        matches = df.loc[df["gene_name"] == gene].copy()
        if matches.empty:
            rows.append({"gene_name": gene, "present_in_filtered_table": False})
            continue
        for _, row in matches.iterrows():
            rows.append(
                {
                    "Geneid": row["Geneid"],
                    "gene_name": gene,
                    "present_in_filtered_table": True,
                    "is_candidate": row["Geneid"] in candidate_ids,
                    "CLIP_enrichment": row["CLIP_enrichment"],
                    "RNA_log2FC": row["RNA_log2FC"],
                    "TE_log2FC": row["TE_log2FC"],
                    "CLIP_group": row["CLIP_group"],
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return table


def plot_candidate_scatter(df: pd.DataFrame, candidates: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.scatter(
        df["CLIP_enrichment"],
        df["TE_log2FC"],
        s=7,
        color="#B8B8B8",
        alpha=0.35,
        edgecolors="none",
        label="filtered genes",
    )
    ax.scatter(
        candidates["CLIP_enrichment"],
        candidates["TE_log2FC"],
        s=18,
        color="#D95F02",
        alpha=0.9,
        edgecolors="none",
        label="candidate targets",
    )
    ax.axhline(CANDIDATE_TE_MIN, color="#D95F02", lw=1.1, ls="--")
    ax.axvline(df["CLIP_enrichment"].quantile(0.80), color="#4C78A8", lw=1.1, ls="--")
    for _, row in candidates.head(10).iterrows():
        ax.annotate(
            row["gene_name"],
            (row["CLIP_enrichment"], row["TE_log2FC"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("LIN28A CLIP enrichment (log2)")
    ax.set_ylabel("TE_log2FC")
    ax.set_title("Week 6: RNA-independent TE derepression candidates", loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_top_candidate_bars(candidates: pd.DataFrame, out_path: Path) -> None:
    top = candidates.head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    colors = plt.cm.viridis(
        (top["CLIP_enrichment"] - top["CLIP_enrichment"].min())
        / (top["CLIP_enrichment"].max() - top["CLIP_enrichment"].min() + 1e-9)
    )
    ax.barh(top["gene_name"], top["TE_log2FC"], color=colors, edgecolor="none")
    ax.axvline(CANDIDATE_TE_MIN, color="#D95F02", lw=1.1, ls="--")
    ax.set_xlabel("TE_log2FC")
    ax.set_ylabel("candidate gene")
    ax.set_title("Top candidate targets ranked by TE derepression", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    sm = plt.cm.ScalarMappable(
        cmap="viridis",
        norm=plt.Normalize(top["CLIP_enrichment"].min(), top["CLIP_enrichment"].max()),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("CLIP enrichment")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_metric_heatmap(candidates: pd.DataFrame, out_path: Path) -> None:
    top = candidates.head(30).copy()
    metrics = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    matrix = top[metrics].to_numpy(dtype=float)
    col_mean = matrix.mean(axis=0)
    col_std = matrix.std(axis=0)
    z = (matrix - col_mean) / np.where(col_std == 0, 1, col_std)

    fig, ax = plt.subplots(figsize=(6.5, 8.0))
    im = ax.imshow(z, cmap="RdBu_r", aspect="auto", vmin=-2.2, vmax=2.2)
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(metrics, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["gene_name"])
    ax.set_title("Top candidates: relative metric pattern", loc="left")
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=6)
    cbar = fig.colorbar(im, ax=ax, shrink=0.75)
    cbar.set_label("within-column z-score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def write_report(
    df: pd.DataFrame,
    candidates: pd.DataFrame,
    selection_steps: pd.DataFrame,
    keyword_summary: pd.DataFrame,
    known_gene_table: pd.DataFrame,
) -> None:
    n_candidates = len(candidates)
    top_names = ", ".join(candidates["gene_name"].head(10).tolist())
    non_other_summary = keyword_summary.loc[keyword_summary["symbol_keyword_category"] != "other"]
    non_other = int(non_other_summary["n_candidates"].sum())
    top_category = non_other_summary.iloc[0] if not non_other_summary.empty else keyword_summary.iloc[0]
    known_present = known_gene_table.loc[
        known_gene_table.get("present_in_filtered_table", False).astype(bool)
        if "present_in_filtered_table" in known_gene_table.columns
        else []
    ]
    known_lines = []
    for _, row in known_present.iterrows():
        known_lines.append(
            f"- {row['gene_name']}: candidate={bool(row['is_candidate'])}, "
            f"CLIP={row['CLIP_enrichment']:.3f}, RNA={row['RNA_log2FC']:.3f}, TE={row['TE_log2FC']:.3f}"
        )
    known_text = "\n".join(known_lines) if known_lines else "- No validation genes were present in the filtered table."

    report = f"""# Week 6 Report: RNA-independent candidate targets

## Goal

Select genes that are strongly bound by LIN28A and show translation efficiency
increase after Lin28a knockdown, while RNA abundance remains comparatively
stable.

## Candidate rule

- CLIP enrichment in the top 20% of week5 filtered genes
- `abs(RNA_log2FC) <= {CANDIDATE_RNA_ABS_MAX}`
- `TE_log2FC >= {CANDIDATE_TE_MIN}`

This produced {n_candidates:,} candidate targets from {len(df):,} filtered genes.
The top candidates by TE_log2FC were: {top_names}.

## Biological interpretation

The candidate list is not a validated direct-target list. It is a computational
set consistent with RNA-independent translational derepression. Still,
{non_other:,} of {n_candidates:,} candidates matched simple gene-symbol keyword
categories related to membrane, vesicle/secretory, ER redox/chaperone,
transport, lipid membrane metabolism, or cell-surface genes. Among these
non-`other` labels, the largest category was `{top_category['symbol_keyword_category']}` with
{int(top_category['n_candidates'])} genes.

## Paper validation genes checked

The paper validated LAMP1, EpCAM, and E-cadherin/Cdh1 by western blot. Their
status in this filtered gene-level analysis was:

{known_text}

## Outputs

- `results/rna_independent_te_targets.tsv`
- `results/top_candidate_targets.tsv`
- `results/candidate_selection_steps.tsv`
- `results/candidate_metric_summary.tsv`
- `results/candidate_keyword_summary.tsv`
- `results/paper_validation_gene_check.tsv`
- `figures/fig1_candidate_scatter.png`
- `figures/fig2_top_candidate_barplot.png`
- `figures/fig3_candidate_metric_heatmap.png`
"""
    (HERE / "report.md").write_text(report)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    df = load_group_metrics()
    candidates = select_candidates(df)

    selection_steps = write_selection_steps(df, candidates, RESULTS / "candidate_selection_steps.tsv")
    write_candidate_summary(candidates, RESULTS / "candidate_metric_summary.tsv")
    keyword_summary = write_keyword_summary(candidates, RESULTS / "candidate_keyword_summary.tsv")
    known_gene_table = write_known_gene_table(df, candidates, RESULTS / "paper_validation_gene_check.tsv")

    candidates.to_csv(
        RESULTS / "rna_independent_te_targets.tsv",
        sep="\t",
        index=False,
        float_format="%.8g",
    )
    candidates.head(30).to_csv(
        RESULTS / "top_candidate_targets.tsv",
        sep="\t",
        index=False,
        float_format="%.8g",
    )

    plot_candidate_scatter(df, candidates, FIGURES / "fig1_candidate_scatter.png")
    plot_top_candidate_bars(candidates, FIGURES / "fig2_top_candidate_barplot.png")
    plot_metric_heatmap(candidates, FIGURES / "fig3_candidate_metric_heatmap.png")
    write_report(df, candidates, selection_steps, keyword_summary, known_gene_table)

    print(f"Week 6 complete: selected {len(candidates):,} candidate targets")
    print(f"Top candidate: {candidates.iloc[0]['gene_name'] if len(candidates) else 'none'}")
    print(f"Output directory: {HERE}")


if __name__ == "__main__":
    main()
