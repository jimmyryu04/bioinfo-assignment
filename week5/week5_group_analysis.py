#!/usr/bin/env python3
"""Week 5: correlation, CLIP-group, and robustness analyses."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ryu/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu, pearsonr, spearmanr


ASSIGNMENTS = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

METRICS_PATH = ASSIGNMENTS / "week4" / "results" / "gene_metrics.tsv"
GROUP_ORDER = ["top 5%", "5-20%", "20-50%", "50-100%"]


def load_metrics() -> pd.DataFrame:
    df = pd.read_csv(METRICS_PATH, sep="\t")
    required = [
        "CLIP_enrichment",
        "RNA_log2FC",
        "RPF_log2FC",
        "TE_log2FC",
        "pass_base_filter",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns from week4 metrics: {missing}")
    return df


def base_filtered(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    out = df.loc[df["pass_base_filter"].astype(bool)].copy()
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=metric_cols)
    return out


def assign_clip_groups(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    q95 = out["CLIP_enrichment"].quantile(0.95)
    q80 = out["CLIP_enrichment"].quantile(0.80)
    q50 = out["CLIP_enrichment"].quantile(0.50)
    conditions = [
        out["CLIP_enrichment"] >= q95,
        (out["CLIP_enrichment"] >= q80) & (out["CLIP_enrichment"] < q95),
        (out["CLIP_enrichment"] >= q50) & (out["CLIP_enrichment"] < q80),
        out["CLIP_enrichment"] < q50,
    ]
    out["CLIP_group"] = np.select(conditions, GROUP_ORDER, default="unassigned")
    out["CLIP_group"] = pd.Categorical(out["CLIP_group"], GROUP_ORDER, ordered=True)
    out["is_top20_clip"] = out["CLIP_enrichment"] >= q80
    return out


def write_correlation_summary(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    rows = []
    for metric in ["RNA_log2FC", "RPF_log2FC", "TE_log2FC"]:
        x = df["CLIP_enrichment"]
        y = df[metric]
        pr = pearsonr(x, y)
        sr = spearmanr(x, y)
        rows.append(
            {
                "x": "CLIP_enrichment",
                "y": metric,
                "n": int(len(df)),
                "pearson_r": pr.statistic,
                "pearson_p": pr.pvalue,
                "spearman_r": sr.statistic,
                "spearman_p": sr.pvalue,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return summary


def write_group_summary(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    metrics = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    rows = []
    for group in GROUP_ORDER:
        part = df.loc[df["CLIP_group"] == group]
        row = {"CLIP_group": group, "n": int(len(part))}
        for metric in metrics:
            row[f"{metric}_mean"] = part[metric].mean()
            row[f"{metric}_median"] = part[metric].median()
            row[f"{metric}_q25"] = part[metric].quantile(0.25)
            row[f"{metric}_q75"] = part[metric].quantile(0.75)
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return summary


def write_group_tests(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    rows = []
    comparisons = {
        "top_5_vs_bottom_50": (df["CLIP_group"] == "top 5%", df["CLIP_group"] == "50-100%"),
        "top_20_vs_bottom_50": (df["is_top20_clip"], df["CLIP_group"] == "50-100%"),
    }
    for metric in ["RNA_log2FC", "RPF_log2FC", "TE_log2FC"]:
        for label, (test_mask, ref_mask) in comparisons.items():
            test_values = df.loc[test_mask, metric].dropna()
            ref_values = df.loc[ref_mask, metric].dropna()
            mw = mannwhitneyu(test_values, ref_values, alternative="two-sided")
            ks = ks_2samp(test_values, ref_values, alternative="two-sided", mode="auto")
            rows.append(
                {
                    "metric": metric,
                    "comparison": label,
                    "n_test": int(test_values.size),
                    "n_reference": int(ref_values.size),
                    "median_test": test_values.median(),
                    "median_reference": ref_values.median(),
                    "median_difference": test_values.median() - ref_values.median(),
                    "mannwhitney_u": mw.statistic,
                    "mannwhitney_p": mw.pvalue,
                    "ks_statistic": ks.statistic,
                    "ks_p": ks.pvalue,
                }
            )
    tests = pd.DataFrame(rows)
    tests.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return tests


def threshold_filter(df: pd.DataFrame, rna_min: int, rpf_min: int) -> pd.Series:
    metric_cols = ["CLIP_enrichment", "RNA_log2FC", "RPF_log2FC", "TE_log2FC"]
    finite = np.isfinite(df[metric_cols]).all(axis=1)
    gene_name = df["gene_name"].fillna("")
    return (
        df["gene_type"].eq("protein_coding")
        & ~gene_name.str.startswith("Hist")
        & (df["RNA-control.bam"] >= rna_min)
        & (df["RNA-siLuc.bam"] >= rna_min)
        & (df["RNA-siLin28a.bam"] >= rna_min)
        & (df["RPF-siLuc.bam"] >= rpf_min)
        & finite
    )


def write_threshold_sensitivity(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    rows = []
    for rna_min in [10, 30, 50]:
        for rpf_min in [50, 80, 100]:
            mask = threshold_filter(df, rna_min, rpf_min)
            part = df.loc[mask].copy()
            if len(part) >= 3:
                clip = part["CLIP_enrichment"]
                te = part["TE_log2FC"]
                rna = part["RNA_log2FC"]
                pearson_te = pearsonr(clip, te).statistic
                spearman_te = spearmanr(clip, te).statistic
                pearson_rna = pearsonr(clip, rna).statistic
                spearman_rna = spearmanr(clip, rna).statistic
                q80 = clip.quantile(0.80)
                candidate_count = int(
                    (
                        (clip >= q80)
                        & (part["RNA_log2FC"].abs() <= 0.5)
                        & (part["TE_log2FC"] >= 0.75)
                    ).sum()
                )
            else:
                pearson_te = spearman_te = pearson_rna = spearman_rna = np.nan
                candidate_count = 0
            rows.append(
                {
                    "rna_min_count": rna_min,
                    "rpf_siLuc_min_count": rpf_min,
                    "n_genes": int(len(part)),
                    "pearson_clip_vs_RNA": pearson_rna,
                    "spearman_clip_vs_RNA": spearman_rna,
                    "pearson_clip_vs_TE": pearson_te,
                    "spearman_clip_vs_TE": spearman_te,
                    "candidate_count_top20_absRNA0.5_TE0.75": candidate_count,
                }
            )
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    return sensitivity


def plot_scatter(df: pd.DataFrame, correlations: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), sharex=True)
    panels = [
        ("RNA_log2FC", "RNA abundance change"),
        ("TE_log2FC", "Translation efficiency change"),
    ]
    for ax, (metric, title) in zip(axes, panels):
        ax.scatter(
            df["CLIP_enrichment"],
            df[metric],
            s=8,
            color="#333333",
            alpha=0.22,
            edgecolors="none",
        )
        x = df["CLIP_enrichment"]
        y = df[metric]
        slope, intercept = np.polyfit(x, y, 1)
        xx = np.linspace(x.quantile(0.005), x.quantile(0.995), 100)
        ax.plot(xx, slope * xx + intercept, color="#D95F02", lw=1.8)
        ax.axhline(0, color="black", lw=0.8, alpha=0.45)
        ax.axvline(0, color="black", lw=0.8, alpha=0.45)
        row = correlations.loc[correlations["y"] == metric].iloc[0]
        ax.text(
            0.04,
            0.96,
            f"Pearson r = {row['pearson_r']:.3f}\nSpearman r = {row['spearman_r']:.3f}",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("LIN28A CLIP enrichment (log2)")
        ax.set_ylabel(metric)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        "Week 5: CLIP enrichment is more strongly linked to TE than RNA abundance",
        x=0.01,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_group_boxplots(df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharex=True)
    panels = [("RNA_log2FC", "RNA abundance"), ("TE_log2FC", "Translation efficiency")]
    colors = ["#4C78A8", "#72B7B2", "#F58518", "#B279A2"]
    positions = np.arange(1, len(GROUP_ORDER) + 1)
    for ax, (metric, title) in zip(axes, panels):
        data = [df.loc[df["CLIP_group"] == group, metric].dropna().values for group in GROUP_ORDER]
        bp = ax.boxplot(
            data,
            positions=positions,
            patch_artist=True,
            widths=0.62,
            showfliers=False,
            medianprops={"color": "black", "lw": 1.2},
        )
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
            patch.set_edgecolor("#333333")
        ax.axhline(0, color="black", lw=0.8, alpha=0.45)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(metric)
        ax.set_xticks(positions)
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        "Week 5: CLIP-high genes show stronger TE derepression",
        x=0.01,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_threshold_sensitivity(sensitivity: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    heatmaps = [
        ("spearman_clip_vs_TE", "Spearman r: CLIP vs TE"),
        ("candidate_count_top20_absRNA0.5_TE0.75", "Candidate count"),
    ]
    for ax, (value_col, title) in zip(axes, heatmaps):
        pivot = sensitivity.pivot(
            index="rna_min_count", columns="rpf_siLuc_min_count", values=value_col
        ).sort_index(ascending=False)
        im = ax.imshow(pivot.values, cmap="viridis", aspect="auto")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("RPF-siLuc min count")
        ax.set_ylabel("RNA min count")
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([str(c) for c in pivot.columns])
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels([str(i) for i in pivot.index])
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                value = pivot.values[i, j]
                label = f"{value:.2f}" if value_col.startswith("spearman") else f"{int(value)}"
                ax.text(j, i, label, ha="center", va="center", color="white", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.78)
    fig.suptitle(
        "Week 5: robustness across count thresholds",
        x=0.01,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def write_report(
    correlations: pd.DataFrame,
    group_summary: pd.DataFrame,
    group_tests: pd.DataFrame,
    sensitivity: pd.DataFrame,
) -> None:
    rna_row = correlations.loc[correlations["y"] == "RNA_log2FC"].iloc[0]
    te_row = correlations.loc[correlations["y"] == "TE_log2FC"].iloc[0]
    top5 = group_summary.loc[group_summary["CLIP_group"] == "top 5%"].iloc[0]
    bottom = group_summary.loc[group_summary["CLIP_group"] == "50-100%"].iloc[0]
    te_test = group_tests.loc[
        (group_tests["metric"] == "TE_log2FC")
        & (group_tests["comparison"] == "top_20_vs_bottom_50")
    ].iloc[0]
    te_p = te_test["mannwhitney_p"]
    te_p_text = "< 1e-300" if te_p == 0 else f"{te_p:.3e}"
    min_te_corr = sensitivity["spearman_clip_vs_TE"].min()
    max_te_corr = sensitivity["spearman_clip_vs_TE"].max()

    report = f"""# Week 5 Report: CLIP enrichment versus RNA and TE changes

## Goal

Use the week4 metric table to test whether LIN28A binding is linked more
strongly to translation efficiency change than to RNA abundance change.

## Main result

Among base-filtered genes, CLIP enrichment showed only a weak negative
relationship with RNA abundance change:

- Pearson r = {rna_row['pearson_r']:.3f}
- Spearman r = {rna_row['spearman_r']:.3f}

By contrast, CLIP enrichment was positively associated with TE change:

- Pearson r = {te_row['pearson_r']:.3f}
- Spearman r = {te_row['spearman_r']:.3f}

This matches the paper's central model: LIN28A-bound transcripts are affected
mainly at the translation level after Lin28a knockdown.

## CLIP group comparison

The top 5% CLIP-enriched genes had median TE_log2FC = {top5['TE_log2FC_median']:.3f},
whereas the bottom 50% had median TE_log2FC = {bottom['TE_log2FC_median']:.3f}.
For top 20% versus bottom 50%, the Mann-Whitney p-value for TE_log2FC was
{te_p_text}.

## Robustness

Across RNA count thresholds 10, 30, 50 and RPF-siLuc thresholds 50, 80, 100,
the Spearman correlation between CLIP enrichment and TE_log2FC ranged from
{min_te_corr:.3f} to {max_te_corr:.3f}. The conclusion is therefore not driven
by a single filtering threshold.

## Outputs

- `results/correlation_summary.tsv`
- `results/clip_group_gene_metrics.tsv`
- `results/clip_group_summary.tsv`
- `results/clip_group_tests.tsv`
- `results/threshold_sensitivity.tsv`
- `figures/fig1_clip_vs_rna_te.png`
- `figures/fig2_clip_group_boxplot.png`
- `figures/fig3_threshold_sensitivity.png`
"""
    (HERE / "report.md").write_text(report)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    all_metrics = load_metrics()
    filtered = assign_clip_groups(base_filtered(all_metrics))

    correlations = write_correlation_summary(filtered, RESULTS / "correlation_summary.tsv")
    group_summary = write_group_summary(filtered, RESULTS / "clip_group_summary.tsv")
    group_tests = write_group_tests(filtered, RESULTS / "clip_group_tests.tsv")
    sensitivity = write_threshold_sensitivity(all_metrics, RESULTS / "threshold_sensitivity.tsv")

    filtered.to_csv(
        RESULTS / "clip_group_gene_metrics.tsv",
        sep="\t",
        index=False,
        float_format="%.8g",
    )
    plot_scatter(filtered, correlations, FIGURES / "fig1_clip_vs_rna_te.png")
    plot_group_boxplots(filtered, FIGURES / "fig2_clip_group_boxplot.png")
    plot_threshold_sensitivity(sensitivity, FIGURES / "fig3_threshold_sensitivity.png")
    write_report(correlations, group_summary, group_tests, sensitivity)

    te_r = correlations.loc[correlations["y"] == "TE_log2FC", "spearman_r"].iloc[0]
    print(f"Week 5 complete: analyzed {len(filtered):,} filtered genes")
    print(f"Spearman(CLIP, TE_log2FC) = {te_r:.3f}")
    print(f"Output directory: {HERE}")


if __name__ == "__main__":
    main()
