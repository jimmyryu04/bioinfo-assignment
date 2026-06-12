#!/usr/bin/env python3
"""Regenerate report figures from already-computed tables, without re-running
the heavy BAM/bedtools steps. Picks up the current style block in run_pipeline
(ggplot + Inter). fig1 (workflow) is intentionally not produced."""
from pathlib import Path

import pandas as pd

import run_pipeline as rp

PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"

paths = rp.Paths(
    project=PROJECT,
    pipeline=PROJECT / "pipeline",
    data_dir=PROJECT / "week1/work/binfo1-datapack1",
    read_counts=PROJECT / "week1/work/read-counts.txt",
    gtf=PROJECT / "week1/work/binfo1-datapack1/gencode.gtf",
    clip_bam=PROJECT / "week1/work/binfo1-datapack1/CLIP-35L33G.bam",
    rna_control_bam=PROJECT / "week1/work/binfo1-datapack1/RNA-control.bam",
    results=RESULTS,
    work=RESULTS / "work",
    tables=RESULTS / "tables",
    figures=RESULTS / "figures",
    logs=RESULTS / "logs",
    samtools="/blaze/ryu/conda/envs/lab/bin/samtools",
    bedtools="/blaze/ryu/conda/envs/lab/bin/bedtools",
)

gene_response = pd.read_csv(paths.tables / "gene_response.tsv", sep="\t")
table = pd.read_csv(paths.work / "analysis_table.tsv", sep="\t")
model_summary = pd.read_csv(paths.tables / "model_summary.tsv", sep="\t")
model_df = pd.read_csv(paths.tables / "model_gene_set.tsv", sep="\t")

rp.cleanup_stale_figures(paths)

print("font.sans-serif ->", rp.mpl.rcParams["font.sans-serif"][:2])
print("style facecolor ->", rp.mpl.rcParams["axes.facecolor"])

rp.plot_baseline(paths, gene_response)
rp.plot_composition(paths, table)
rp.plot_models(paths, model_summary, model_df)
rp.plot_dose(paths, table)
print("figures regenerated:", sorted(p.name for p in paths.figures.glob("*.png")))
