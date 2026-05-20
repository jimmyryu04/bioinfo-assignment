# Week 4 Report: metric table preparation

## Goal

Prepare a reusable gene-level table for the own analysis project. The table
separates RNA abundance change from ribosome footprint change so that later
weeks can test translation efficiency derepression.

## Inputs

- Count table: `../week1/work/read-counts.txt`
- Annotation: `../week1/work/binfo1-datapack1/gencode.gtf`
- Pseudocount: `0.1` CPM

## Metrics

- `CLIP_enrichment = log2((CPM_CLIP + pc) / (CPM_RNA_control + pc))`
- `RNA_log2FC = log2((CPM_RNA_siLin28a + pc) / (CPM_RNA_siLuc + pc))`
- `RPF_log2FC = log2((CPM_RPF_siLin28a + pc) / (CPM_RPF_siLuc + pc))`
- `TE_log2FC = RPF_log2FC - RNA_log2FC`

## Base filter

The main analysis keeps protein-coding, non-histone genes with RNA counts at
least 30 in `RNA-control`, `RNA-siLuc`, and `RNA-siLin28a`, plus at least 80
RPF reads in `RPF-siLuc`.

Starting from 55,359 featureCounts rows, 7,840 genes passed the
base filter. The median CLIP enrichment among filtered genes was -0.595,
and the median TE change was -0.278.

## Outputs

- `results/gene_metrics.tsv`: all genes with raw counts, CPM values, metrics,
  and filter flags.
- `results/filter_summary.tsv`: number of genes retained at each filter step.
- `results/metric_summary.tsv`: summary statistics for the filtered genes.
- `figures/week4_metric_distributions.png`: distributions of the four core
  metrics.
