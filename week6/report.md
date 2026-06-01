# Week 6 Report: RNA-independent candidate targets

## Goal

Select genes that are strongly bound by LIN28A and show translation efficiency
increase after Lin28a knockdown, while RNA abundance remains comparatively
stable.

## Candidate rule

- CLIP enrichment in the top 20% of week5 filtered genes
- `abs(RNA_log2FC) <= 0.5`
- `TE_log2FC >= 0.75`

This produced 177 candidate targets from 7,840 filtered genes.
The top candidates by TE_log2FC were: Dram2, Tmx1, Ltbp1, Lamp2, Tmem87a, Cspg4b, Tmem67, Armh4, Elovl7, Tm2d3.

## Biological interpretation

The candidate list is not a validated direct-target list. It is a computational
set consistent with RNA-independent translational derepression. Still,
35 of 177 candidates matched simple gene-symbol keyword
categories related to membrane, vesicle/secretory, ER redox/chaperone,
transport, lipid membrane metabolism, or cell-surface genes. Among these
non-`other` labels, the largest category was `transporter_or_channel_symbol` with
17 genes.

## Paper validation genes checked

The paper validated LAMP1, EpCAM, and E-cadherin/Cdh1 by western blot. Their
status in this filtered gene-level analysis was:

- Lamp1: candidate=True, CLIP=0.786, RNA=0.054, TE=0.797
- Epcam: candidate=False, CLIP=0.601, RNA=0.094, TE=0.663
- Cdh1: candidate=False, CLIP=1.623, RNA=0.085, TE=0.700

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
