# Week 5 Report: CLIP enrichment versus RNA and TE changes

## Goal

Use the week4 metric table to test whether LIN28A binding is linked more
strongly to translation efficiency change than to RNA abundance change.

## Main result

Among base-filtered genes, CLIP enrichment showed only a weak negative
relationship with RNA abundance change:

- Pearson r = -0.272
- Spearman r = -0.295

By contrast, CLIP enrichment was positively associated with TE change:

- Pearson r = 0.488
- Spearman r = 0.481

This matches the paper's central model: LIN28A-bound transcripts are affected
mainly at the translation level after Lin28a knockdown.

## CLIP group comparison

The top 5% CLIP-enriched genes had median TE_log2FC = 0.539,
whereas the bottom 50% had median TE_log2FC = -0.513.
For top 20% versus bottom 50%, the Mann-Whitney p-value for TE_log2FC was
< 1e-300.

## Robustness

Across RNA count thresholds 10, 30, 50 and RPF-siLuc thresholds 50, 80, 100,
the Spearman correlation between CLIP enrichment and TE_log2FC ranged from
0.462 to 0.487. The conclusion is therefore not driven
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
