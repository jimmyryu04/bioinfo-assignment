# Final project pipeline

This folder contains the MVP analysis code for:

> LIN28A binding amount and regional position composition as predictors of translational derepression after Lin28a knockdown.

Run from the repository/workspace root:

```bash
python etc/bioinformatics/final_project/pipeline/run_pipeline.py
```

To regenerate only figures from existing tables/intermediate outputs:

```bash
cd etc/bioinformatics/final_project/pipeline
python regen_figures.py
```

Outputs are written to:

```text
etc/bioinformatics/final_project/results/
  figures/     # PNG and PDF figure outputs
  tables/
  work/
  logs/
  result_summary.md
```

The implementation intentionally uses only packages available in the current environment. OLS and nested F-tests are computed with `numpy`/`scipy`; BAM parsing is done with `samtools view` plus a CIGAR-aware Python parser.
