# Amount- and position-aware analysis of LIN28A binding as a predictor of translational derepression after *Lin28a* knockdown

Dongmin Ryu (2025-25894)

Bioinformatics Final Project

**LIN28A binds thousands of mature mRNAs and has been proposed to suppress translation of targets associated with the endoplasmic reticulum. Here I reanalyse the Cho et al. LIN28A CLIP-seq and ribosome-profiling data to ask whether translational derepression after *Lin28a* knockdown is explained primarily by total LIN28A binding amount or by the regional composition of binding across the 5-prime UTR, coding sequence and 3-prime UTR. The analysis recapitulates the original transcript-level association between LIN28A CLIP enrichment and ribosome-density change, then decomposes binding into a stranded, RNA-control-adjusted point-count amount signal and a regional composition signal. Total binding amount is the dominant predictor: adding total CLIP amount to an expression- and transcript-length model increases explained variance by DeltaR2 = 0.2142. Regional composition is statistically detectable but small, adding only DeltaR2 = 0.0080. These results support an amount-dominant interpretation of LIN28A-linked translational repression and do not support a strong position-specific mechanism in this dataset.**

LIN28 proteins are best known for repressing let-7 microRNA maturation: LIN28 blocks Microprocessor processing of pri-let-7 and promotes terminal uridylation of pre-let-7 [1,2]. This canonical axis does not exhaust LIN28 function. CLIP-based work showed that LIN28 also binds messenger RNAs at GGAGA-containing motifs and regulates splicing-factor abundance [3]. In mouse embryonic stem cells, Cho and colleagues then combined LIN28A CLIP-seq with ribosome footprinting to show that LIN28A binds many mature mRNAs and that LIN28A-bound transcripts tend to show increased ribosome occupancy after *Lin28a* knockdown [4]. This result reframed LIN28A as a regulator of mRNA translation, particularly for transcripts encoding endoplasmic-reticulum-associated, membrane and secretory proteins. The original study established a strong transcript-level connection between LIN28A binding and translation, but it left a more architectural question open: does the translational response depend mostly on how much LIN28A binds a transcript, or does the location of binding within the transcript provide additional explanatory information?

This question matters because total CLIP enrichment and binding position are not equivalent biological quantities. A high CLIP signal could reflect many binding events distributed across a transcript, local enrichment in a particular region, increased RNA abundance or region length, or some combination of these features. A transcript with relatively CDS-enriched binding might suggest a connection to elongation or ribosome traffic, whereas a transcript with more 5-prime UTR or 3-prime UTR enrichment could suggest effects on initiation, mRNA localization or regulatory-complex recruitment. Those interpretations are plausible only if regional binding preference is separated from total binding amount and from region-specific RNA coverage. I therefore framed the final project around a single testable question: after controlling for total LIN28A binding amount, does RNA-control-adjusted regional binding composition add explanatory power for translational derepression after *Lin28a* knockdown?

The analysis was deliberately designed as a conservative association study. The local datapack contains one library per condition rather than biological replicates, so it cannot support gene-level differential translation tests or uncertainty estimates for individual genes. All inferential statements below refer to across-gene associations: genes are data points, DeltaRD summarizes the change in ribosome density after knockdown, and model p values describe regression comparisons, not replicate-based differential translation. This restriction is central to the interpretation of every result.

## LIN28A binding amount recapitulates translational derepression

As a baseline check, I first reproduced the transcript-level relationship between LIN28A binding and translational derepression using the week 1 featureCounts-derived gene table. The baseline predictor was total LIN28A CLIP enrichment, defined as log2(CPM_CLIP / CPM_RNA-control), and the response was DeltaRD, defined as log2[(RPF/RNA)_siLin28a / (RPF/RNA)_siLuc]. Among 8,986 baseline-eligible genes, total CLIP enrichment was positively associated with DeltaRD (Pearson r = 0.4637, Spearman rho = 0.4578, p < 1e-300; Fig. 1). This result does not establish causality, but it confirms that the working analysis table retains the central signal reported by Cho et al. [4]: LIN28A-bound transcripts tend to show greater ribosome-density increase when LIN28A is depleted.

This baseline association is important for two reasons. First, it links the final project directly to the original paper rather than creating an unrelated analysis on the same files. Second, it provides a quality-control gate for the subsequent decomposition of LIN28A binding architecture. If the baseline CLIP enrichment failed to associate with DeltaRD, then any amount-versus-position analysis would be difficult to interpret. Instead, the baseline result was strong and directionally consistent, justifying a more detailed model in which total binding and regional binding composition were handled separately.

![Fig. 1](results/figures/fig1_baseline_total_clip_vs_delta_rd.pdf)

**Fig. 1 | LIN28A binding amount correlates with translational derepression after *Lin28a* knockdown.** Each point is a gene (n = 8,986 baseline-eligible genes). The x axis shows transcript-level LIN28A CLIP enrichment, log2(CPM_CLIP / CPM_RNA-control), from featureCounts. The y axis shows DeltaRD, log2[(RPF/RNA)_siLin28a / (RPF/RNA)_siLuc]. The red line is an ordinary least-squares fit. Binding and derepression are positively correlated (Pearson r = 0.4637, p < 1e-300; Spearman rho = 0.4578), recapitulating the transcript-level relationship in Cho et al. [4].

## Regional composition was defined from RNA-normalized point signals

The central methodological step was to construct regional predictors that were not merely raw regional read counts. I selected one representative protein-coding transcript per gene, derived strand-aware 5-prime UTR, CDS and 3-prime UTR intervals from the GENCODE annotation, and counted LIN28A CLIP 5-prime-end points within those regions using sense-stranded overlap. The CLIP signal was restricted to NH=1 alignments to reduce ambiguity from multimapping reads. RNA-control reads were counted over the same transcript regions and used as the regional denominator.

For each gene and region, the regional enrichment was defined as E_region = (CLIP_region + pseudocount) / (RNA-control_region + pseudocount), using CPM-scaled point and RNA-control signals. The three regional enrichments were then normalized to fractions, f_5utr, f_cds and f_3utr, which sum to one within each gene. This construction makes the position predictor a relative regional preference rather than a second measurement of total CLIP burden. Because the three fractions are compositional, the regression model used f_cds and f_3utr with f_5utr as the implicit reference.

Fig. 2 shows why this adjustment is necessary. Raw CLIP 5-prime point fractions are dominated by CDS signal, which is expected because CDS intervals are long and often better covered. After RNA-control adjustment, the mean composition shifts, with a larger apparent contribution from the 3-prime UTR and a small but measurable 5-prime UTR component. The final model set contained 5,723 genes that passed response filters, had a matched representative transcript with all three regions, contained at least ten exon-assigned CLIP points and had at least five RNA-control reads in each region.

![Fig. 2](results/figures/fig2_region_composition_summary.pdf)

**Fig. 2 | Defining an RNA-normalized regional binding composition.** **a**, Distribution of 5-prime UTR, CDS and 3-prime UTR lengths for representative transcripts. **b**, RNA-control-adjusted binding composition for composition-eligible genes (n = 5,723), where regional enrichment is normalized to a fraction across the three transcript regions. **c**, Mean raw CLIP 5-prime point fractions compared with RNA-adjusted fractions, illustrating why regional RNA-control normalization is needed. **d**, Filtering flow from all count-table genes to the final model set.

## Binding amount dominated the nested regression models

The primary statistical test compared three nested ordinary least-squares models fit on exactly the same 5,723 composition-eligible genes. The covariate model M0 used log RNA expression and log transcript length. The amount model M1 added total_CLIP_amount_point, the log2 ratio of total exon-assigned CLIP 5-prime points to total RNA-control reads. The position model M2 added f_cds and f_3utr to M1. The primary endpoint was DeltaR2_position = R2(M2) - R2(M1).

The covariate-only model explained a modest fraction of DeltaRD variation (M0 R2 = 0.0701). Adding total binding amount had a large effect, increasing R2 to 0.2843 (DeltaR2_amount = 0.2142; Fig. 3a,b). Adding regional composition increased R2 further to 0.2923, but the increment was much smaller (DeltaR2_position = 0.0080). The M1-versus-M2 nested F-test was statistically significant (p = 1.2e-14), as expected with thousands of genes, but the magnitude of the improvement was small. Thus, the model comparison supports the core claim that total LIN28A binding amount, not regional position composition, carries most of the explanatory signal for translational derepression in this analysis.

The coefficient pattern reinforced this conclusion. On the predictor-standardized DeltaRD scale recorded in the coefficient table, total_CLIP_amount_point had a coefficient of 0.2877, whereas f_cds and f_3utr had coefficients of -0.0544 and -0.0050, respectively. The corresponding fully standardized betas plotted in Fig. 3c, after also scaling DeltaRD, were 0.5117, -0.0968 and -0.0089. Both coefficient scales give the same interpretation: once total binding amount is controlled, higher CDS-adjusted composition is not associated with stronger derepression. The negative f_cds coefficient argues against a simple CDS-positive mechanistic story. A stronger version of the position hypothesis would have predicted a sizeable positive regional term, especially if CDS-localized LIN28A binding were directly linked to ribosome movement. The observed data instead indicate a small residual compositional signal whose biological interpretation should remain cautious.

![Fig. 3](results/figures/fig3_nested_model.pdf)

**Fig. 3 | Total binding amount is the larger predictor than regional position.** Nested OLS models were fit on the same 5,723 composition-eligible genes. **a**, Variance explained by M0 (log RNA expression and log transcript length), M1 (M0 plus total CLIP amount) and M2 (M1 plus CDS and 3-prime UTR composition). **b**, Incremental variance explained: total binding amount adds DeltaR2 = 0.2142, whereas regional composition adds DeltaR2 = 0.0080 (M1 versus M2 nested F-test, p = 1.2e-14). **c**, Fully standardized M2 coefficients, showing that total binding is the largest term and that the CDS composition term is not positive.

## Binding dose showed only a weak trend

Because the position-composition signal was small, I also examined whether a simpler binding-dose measure retained any association with DeltaRD. Genes in the composition-eligible set were binned into quartiles by total exon-assigned CLIP 5-prime point count. Median DeltaRD increased weakly across the dose bins, and the overall Spearman correlation between total point count and DeltaRD was positive but very small (rho = 0.0662, p = 5.4e-07; Fig. 4). The large sample size makes this trend statistically detectable, but the effect is weak in magnitude and should not be interpreted as a strong dose-response curve.

The dose result is nevertheless useful because it agrees with the nested-model interpretation. The most informative amount variable was not raw point count but RNA-control-normalized total CLIP amount. Raw binding count alone was only weakly associated with response, suggesting that normalization and expression context matter. Taken together, the baseline enrichment, nested regression and dose analysis converge on a measured conclusion: binding amount contains the main signal, but no single simple count or regional fraction is sufficient to explain the response.

![Fig. 4](results/figures/fig4_dose_response.pdf)

**Fig. 4 | Higher total binding shows a weak positive dose relationship with derepression.** Composition-eligible genes were binned into quartiles by total exon-assigned CLIP 5-prime point count. Violins show DeltaRD distributions, white boxes show the interquartile range and median, and the red line connects group medians. The dose trend is weak but statistically detectable (Spearman rho = 0.0662, p = 5.4e-07).

## Discussion and limitations

The main result is therefore not a claim that transcript position is irrelevant. Regional composition added a statistically significant increment to the model, and the composition construction likely captures real variation in the distribution of LIN28A binding across transcript regions. The result is instead a statement about scale. In these data, under this course-project design, the additional positional signal is much smaller than the total binding-amount signal. The data support an amount-dominant model in which transcripts with more RNA-normalized LIN28A binding are more likely to show translational derepression after *Lin28a* knockdown, while regional binding composition contributes only a small residual association.

This conclusion extends the original study rather than reproducing a single figure. Cho et al. showed that LIN28A binds mature mRNAs, that average CLIP density varies across transcript regions, and that transcript-level binding is associated with ribosome-density change [4]. The present analysis connects those axes directly by asking whether the positional architecture of binding adds explanatory information beyond total binding. The answer is qualified: yes in a statistical sense, but only weakly in explanatory magnitude. This makes the project useful as a mechanistic narrowing exercise. It argues against over-interpreting region-specific coefficients and points instead toward broader binding burden, RNA context, ER-associated localization or other transcript-level properties as more plausible drivers of the observed derepression.

Several limitations bound this interpretation. The available datapack lacks biological replicates, so the analysis cannot estimate condition-level variance or call individual genes differentially translated. DeltaRD was computed from the featureCounts gene-level table and is therefore an approximation to the CDS-limited ribosome-density metric used in the original study. CLIP 5-prime-end points were used as a practical point-signal proxy rather than genome-wide CIMS calls, even though CIMS-based site calling is more directly tied to crosslink positions at nucleotide resolution [7]. The position effect is small, and the dose trend is weak. These caveats do not invalidate the amount-dominant conclusion, but they define it: the result is an across-gene association analysis, not a causal or replicate-supported mechanism for a specific transcript region.

## Main references

1. Viswanathan, S. R., Daley, G. Q. & Gregory, R. I. Selective blockade of microRNA processing by Lin28. *Science* **320**, 97-100 (2008). https://doi.org/10.1126/science.1154040

2. Heo, I. *et al.* Lin28 mediates the terminal uridylation of let-7 precursor MicroRNA. *Mol. Cell* **32**, 276-284 (2008). https://doi.org/10.1016/j.molcel.2008.09.014

3. Wilbert, M. L. *et al.* LIN28 binds messenger RNAs at GGAGA motifs and regulates splicing factor abundance. *Mol. Cell* **48**, 195-206 (2012). https://doi.org/10.1016/j.molcel.2012.08.004

4. Cho, J. *et al.* LIN28A is a suppressor of ER-associated translation in embryonic stem cells. *Cell* **151**, 765-777 (2012). https://doi.org/10.1016/j.cell.2012.10.019

5. Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S. & Weissman, J. S. Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling. *Science* **324**, 218-223 (2009). https://doi.org/10.1126/science.1168978

6. Licatalosi, D. D. *et al.* HITS-CLIP yields genome-wide insights into brain alternative RNA processing. *Nature* **456**, 464-469 (2008). https://doi.org/10.1038/nature07488

7. Zhang, C. & Darnell, R. B. Mapping in vivo protein-RNA interactions at single-nucleotide resolution from HITS-CLIP data. *Nat. Biotechnol.* **29**, 607-614 (2011). https://doi.org/10.1038/nbt.1873

## Methods

### Data source and study context

This project reanalysed the local final-project datapack derived from Cho et al. [4], which contains LIN28A CLIP-seq, RNA-seq and ribosome-footprinting libraries from mouse embryonic stem cells. Ribosome profiling was originally introduced as a nucleotide-resolution method for measuring translation in vivo [5], and HITS-CLIP established a sequencing-based framework for mapping protein-RNA interactions in living tissue [6]. The original Cho study reported that LIN28A binds many mature mRNAs and that LIN28A-bound transcripts show increased ribosome occupancy after *Lin28a* knockdown. The present analysis used the local BAM files, the GENCODE GTF annotation and the week 1 featureCounts table supplied with the course project. No external data were downloaded for the course-project analysis.

### Gene-level response calculation

Gene-level counts were read from `week1/work/read-counts.txt`. Counts were converted to CPM using library-size normalization. Baseline total CLIP enrichment was calculated as log2((CPM_CLIP-35L33G + pc) / (CPM_RNA-control + pc)). DeltaRD was calculated as log2(((CPM_RPF-siLin28a + pc) / (CPM_RNA-siLin28a + pc)) / ((CPM_RPF-siLuc + pc) / (CPM_RNA-siLuc + pc))). The pseudocount was 0.1 CPM. Genes used in the response analysis were filtered by minimum RNA/RPF counts as described in the pipeline. The baseline figure used 8,986 baseline-eligible genes.

### Representative transcript selection

To connect gene-level DeltaRD values to transcript-level region annotation, a single representative protein-coding transcript was selected for each gene. The selection prioritized protein-coding genes and transcripts, transcript support level 1, APPRIS principal, CCDS and basic tags, then longer CDS and transcript lengths. This avoided assigning multiple transcript annotations to one gene-level response value.

### UTR and CDS region annotation

CDS and UTR features were extracted from the GTF for each representative transcript. Because the annotation did not provide separate 5-prime and 3-prime UTR features in the required form, UTR intervals were classified relative to the strand-aware CDS boundaries. For plus-strand transcripts, UTR intervals upstream of the CDS were assigned to the 5-prime UTR and downstream intervals to the 3-prime UTR. For minus-strand transcripts, the assignment was reversed. Ambiguous UTR intervals overlapping the CDS boundary were not used in the regional count table.

### CLIP 5-prime-end point extraction

LIN28A CLIP signal was represented by strand-aware 5-prime-end points from `CLIP-35L33G.bam`. Primary, non-secondary and non-supplementary alignments were extracted with `samtools view`, and only NH=1 alignments were retained. For forward reads, the point was the leftmost reference coordinate. For reverse reads, the point was the rightmost reference-consuming coordinate, computed with a CIGAR-aware parser. The pipeline observed 21,877,250 CLIP alignments and wrote 18,802,078 NH=1 5-prime points. Actb strandedness was used as a sanity check; the same-strand fraction was 0.9996, supporting sense-stranded counting. This point representation is simpler than CIMS-based CLIP analysis, which uses crosslink-induced mutation sites to localize protein-RNA contacts at single-nucleotide resolution [7].

### RNA-control-adjusted regional composition

CLIP points and RNA-control reads were counted over the representative transcript 5-prime UTR, CDS and 3-prime UTR intervals using same-strand overlap. Region-specific CLIP and RNA-control signals were CPM-normalized. For each region, enrichment was calculated as E_region = (CLIP_region_CPM + pc) / (RNA-control_region_CPM + pc). Regional composition fractions were then defined as f_region = E_region / (E_5utr + E_cds + E_3utr). These fractions describe relative regional preference after RNA-control adjustment rather than raw point allocation.

### Eligibility filtering

The composition model used genes that were response eligible, matched to a representative transcript, had all three annotated regions, had at least ten total exon-assigned CLIP 5-prime points and had at least five RNA-control reads in each region. The filtering flow was: 55,359 gene-response rows; 9,149 response-eligible genes; 21,743 genes matched to representative transcript regions; 20,465 genes with all three regions; 12,641 genes with at least ten total CLIP points; 6,230 genes passing the RNA-control threshold in each region; and 5,723 genes in the final composition-eligible model set.

### Nested OLS models

All nested models were fit on the same 5,723 composition-eligible genes. M0 used log RNA expression and log transcript length. M1 added total_CLIP_amount_point, defined as the log2 ratio of total exon-assigned CLIP point CPM to total RNA-control CPM. M2 added f_cds and f_3utr, leaving f_5utr as the implicit reference because the three composition fractions sum to one. OLS models were fit with NumPy. Nested F tests were computed from residual sums of squares and evaluated with the SciPy F survival function.

### Association and dose-response analysis

Pearson and Spearman correlations were calculated for baseline CLIP enrichment, point-based total CLIP amount, raw point count and regional fractions against DeltaRD. For the dose analysis, genes were binned into quartiles by total exon-assigned CLIP 5-prime point count. DeltaRD distributions were shown with violin and box plots, and the overall Spearman correlation between point count and DeltaRD was reported.

### Software and reproducibility

The analysis pipeline is implemented in `pipeline/run_pipeline.py`, with figure regeneration supported by `pipeline/regen_figures.py`. Outputs are written to `results/tables`, `results/work`, `results/figures` and `results/result_summary.md`. The current manuscript was built from the result tables and figure PDFs in the project directory. The standard manuscript compilers requested for this report (`typst`, `pandoc`, `xelatex`, `lualatex`, `pdflatex` and `tectonic`) were not available on the current `PATH`; therefore, a Typst source file is provided for future compilation, and the submitted PDF in this environment was generated by a local PyMuPDF fallback builder.

## Data availability

Cho et al. LIN28A CLIP-seq and ribosome-profiling study [4]. The manuscript does not introduce new external datasets.
