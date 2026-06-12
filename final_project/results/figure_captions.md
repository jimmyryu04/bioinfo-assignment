# Figure captions

> 보고서 본문에 그대로 붙여 쓸 수 있는 캡션. 모든 figure는 제목(title) 없이 출력되며, 패널은 A/B/C/D 로 식별하고 내용은 아래 캡션이 설명한다. 수치는 `results/tables/`, `results/result_summary.md` 와 일치한다.

---

## Figure 1. LIN28A binding amount correlates with translational derepression after *Lin28a* knockdown.

Each point is a gene (n = 8,986 baseline-eligible genes). The x-axis is transcript-level LIN28A CLIP enrichment, log₂(CPM_CLIP / CPM_RNA-control) from featureCounts. The y-axis is the change in ribosome density after *Lin28a* knockdown, log₂[(RPF/RNA)_siLin28a / (RPF/RNA)_siLuc]. The red line is an ordinary least-squares fit. Binding and derepression are positively correlated (Pearson r = 0.464, p < 1×10⁻³⁰⁰; Spearman ρ = 0.458), reproducing the transcript-level relationship of Cho et al. (2012), Fig. 4D.

*KR: 유전자별 총 LIN28A 결합량(x)과 knockdown 후 번역 탈억제(y)의 양의 상관 — 원 논문 Fig 4D 재현(baseline 검증).*

---

## Figure 2. Defining an RNA-normalized regional binding composition.

**(A)** Distribution of 5′ UTR, CDS, and 3′ UTR lengths (log₁₀) of the representative transcripts used for region assignment. **(B)** RNA-control-adjusted binding composition per region for composition-eligible genes (n = 5,723); for each gene the regional enrichment E_region = CLIP/RNA-control is normalized to a fraction across the three regions. Dashed lines are quartiles. **(C)** Mean composition computed from raw CLIP 5′-end point fractions versus the RNA-adjusted fractions; the raw fraction tracks region length/coverage (CDS-dominated), whereas the RNA-adjusted fraction is the predictor used in the model. **(D)** Gene-count filtering flow from all annotated genes to the composition-eligible model set.

*KR: (A) 대표 transcript의 region 길이, (B) RNA 보정 결합 구성비, (C) raw vs RNA 보정 구성비(보정이 왜 필요한지), (D) 필터링 단계별 유전자 수.*

---

## Figure 3. Total binding amount, not regional position, drives translational derepression.

Nested OLS models fit on the same composition-eligible genes (n = 5,723). **(A)** Variance explained (R²) by M0 (covariates: log RNA expression, log transcript length), M1 (M0 + total CLIP amount), and M2 (M1 + CDS and 3′ UTR composition). **(B)** Incremental ΔR²: total binding amount adds 0.214 over covariates, while regional composition adds only 0.008 (nested F-test M1 vs M2, F = 32.21, p = 1.2×10⁻¹⁴ — significant but small). **(C)** M2 standardized β values after standardizing both predictors and response: total CLIP amount +0.512, CDS fraction −0.097, 3′ UTR fraction −0.009 (5′ UTR is the reference). Amount is the dominant predictor; positional composition contributes a small, statistically significant increment in which CDS-concentrated binding is associated with slightly *less* derepression.

*KR: 같은 유전자 집합의 nested model — (A) R², (B) 증분 ΔR²(amount 0.214 vs position 0.008), (C) M2 표준화 계수. 위치보다 총 결합량이 주효과.*

---

## Figure 4. Higher total binding shows a weak positive dose relationship with derepression.

Composition-eligible genes binned into quartiles (Q1–Q4) by total exon-assigned CLIP 5′-end point count (n per group shown on the x-axis). Violins show the distribution of Δ ribosome density, white boxes the interquartile range and median, and the red line connects group medians. The dose trend is weak but significant (Spearman ρ = 0.066, p = 5.4×10⁻⁷), consistent with the small position effect in Figure 3 and with binding burden contributing modestly to translational repression.

*KR: CLIP 5′ point 수 분위수별 Δribosome density — 약하지만 유의한 dose 추세(rho=0.066).*

---

### 본문 연결 한 줄
원 논문(Fig 4D)의 transcript-level 결합–번역 관계를 재현(Fig 1)한 뒤, 결합을 총량과 RNA 보정 위치 구성으로 분해(Fig 2)하여 nested model로 검정한 결과(Fig 3), 번역 탈억제는 **결합 위치보다 총 결합량**으로 더 잘 설명되며 dose 관계도 약하게 존재한다(Fig 4).
