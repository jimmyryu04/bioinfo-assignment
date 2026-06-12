# Pipeline result summary

## 실행 산출물

- `results/tables/`: gene-level response, region counts, model summaries
- `results/work/`: BED/intermediate count files
- `results/figures/`: final report figures

## QC

- Actb strandedness check: same-strand fraction = 0.9996. 이 값이 높으면 sense-stranded counting(`-s 1`, `bedtools -s`)이 타당하다.
- CLIP alignments seen = 21,877,250; NH=1 5-prime points written = 18,802,078.
- Composition-eligible model set n = 5,723. M0/M1/M2는 동일 n에서 적합했다.

## Baseline 재현

- baseline featureCounts CLIP enrichment vs DeltaRD: n = 8,986, Pearson r = 0.4637, p = <1e-300; Spearman rho = 0.4578.
- 이 그림은 원 논문/W1의 transcript-level total binding signal이 knockdown 후 ribosome density change와 연결되는지 확인하는 QC 역할이다.

## Nested model 결과

- M0 R2 = 0.0701
- M1 R2 = 0.2843; DeltaR2_amount = 0.2142; F-test p = <1e-300
- M2 R2 = 0.2923; DeltaR2_position = 0.0080; F-test p = 1.2e-14
- M2 standardized beta: total_CLIP_amount_point = 0.5117, f_cds = -0.0968, f_3utr = -0.0089.
- 방향 해석: M2에서 `f_cds` 계수는 음수이므로, total amount를 통제한 뒤 CDS 쪽 adjusted composition이 높을수록 derepression이 강해진다는 방향은 아니다. 위치 효과는 작고 보조적이며, 특히 CDS-positive story로 과장하면 안 된다.
- 해석: total CLIP amount는 covariate-only model보다 유의하게 DeltaRD를 더 설명했다. position composition은 통계적으로 유의하지만 증분 설명력은 작다. 따라서 주효과는 total amount이고, 위치는 보조적인 신호로 해석하는 것이 안전하다.

## Dose 분석

- total exon-assigned CLIP 5-prime point count vs DeltaRD Spearman rho = 0.0662, p = 5.4e-07.
- dose quartile별 ribosome density 변화 분포는 `fig4_dose_response.png`와 `tables/dose_summary.tsv`에서 확인한다.

## 변수별 상관 요약

| variable                            | n    | pearson_r           | pearson_p              | spearman_rho        | spearman_p              |
| ----------------------------------- | ---- | ------------------- | ---------------------- | ------------------- | ----------------------- |
| total_CLIP_enrichment_featureCounts | 5723 | 0.4713498910837573  | 1.675254159e-314       | 0.453457220491219   | 2.465561316873793e-288  |
| total_CLIP_amount_point             | 5723 | 0.4138658985127801  | 1.009392593505742e-235 | 0.4201961073451586  | 1.1413175781070845e-243 |
| clip_points_total                   | 5723 | 0.0670032391634288  | 3.907207888110004e-07  | 0.0661821686166647  | 5.411537211403647e-07   |
| f_5utr                              | 5723 | -0.096808485890197  | 2.147861499529035e-13  | -0.1022681192559923 | 8.814473171398101e-15   |
| f_cds                               | 5723 | -0.0395600654748243 | 0.0027601707800168     | -0.0403956140964696 | 0.002239159031478       |
| f_3utr                              | 5723 | 0.1043533782725914  | 2.486938246158333e-15  | 0.0987592231855894  | 7.000334438439961e-14   |

## 표준화 회귀계수

| model | term                    | effect_per_1sd_predictor_delta_rd | standardized_beta   | n    |
| ----- | ----------------------- | --------------------------------- | ------------------- | ---- |
| M2    | intercept               | -0.2763308476293198               |                     | 5723 |
| M2    | log_rna_expression      | -0.1251601997726781               | -0.2226062973326539 | 5723 |
| M2    | log_transcript_length   | -0.1256565497583016               | -0.223489091005734  | 5723 |
| M2    | total_CLIP_amount_point | 0.2877273603516839                | 0.511743529057297   | 5723 |
| M2    | f_cds                   | -0.0544155374878837               | -0.0967818950396745 | 5723 |
| M2    | f_3utr                  | -0.0049929178097547               | -0.0088802586487913 | 5723 |

## Eligibility flow

| stage                          | n     |
| ------------------------------ | ----- |
| gene_response_rows             | 55359 |
| response_eligible              | 9149  |
| representative_region_matched  | 21743 |
| has_all_regions                | 20465 |
| clip_points_total_ge_10        | 12641 |
| rna_each_region_ge_5           | 6230  |
| composition_eligible_model_set | 5723  |

## 결론 가이드

- 유의미한 분석인가: 같은 gene set에서 nested comparison을 수행했고, amount와 RNA-control adjusted position composition을 같은 stranded point signal에서 유도했으므로 설계상 주제 질문에 직접 답한다.
- 단, 조건당 생물학적 반복이 없으므로 gene별 differential translation 유의성은 주장하지 않는다. 결과는 across-gene association이다.
- `DeltaR2_position`이 작거나 비유의적이면 negative result도 의미가 있다. 이는 LIN28A의 번역 억제 효과가 region-specific preference보다 전체 결합량에 더 잘 설명될 수 있음을 뜻한다.
