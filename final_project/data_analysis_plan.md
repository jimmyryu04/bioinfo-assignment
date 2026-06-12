# Data analysis plan

## Project title

**LIN28A 결합의 총량과 위치 구성이 Lin28a knockdown 후 번역 탈억제를 얼마나 설명하는가**

English working title:

> Amount- and position-aware analysis of LIN28A binding as a predictor of translational derepression after Lin28a knockdown

## 1. 분석 목표

원 논문 Cho et al. (2012)은 LIN28A의 transcript-level CLIP enrichment가 Lin28a knockdown 후 ribosome density 증가와 연결된다는 것을 보였다. 본 분석은 이 결과를 확장해, LIN28A 결합을 다음 두 축으로 분리한다.

1. **Amount**: 한 gene/transcript에 LIN28A가 얼마나 많이 결합하는가?
2. **Position composition**: 같은 총 결합량이라도, LIN28A 결합이 5-prime UTR, CDS, 3-prime UTR 중 어디에 상대적으로 농축되어 있는가?

핵심 검정은 다음이다.

```text
M0: DeltaRD ~ expression + gene_length
M1: DeltaRD ~ expression + gene_length + total_CLIP_amount
M2: DeltaRD ~ expression + gene_length + total_CLIP_amount + position_composition
```

주요 endpoint:

```text
Delta R2_position = R2(M2) - R2(M1)
```

즉, total CLIP amount를 통제한 뒤에도 위치 구성이 Lin28a knockdown 후 번역 탈억제를 추가로 설명하는지 평가한다.

## 2. 입력 데이터

사용 데이터는 모두 로컬 `final_project` 폴더 안에 있다.

| 데이터 | 경로 | 용도 |
|---|---|---|
| LIN28A CLIP-seq | `week1/work/binfo1-datapack1/CLIP-35L33G.bam` | LIN28A 결합 amount/position 계산 |
| RNA control | `week1/work/binfo1-datapack1/RNA-control.bam` | region별 RNA coverage 보정, total CLIP enrichment denominator |
| RNA siLuc | `week1/work/binfo1-datapack1/RNA-siLuc.bam` | control condition mRNA abundance |
| RNA siLin28a | `week1/work/binfo1-datapack1/RNA-siLin28a.bam` | Lin28a knockdown mRNA abundance |
| RPF siLuc | `week1/work/binfo1-datapack1/RPF-siLuc.bam` | control condition ribosome footprint |
| RPF siLin28a | `week1/work/binfo1-datapack1/RPF-siLin28a.bam` | Lin28a knockdown ribosome footprint |
| GENCODE GTF | `week1/work/binfo1-datapack1/gencode.gtf` | representative transcript, UTR/CDS annotation |
| featureCounts table | `week1/work/read-counts.txt` | gene-level Delta ribosome density와 baseline CLIP enrichment |

주의:

- 이 datapack은 조건당 생물학적 반복이 없다. 따라서 gene별 differential translation p-value가 아니라, gene을 데이터포인트로 한 across-gene association 분석으로 해석한다.
- 원 논문 ribosome density는 CDS-limited RPF/RNA count로 계산했다. MVP에서는 `read-counts.txt`의 gene-level count를 사용하므로 원 논문 정의의 근사값이다. 시간이 남으면 representative CDS에 대해 RPF/RNA를 재카운트하는 sensitivity analysis를 추가한다.

## 3. 사용할 도구

| 도구 | 경로/패키지 | 용도 |
|---|---|---|
| `samtools` | `/blaze/ryu/conda/envs/lab/bin/samtools` | BAM filtering, indexing, region extraction |
| `bedtools` | `/blaze/ryu/conda/envs/lab/bin/bedtools` | BED intersection, region count |
| `featureCounts` | `/blaze/ryu/conda/envs/lab/bin/featureCounts` | 필요 시 CDS/region count 재계산 |
| Python `pandas`, `numpy` | installed | table parsing, normalization |
| Python `scipy` | installed | Pearson/Spearman, F-test p-value, Mann-Whitney |
| Python `sklearn` | installed | linear regression / model fitting helper |
| Python `matplotlib`, `seaborn` | installed | figure generation |
| Python `statsmodels` | optional | 설치되어 있으면 OLS/anova_lm 사용, 없으면 `numpy`/`scipy`로 F-test 직접 계산 |

권장 실행 환경:

```bash
export MPLCONFIGDIR=/tmp/ryu/matplotlib
```

## 4. 결과 디렉토리 구조

분석은 새 작업 폴더에서 진행한다.

```text
your_own_analysis/
  scripts/
    01_gene_response.py
    02_build_representative_regions.py
    03_count_clip_regions.py
    04_model_and_plot.py
  work/
    representative_transcripts.tsv
    regions.5utr.bed
    regions.cds.bed
    regions.3utr.bed
    clip_5end_points.bed
    region_counts.tsv
    analysis_table.tsv
  results/
    fig1_baseline_total_clip_vs_delta_rd.png
    fig2_region_composition_summary.png
    fig3_nested_model.png
    fig4_dose_response.png
  tables/
    gene_response.tsv
    model_summary.tsv
    eligibility_summary.tsv
```

## 5. MVP 분석 흐름

보고서 본문은 MVP만으로 완결되도록 한다. CIMS, permutation, CDS-limited re-counting 등은 시간이 남을 때 supplementary로 수행한다.

실행 전 반드시 확인할 구현 조건:

- 이 datapack의 CLIP/RNA/RPF library는 strand-specific로 보인다. region count는 `bedtools intersect -s` 또는 `featureCounts -s 1`에 해당하는 **sense-stranded** 방식으로 수행한다.
- 기존 Actb locus 확인에서는 read strand가 gene strand와 일치해 `-s 1`이 맞을 가능성이 높다. 구현 직전 Actb 같은 known gene에서 `-s 1` vs `-s 2`를 한 번 더 재확인하고, 반대이면 `-S`/`-s 2`로 바꾼다.
- `read-counts.txt`는 week1에서 `featureCounts`를 `-s` 없이 실행한 unstranded baseline이다. 따라서 이 파일은 Fig 4D baseline 재현과 gene-level DeltaRD 근사에 사용하되, model용 LIN28A amount/position predictor는 새로 만든 stranded point count에서 유도한다.
- `read-counts.txt`의 Ensembl gene ID는 version이 포함되어 있다. GTF-derived table과 merge하기 전에 모든 gene ID에서 version suffix, 예: `.2`, `.15`, 를 제거한다.

### Step 1. gene-level response 계산

입력:

- `week1/work/read-counts.txt`

계산:

```text
CPM = raw_count / library_size * 1e6

total_CLIP_enrichment =
log2((CPM_CLIP-35L33G + pc) / (CPM_RNA-control + pc))

DeltaRD =
log2((CPM_RPF-siLin28a / CPM_RNA-siLin28a) /
     (CPM_RPF-siLuc    / CPM_RNA-siLuc))
```

필터:

- 기본: RNA/RPF/CLIP count가 너무 낮은 gene 제거
- 시작값: week1과 맞추기 위해 all sample count >= 30 또는 RPF-siLuc >= 80 기준을 비교
- 최종 보고서에는 필터 통과 gene 수를 명시

산출물:

- `tables/gene_response.tsv`

주요 컬럼:

```text
gene_id
gene_name
gene_length
clip_cpm
rna_control_cpm
total_CLIP_enrichment
delta_rd
rna_expression
```

시각화:

- `fig1_baseline_total_clip_vs_delta_rd.png`
- x축: total CLIP enrichment
- y축: Delta ribosome density
- 표시: Pearson/Spearman r, regression line, n

보여주려는 것:

- 원 논문의 Fig 4D/W1 결과가 이 분석 테이블에서도 재현되는지 확인한다.
- 기존 featureCounts 기반 total CLIP enrichment가 실제로 translation derepression과 연결되는지 검증한다.
- 단, 이 total CLIP enrichment는 unstranded full-read gene count에서 온 baseline 검증용 값이다. M1/M2 model에 들어가는 `total_CLIP_amount`는 Step 4-5의 stranded point signal에서 새로 계산한다.

### Step 2. representative transcript set 구축

입력:

- `gencode.gtf`

목표:

gene-level DeltaRD와 transcript-level UTR/CDS annotation을 1:1로 연결하기 위해 gene당 대표 transcript 1개를 선택한다.

대표 transcript 선택 우선순위:

1. `gene_type == protein_coding`
2. `transcript_type == protein_coding`
3. `transcript_support_level == "1"` 우선
4. `tag`에 `basic`, `appris_principal`, `CCDS`가 있으면 우선
5. CDS 길이가 가장 긴 transcript
6. tie가 있으면 transcript 길이가 가장 긴 transcript

산출물:

- `work/representative_transcripts.tsv`

주요 컬럼:

```text
gene_id
gene_id_base
gene_name
transcript_id
chrom
strand
transcript_length
cds_length
utr_length
selection_reason
```

시각화:

- `fig2_region_composition_summary.png`의 일부 panel
- representative transcript의 5-prime UTR/CDS/3-prime UTR length distribution
- stacked bar: total annotated bases by region

보여주려는 것:

- 분석 단위가 gene과 transcript 사이에서 중복되지 않도록 정리되었음을 보인다.
- UTR/CDS region annotation이 충분한 gene 수를 확보했는지 확인한다.

### Step 3. 5-prime UTR / CDS / 3-prime UTR BED 생성

입력:

- `gencode.gtf`
- `representative_transcripts.tsv`

방법:

1. representative transcript의 `CDS`와 `UTR` feature만 추출한다.
2. transcript별 CDS genomic start/end를 계산한다.
3. strand-aware하게 UTR를 분류한다.

분류 규칙:

```text
+ strand:
  UTR end < CDS start   -> 5-prime UTR
  UTR start > CDS end   -> 3-prime UTR

- strand:
  UTR start > CDS end   -> 5-prime UTR
  UTR end < CDS start   -> 3-prime UTR
```

산출물:

- `work/regions.5utr.bed`
- `work/regions.cds.bed`
- `work/regions.3utr.bed`
- `work/regions.all.bed`

BED 컬럼:

```text
chrom
start
end
gene_id
score
strand
transcript_id
region
region_length
gene_name
```

시각화:

- `fig2_region_composition_summary.png`
- panel A: region별 length distribution
- panel B: gene별 region 존재 여부, 예: 5UTR/CDS/3UTR 모두 가진 gene 수

보여주려는 것:

- GTF가 `UTR`만 제공하므로, 5-prime/3-prime UTR를 직접 유도했다는 Methods 핵심을 명확히 보여준다.

### Step 4. CLIP point signal과 RNA-control region signal 계산

입력:

- `CLIP-35L33G.bam`
- `RNA-control.bam`
- `regions.*.bed`

MVP CLIP signal:

- `CLIP-35L33G.bam`에서 `NH=1` read만 사용한다.
- 각 read의 5-prime end를 strand-aware하게 1 bp BED point로 변환한다.
- point가 같은 strand의 representative transcript 5-prime UTR/CDS/3-prime UTR 중 어디에 들어가는지 count한다.

5-prime end 좌표:

```text
forward read: BED start = reference_start, BED end = reference_start + 1
reverse read: BED start = reference_end - 1, BED end = reference_end
```

주의:

- reverse read의 5-prime end는 leftmost POS가 아니라 rightmost aligned reference position이다.
- CIGAR의 deletion/skipped region을 고려해 read length가 아니라 reference-consuming CIGAR 기준 `reference_end`를 사용한다.
- `samtools view` SAM 출력과 Python 표준 파서로 구현 가능하다. `pysam`이 설치되어 있으면 `read.reference_start`, `read.reference_end`, `read.is_reverse`를 사용한다.

MVP에서 5-prime end를 쓰는 이유:

- genome-wide CIMS calling은 더 원칙적이지만 무겁다.
- 본 분석의 핵심은 amount vs position composition 모델이므로, CIMS는 시간이 남을 때 sensitivity로 수행한다.
- 보고서에서는 5-prime end point가 practical proxy이며, CIMS site가 principled point signal임을 limitation/supplementary plan에 명시한다.

RNA-control signal:

- `RNA-control.bam` read overlap을 같은 strand 기준으로 region별 count한다.
- region별 RNA-control count는 region length와 local RNA abundance를 반영하므로, CLIP point signal을 보정하는 denominator로 사용한다.

도구:

- `samtools view`: SAM/BAM read extraction
- Python 표준 파서 또는 `pysam` optional: NH=1 filtering, CIGAR-aware 5-prime end point 생성
- `bedtools intersect -s -c`: region별 CLIP point count
- `bedtools multicov -s` 또는 `bedtools coverage -s`: RNA-control stranded region count

산출물:

- `work/clip_5end_points.bed`
- `work/region_counts.tsv`

주요 컬럼:

```text
gene_id_base
region
clip_points
rna_control_reads
region_length
clip_points_cpm
rna_control_cpm
```

시각화:

- `fig2_region_composition_summary.png`
- panel C: raw CLIP point fraction by region
- panel D: RNA-control-adjusted enrichment composition by region

보여주려는 것:

- raw point fraction이 단순히 region length를 따를 수 있으므로, RNA-control 보정 composition이 왜 필요한지 보여준다.

### Step 5. RNA-control 보정 region composition 계산

입력:

- `region_counts.tsv`
- `gene_response.tsv`

region enrichment:

```text
E_5utr = (CLIP_points_5utr_cpm + pc) / (RNA_control_5utr_cpm + pc)
E_cds  = (CLIP_points_cds_cpm  + pc) / (RNA_control_cds_cpm  + pc)
E_3utr = (CLIP_points_3utr_cpm + pc) / (RNA_control_3utr_cpm + pc)
```

amount predictor:

```text
clip_points_total = clip_points_5utr + clip_points_cds + clip_points_3utr
rna_control_total = rna_control_5utr + rna_control_cds + rna_control_3utr

total_CLIP_amount_point =
log2((clip_points_total_cpm + pc) / (rna_control_total_cpm + pc))
```

`clip_points_total`은 representative transcript의 exon region(5-prime UTR/CDS/3-prime UTR)에 배정된 point pool만 사용한다. 따라서 M1의 amount와 M2의 composition은 같은 stranded point-count signal에서 나온다.

composition:

```text
f_5utr = E_5utr / (E_5utr + E_cds + E_3utr)
f_cds  = E_cds  / (E_5utr + E_cds + E_3utr)
f_3utr = E_3utr / (E_5utr + E_cds + E_3utr)
```

모델에는 세 fraction을 모두 넣지 않는다.

MVP model:

```text
M2: DeltaRD ~ covariates + total_CLIP_amount_point + f_cds + f_3utr
```

이때 5-prime UTR fraction은 implicit reference가 된다.

composition eligibility:

- total CLIP point 수가 충분한 gene만 position model에 사용한다.
- 시작 기준: total CLIP points >= 10
- 각 region의 RNA-control coverage가 최소 기준을 통과해야 한다.
- 시작 기준: RNA-control reads >= 5 in each required region
- 최종 기준은 실제 남는 n을 보고 조정하되, 본문에 n을 명시한다.

산출물:

- `work/analysis_table.tsv`
- `tables/eligibility_summary.tsv`

주요 컬럼:

```text
gene_id_base
delta_rd
total_CLIP_enrichment_featureCounts
total_CLIP_amount_point
clip_points_total
clip_points_5utr
clip_points_cds
clip_points_3utr
rna_control_5utr
rna_control_cds
rna_control_3utr
E_5utr
E_cds
E_3utr
f_5utr
f_cds
f_3utr
composition_eligible
```

시각화:

- `fig2_region_composition_summary.png`
- panel E: `f_5utr`, `f_cds`, `f_3utr` distribution
- panel F: composition-eligible gene 수와 filtering flow

보여주려는 것:

- 위치 변수는 raw count fraction이 아니라 RNA-control로 보정한 relative preference임을 보여준다.
- M1의 amount predictor와 M2의 composition predictor가 같은 stranded point-count signal에서 나온다는 것을 보장한다.
- 분석 대상 gene 수가 충분한지 투명하게 제시한다.

### Step 6. nested model로 position effect 검정

입력:

- `analysis_table.tsv`

model fitting set:

- M0/M1/M2는 모두 동일한 composition-eligible gene set에서 적합한다.
- 즉, `DeltaRD`, covariates, `total_CLIP_amount_point`, `f_cds`, `f_3utr`가 모두 정의된 gene만 사용한다.
- Step 1의 baseline 재현 figure는 더 큰 expression/count-filtered gene set에서 따로 보고할 수 있지만, nested model의 R2와 F-test에는 섞지 않는다.
- 모델 실행 전 `n_M0 == n_M1 == n_M2`를 assert하고 `model_summary.tsv`에 n을 기록한다.

covariates:

- `rna_expression`: RNA-control 또는 평균 RNA CPM
- `gene_length` 또는 representative transcript length

model:

```text
M0: DeltaRD ~ log_rna_expression + log_gene_length
M1: DeltaRD ~ log_rna_expression + log_gene_length + total_CLIP_amount_point
M2: DeltaRD ~ log_rna_expression + log_gene_length + total_CLIP_amount_point + f_cds + f_3utr
```

통계:

- OLS regression
- `R2(M0)`, `R2(M1)`, `R2(M2)`
- `Delta R2_amount = R2(M1) - R2(M0)`
- `Delta R2_position = R2(M2) - R2(M1)`
- nested F-test: primary test for M1 vs M2

도구:

- Python `statsmodels.OLS` + `anova_lm`을 사용할 수 있으면 가장 간단하다.
- 현재 환경에서 `statsmodels`가 없으면 Python `numpy`로 OLS를 적합하고, F-statistic은 SSE와 자유도로 직접 계산한다.
- p-value는 `scipy.stats.f.sf`로 계산한다.

산출물:

- `tables/model_summary.tsv`

시각화:

- `fig3_nested_model.png`
- panel A: barplot of R2 for M0/M1/M2
- panel B: barplot of Delta R2 amount vs Delta R2 position
- panel C: M2 standardized beta for total amount and composition terms

보여주려는 것:

- total CLIP amount가 translation derepression을 설명하는 baseline인지 확인한다.
- 그 위에 region composition이 설명력을 실제로 추가하는지 검정한다.
- 주장의 핵심은 계수 하나가 아니라 `Delta R2_position`이다.

### Step 7. binding dose 분석

입력:

- `analysis_table.tsv`

dose 변수:

- MVP: `clip_points_total`
- 보조: `clip_points_cds`, `clip_points_3utr`
- optional: CIMS/high-confidence site count

분석:

1. total CLIP point count로 gene을 분위수 그룹으로 나눈다.
   - 예: Q1, Q2, Q3, Q4 또는 bottom 50%, 50-80%, top 20%, top 5%
2. 각 그룹의 DeltaRD 분포를 비교한다.
3. Spearman correlation과 Mann-Whitney U test를 보조로 계산한다.

산출물:

- `fig4_dose_response.png`

시각화:

- violin/box plot: dose group별 DeltaRD
- line plot: dose group별 median DeltaRD
- annotation: group n, Spearman rho

보여주려는 것:

- region composition 결과가 약하더라도, LIN28A binding burden 자체가 translation derepression과 연결되는지 확인한다.
- 원 논문에서 언급한 multiple binding sites/mRNA 모델과 연결한다.

## 6. Optional / supplementary analyses

MVP 이후 시간이 남으면 아래 순서로 추가한다.

### Optional 0. stranded gene/CDS count 재계산

목적:

- `read-counts.txt`가 unstranded featureCounts 결과이므로, model response인 DeltaRD도 stranded count 또는 CDS-limited count로 재계산했을 때 결과가 유지되는지 확인한다.

방법:

- `featureCounts -s 1`로 RNA/RPF BAM을 gene 또는 representative CDS에 다시 count한다.
- Actb 등 known gene으로 `-s 1`과 `-s 2` 중 어느 방향이 sense인지 먼저 검증한다.
- 기존 `read-counts.txt` 기반 DeltaRD와 stranded/CDS-limited DeltaRD의 correlation을 보고한다.

### Optional 1. raw point composition vs adjusted composition 비교

목적:

- raw point composition 효과가 region length/RNA coverage confound였는지 확인한다.

시각화:

- raw `f_region`과 adjusted `f_region` scatter
- raw composition model과 adjusted composition model의 Delta R2 비교

### Optional 2. permutation robustness check

목적:

- `Delta R2_position`의 robustness를 검증한다.

주의:

- 단순 gene-level shuffle은 composition-covariate 관계를 깨므로 과낙관적 null을 만들 수 있다.
- 가능하면 expression/length bin 안에서 permutation하거나, composition predictor를 covariate에 대해 residualize한 뒤 residual을 섞는다.

시각화:

- null distribution of Delta R2_position
- observed Delta R2_position vertical line

### Optional 3. CIMS-based point signal

목적:

- HITS-CLIP에서 원칙적으로 더 적합한 crosslink point인 CIMS site를 사용해 MVP 결과가 유지되는지 확인한다.

방법:

- week3 pileup/entropy logic을 expressed gene 또는 strong binder subset에 확장한다.
- high-confidence entropy site를 region에 배정한다.
- CIMS site count 기반 dose/composition 분석을 MVP와 비교한다.

### Optional 4. CDS-restricted Delta ribosome density

목적:

- 원 논문 정의에 더 가까운 DeltaRD를 만든다.

방법:

- representative CDS BED에 대해 RPF/RNA reads를 재카운트한다.
- gene-level featureCounts DeltaRD와 CDS-limited DeltaRD의 correlation을 확인한다.
- 주요 model을 CDS-limited DeltaRD로 재실행한다.

## 7. 최종 보고서에서 각 그림이 담당할 주장

| Figure | 내용 | 보여주려는 주장 |
|---|---|---|
| Fig 1 | total CLIP enrichment vs DeltaRD | 원 논문/W1 baseline이 재현된다 |
| Fig 2 | representative transcript region과 adjusted composition 요약 | 5-prime UTR/CDS/3-prime UTR 위치 변수를 신뢰성 있게 만들었다 |
| Fig 3 | M0/M1/M2 nested model Delta R2 및 standardized beta | total amount 대비 position composition의 추가 설명력을 검정했다 |
| Fig 4 | binding dose group별 DeltaRD | LIN28A 결합량/dose가 번역 탈억제와 약하게 연결되는지 확인했다 |

## 8. 예상 결론 시나리오

### Scenario A. `Delta R2_position`이 유의하고 CDS fraction이 중요

해석:

- 같은 total CLIP amount라도 CDS 쪽에 상대적으로 결합이 농축된 gene에서 Lin28a knockdown 후 번역 탈억제가 더 크다.
- LIN28A가 elongation 또는 CDS 상 ribosome traffic에 영향을 줄 가능성을 제시한다.

주의:

- 상관 기반 분석이므로 causality로 과장하지 않는다.
- DeltaRD 자체가 CDS translation을 반영하므로 해석에 한계를 둔다.

### Scenario B. `Delta R2_position`은 작고 total amount/dose만 중요

해석:

- LIN28A 억제 효과는 특정 region 위치보다 전체 결합량 또는 binding burden에 더 의존할 수 있다.
- 원 논문의 transcript-level model을 지지한다.

### Scenario C. raw composition만 효과가 있고 adjusted composition에서는 사라짐

해석:

- 위치 효과로 보였던 신호는 region length 또는 RNA coverage 차이에서 온 confound일 가능성이 크다.
- 최종 결론은 adjusted composition 결과를 기준으로 낸다.

### Scenario D. total CLIP baseline도 재현되지 않음

해석:

- biological conclusion을 내기 전에 filtering, normalization, gene_id merge, count 계산을 점검한다.
- W1/Fig 4D 재현 실패는 pipeline 문제일 가능성이 높다.

## 9. 최종 MVP 체크리스트

- [ ] `gene_response.tsv` 생성
- [ ] `read-counts.txt`와 GTF-derived table의 `gene_id` version suffix 제거 및 merge 확인
- [ ] total CLIP enrichment vs DeltaRD baseline plot 생성
- [ ] representative transcript set 생성
- [ ] 5-prime UTR/CDS/3-prime UTR BED 생성
- [ ] known gene으로 strandedness 방향 확인 (`-s 1` vs `-s 2`)
- [ ] NH=1 CLIP 5-prime-end point BED 생성, reverse read는 rightmost coordinate 사용
- [ ] strand-aware region별 CLIP point count와 RNA-control count 생성
- [ ] adjusted composition 계산
- [ ] composition-eligible gene 수 보고
- [ ] `total_CLIP_amount_point = log2((clip_points_total_cpm + pc) / (rna_control_total_cpm + pc))` 계산
- [ ] M0/M1/M2가 동일한 composition-eligible gene set에서 적합되는지 확인 (`n_M0 == n_M1 == n_M2`)
- [ ] M0/M1/M2 nested model 실행
- [ ] `Delta R2_position` nested F-test 계산
- [ ] dose-response figure 생성
- [ ] 모든 figure에 n과 필터 기준 표시
