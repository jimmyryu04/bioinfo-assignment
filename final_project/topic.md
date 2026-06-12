# Your Own Analysis 최종 주제

## 결론

최종 주제는 다음으로 정한다.

> **LIN28A 결합의 총량과 위치 구성이 Lin28a knockdown 후 번역 탈억제를 얼마나 설명하는가**

영문 제목 후보:

> **Amount- and position-aware analysis of LIN28A binding as a predictor of translational derepression after Lin28a knockdown**

핵심은 다음 한 문장이다.

> 원 논문은 LIN28A의 transcript-level 결합량이 ribosome density 변화와 연결된다는 것을 보였지만, 본 분석은 LIN28A 결합을 **총량(amount)**과 **RNA abundance/region length로 보정한 위치 구성(enrichment-adjusted composition)**으로 분리하여, 같은 결합량에서도 어느 region에 상대적으로 농축되어 결합하는지가 번역 탈억제 설명력을 추가하는지 검정한다.

이 문서는 기존 `topic.md`의 장점인 단일 명제, 예측/회귀 프레이밍, 사전 가설을 유지하면서, 추가 피드백의 핵심 4가지를 반영한 최종 버전이다.

- 위치 효과를 총 결합량에서 분리하기 위해 **region별 raw enrichment**가 아니라 **region composition**을 주 predictor로 사용한다.
- 조건당 생물학적 반복이 없는 datapack의 한계를 명시하고, 분석을 **gene 간 association**으로 제한한다.
- region composition은 raw CLIP point fraction이 아니라 **region별 CLIP/RNA-control enrichment를 정규화한 composition**으로 정의해 gene 구조와 region length confound를 줄인다.
- 결합 위치는 CLIP read 전체 overlap이 아니라 **CIMS site** 같은 point signal로 세는 것을 원칙으로 한다. 일정상 genome-wide CIMS calling이 부담되면 5-prime end 1 bp를 primary로 쓰고 CIMS를 sensitivity로 둔다.
- gene-level response와 transcript-level region annotation을 맞추기 위해 **gene당 대표 transcript 1개**를 사용하고, multimapper 처리를 명시한다.
- headline endpoint는 `M1(amount)` 대비 `M2(amount + position)`의 `Delta R2_position`이며, nested F-test 또는 permutation test로 검정한다.

---

## 1. 최종 연구 질문

### Main question

LIN28A 결합의 **총량**이 번역 탈억제를 설명하는가, 아니면 총 결합량을 통제한 뒤에도 **5-prime UTR, CDS, 3-prime UTR 중 어느 region에 결합이 상대적으로 농축되는지**가 Lin28a knockdown 후 ribosome density 증가를 추가로 설명하는가?

### 세부 질문

1. **Amount question**
   total CLIP enrichment 또는 total CLIP signal이 Delta ribosome density를 재현성 있게 예측하는가?

2. **Position composition question**
   total CLIP signal을 통제한 상태에서, RNA-control로 보정한 LIN28A 결합 농축도의 region별 구성비가 추가 설명력을 가지는가?

   ```text
   E_5utr = (CLIP_signal_5utr + pc) / (RNA_control_signal_5utr + pc)
   E_cds  = (CLIP_signal_cds  + pc) / (RNA_control_signal_cds  + pc)
   E_3utr = (CLIP_signal_3utr + pc) / (RNA_control_signal_3utr + pc)

   f_5utr = E_5utr / (E_5utr + E_cds + E_3utr)
   f_cds  = E_cds  / (E_5utr + E_cds + E_3utr)
   f_3utr = E_3utr / (E_5utr + E_cds + E_3utr)
   ```

   단, `f_5utr + f_cds + f_3utr = 1`이므로 세 변수를 모두 회귀식에 넣지 않는다. 두 fraction만 넣거나, additive log-ratio/alr 또는 isometric log-ratio/ilr 변환을 사용한다.

3. **Dose question**
   mRNA당 LIN28A 결합량 또는 high-confidence binding site 수가 많을수록 Lin28a knockdown 후 번역 탈억제가 더 커지는가?

4. **Region-independent question**
   위치 구성비의 추가 설명력이 작고 total CLIP amount 또는 dose가 대부분의 신호를 설명하는가?

이 설계에서 진짜 새 질문은 `M2 - M1`, 즉 **총 결합량을 통제한 뒤 RNA-control/길이 효과를 보정한 위치 구성비가 설명력을 추가하는가**이다.

---

## 2. 왜 7-10쪽 보고서가 가능한가?

**가능하다.** 단순 figure 재현이 아니라, 원 논문의 두 결과를 연결해 새 가설을 검정하는 분석이기 때문이다.

보고서 구성은 다음처럼 잡는다.

| 섹션 | 내용 | 예상 분량 |
|---|---|---|
| Abstract | LIN28A binding amount/composition과 translation derepression 관계 요약 | 0.5쪽 |
| Introduction | LIN28A, CLIP-seq, ribosome profiling, 원 논문 결론 2/3, 남은 mechanistic question | 1.0-1.5쪽 |
| Research question | transcript-level total binding 분석의 한계와 amount/position 분리 필요성 | 0.5쪽 |
| Methods | 대표 transcript, UTR 유도, CIMS/point-based CLIP counting, RNA-control 보정 composition predictor, nested/permutation model | 2.0쪽 |
| Results 1 | total CLIP amount와 Delta ribosome density의 baseline association | 0.8쪽 |
| Results 2 | region composition과 predictor correlation 구조 | 1.0쪽 |
| Results 3 | M0/M1/M2 nested model Delta R2, partial correlation | 1.2쪽 |
| Results 4 | binding dose/site count와 Delta ribosome density dose-response | 1.0쪽 |
| Discussion | 원 논문과의 연결, null result 해석, 반복 부재와 coordinate 한계 | 1.0-1.5쪽 |

본문 figure는 4-5개로 압축한다.

1. 분석 개요도: amount predictor, enrichment-adjusted position composition predictor, dose predictor
2. 대표 transcript 기반 5-prime UTR/CDS/3-prime UTR annotation 및 CLIP point signal 분포
3. total CLIP amount vs Delta ribosome density baseline plot
4. region composition + nested model Delta R2/permutation summary
5. binding dose/site count에 따른 Delta ribosome density dose-response

### 실행 범위: MVP와 optional

보고서 본문은 아래 MVP 분석만으로 완성되도록 한다.

| 구분 | 분석 | 보고서 내 역할 |
|---|---|---|
| MVP 1 | `read-counts.txt` 기반 Delta ribosome density와 total CLIP enrichment 계산 | Fig 4D baseline 재현/검증 |
| MVP 2 | representative transcript set 구축 및 5-prime UTR/CDS/3-prime UTR BED 생성 | Methods 핵심 |
| MVP 3 | NH=1 CLIP 5-prime-end point와 RNA-control region count | primary region signal |
| MVP 4 | RNA-control 보정 composition 계산 후 M0/M1/M2 nested model | 핵심 novelty |
| MVP 5 | `Delta R2_position` nested F-test | primary test |
| MVP 6 | total point dose 분위수별 Delta ribosome density | robustness/result 확장 |

아래 항목은 시간이 남을 때 supplementary 또는 sensitivity analysis로 수행한다.

- predictor correlation matrix/VIF
- genome-wide CIMS 기반 counting
- pseudocount sensitivity
- low-count filter sensitivity
- permutation test
- raw point composition vs RNA-control-adjusted composition 비교
- membrane gene subset 분석
- full-read overlap 방식과 point-based counting 비교
- CDS-restricted Delta ribosome density 재카운트
- alr/ilr composition 변환

---

## 3. 기존 논문과 무엇이 다른가?

### 원 논문에서 한 것

Cho et al. (2012)은 다음을 보였다.

- LIN28A CLIP tag는 mRNA의 여러 region에 분포한다.
  - Fig 3C: mRNA region별 CLIP tag density
- LIN28A의 transcript-level CLIP enrichment는 Lin28a knockdown 후 ribosome density 증가와 연결된다.
  - Fig 4D/4E: 전체 transcript 단위의 CLIP enrichment vs ribosome density change
- LIN28A 표적은 ER-associated 또는 membrane protein mRNA에 풍부하다.
  - Fig 5B, Fig S6A
- 논문 Discussion에서는 LIN28A 결합 site가 CDS와 3-prime UTR 전반에 여러 개 존재하므로, translation initiation 이전 단계 또는 elongation 중간 단계에 영향을 줄 수 있다고 제안했지만, 그 연결을 정량적으로 풀지는 않았다.

### 이번 분석에서 새로 하는 것

본 분석은 원 논문이 따로 제시한 두 축을 직접 연결한다.

| 항목 | 원 논문 | 본 분석 |
|---|---|---|
| 결합량 단위 | transcript/gene 전체 CLIP enrichment | total amount와 RNA-control 보정 region composition을 분리 |
| 위치 변수 | region별 CLIP density를 서술 | 같은 total amount에서 5-prime UTR/CDS/3-prime UTR 상대 농축도가 다른지 검정 |
| 결합 구조 | 평균적으로 많은 sites/mRNA가 존재한다고 언급 | mRNA별/region별 결합량 또는 site 수를 dose 변수로 사용 |
| 주요 질문 | LIN28A가 결합한 mRNA는 KD 후 번역이 증가하는가? | 총량을 통제한 뒤 결합 위치 구성이 번역 탈억제를 추가 설명하는가? |
| 통계 설계 | transcript-level correlation / distribution 비교 | nested model Delta R2, partial correlation, composition-aware predictor |
| 해석 | LIN28A는 translation suppressor | 번역 억제의 amount 의존성, position 의존성, dose 의존성 평가 |

따라서 이 주제는 Fig 3C나 Fig 4D를 그대로 재현하는 것이 아니다. Fig 3C의 "where LIN28A binds"와 Fig 4D의 "translation derepression after knockdown"을 결합하고, 특히 **total binding amount와 region length/expression을 통제한 위치 효과**를 새로 묻는다.

---

## 4. 기존 논문을 어떻게 확장하는가?

원 논문의 중심 결론:

> LIN28A binds many spliced mRNAs and suppresses translation, especially for ER-associated targets.

본 분석의 확장:

> LIN28A-mediated translational suppression may depend on both the amount and RNA-normalized positional composition of LIN28A binding within each mRNA.

확장 논리는 다음과 같다.

1. 원 논문은 LIN28A가 mRNA에 광범위하게 결합한다고 보였다.
2. 원 논문은 LIN28A 결합이 mRNA abundance보다 ribosome density 변화와 더 잘 연결된다고 보였다.
3. 그러나 transcript 전체 CLIP enrichment만으로는 LIN28A가 translation의 어느 단계에 영향을 주는지 구체적으로 알기 어렵다.
4. 본 분석은 `amount predictor`와 `RNA-control-adjusted position composition predictor`를 분리하여, **총량이 많은 gene이 더 억제되는지**와 **같은 총량에서도 결합이 특정 region에 상대적으로 농축되면 더 억제되는지**를 구분한다.
5. CDS 구성비가 추가 설명력을 가지면 elongation 또는 CDS 상 ribosome traffic과 연결해 해석할 수 있다.
6. 5-prime UTR 구성비가 중요하면 initiation 조절 가능성을 논의할 수 있다.
7. 3-prime UTR 구성비가 중요하면 mRNA localization, closed-loop translation regulation, regulatory complex recruitment 가능성을 논의할 수 있다.
8. 위치 구성비의 추가 설명력이 작고 amount/dose가 대부분을 설명하면, LIN28A 억제 효과는 특정 위치보다 다중 결합량, 반복 결합, 또는 peri-ER localization에 의해 좌우될 가능성이 있다.

즉, 이 주제는 원 논문의 결론을 반복하지 않고, **LIN28A binding architecture가 translation repression mechanism에 주는 단서**를 찾는 확장 분석이다.

---

## 5. 사전 가설

### Main hypothesis

LIN28A의 total binding amount가 translation derepression을 설명하는 baseline predictor이고, RNA-control로 보정한 region composition 또는 binding dose가 그 위에 추가 설명력을 제공할 수 있다.

### Alternative hypotheses

| 가설 | 예상 결과 | 해석 |
|---|---|---|
| Amount model | total CLIP amount가 Delta ribosome density와 연결되고 adjusted region composition Delta R2는 작음 | 특정 위치보다 전체 결합량 또는 ER localization 중심 모델 |
| CDS-composition model | total amount 통제 후 CDS adjusted fraction 또는 CDS-vs-reference log-ratio가 추가 설명력 제공 | elongation 또는 CDS 상 ribosome movement 방해 가능성 |
| UTR-composition model | 5-prime UTR 또는 3-prime UTR adjusted fraction/log-ratio가 추가 설명력 제공 | initiation, mRNA localization, regulatory complex 가능성 |
| Dose model | binding count/site count가 증가할수록 derepression 증가 | 다중 결합 또는 binding burden이 중요 |
| Pipeline-check null | total CLIP도 Delta ribosome density를 거의 예측하지 못함 | W1/Fig 4D 재현과 맞지 않으므로 pipeline/filtering 점검 필요 |
| Length-confounded model | raw point composition에서는 효과가 보이나 RNA-control/length 보정 후 사라짐 | 위치 효과가 아니라 gene 구조 또는 region 길이 차이가 만든 confound |
| True null | total amount는 재현되지만 adjusted composition/dose의 추가 설명력이 작음 | 위치 특이 모델보다 transcript-level 또는 다른 생물학적 요인이 중요 |

negative result도 결론이 되도록 설계한다. 특히 **total amount는 재현되지만 adjusted composition의 추가 설명력이 작다**는 결과는 실패가 아니라, LIN28A 번역 억제가 특정 mRNA region 하나로 설명되지 않는다는 의미 있는 결론이다.

---

## 6. 핵심 리스크와 대응 설계

### 6.1 region length와 expression confound

raw CLIP point fraction은 LIN28A의 위치 선호만 반영하지 않는다. 3-prime UTR가 긴 gene은 LIN28A가 균일하게 결합해도 3-prime UTR point가 많아질 수 있고, region별 RNA abundance/coverage 차이도 composition에 섞인다.

대응:

- primary composition은 raw CLIP point fraction이 아니라 **region별 CLIP/RNA-control enrichment**에서 만든다.
- RNA-control region signal은 region length와 local expression/coverage를 함께 반영하므로, CLIP signal을 RNA-control로 나누어 gene 구조 차이를 줄인다.
- raw point composition은 보조 분석으로만 사용하고, 보정 전후 결과가 바뀌면 length/expression confound로 해석한다.
- region length는 covariate 또는 sensitivity analysis 변수로 기록한다.

### 6.2 위치 효과와 총량 효과의 분리

region별 raw CLIP enrichment는 총량과 위치 정보를 동시에 담는다. 예를 들어 CDS enrichment가 큰 gene은 단순히 전체 CLIP amount가 큰 gene일 수 있다.

대응:

- primary position predictor는 raw point composition이 아니라 **RNA-control-adjusted region composition**으로 정의한다.
- amount predictor는 total CLIP enrichment 또는 total point signal로 따로 둔다.
- composition 변수는 합이 1이므로 세 fraction을 모두 회귀식에 넣지 않는다.
- 기본 구현은 두 fraction만 투입하거나 alr/ilr 변환을 사용한다.
- 핵심 비교는 `M1 = amount`와 `M2 = amount + composition`의 Delta R2이다.

### 6.3 다중공선성

LIN28A는 mRNA 전체에 광범위하게 결합하므로 region별 signal은 서로 강하게 상관될 수 있다. 단순 multiple regression의 계수 부호나 크기를 그대로 해석하면 위험하다.

대응:

- predictor 간 Pearson/Spearman correlation matrix를 먼저 제시한다.
- VIF 또는 유사한 공선성 지표를 계산한다.
- raw regression coefficient 중심 해석을 피한다.
- nested model Delta R2, partial correlation, variance partitioning 중심으로 해석한다.
- "어느 region의 계수가 유의하다"보다 "composition 정보가 total amount 대비 추가 설명력을 주는가"를 주요 결론으로 삼는다.

### 6.4 생물학적 반복 부재

현재 로컬 datapack은 조건당 라이브러리 1개만 포함한다. 예를 들어 RPF-siLuc, RPF-siLin28a, RNA-siLuc, RNA-siLin28a, CLIP-35L33G가 각각 1개씩이다. 따라서 condition 내 분산을 추정할 수 없고, gene별 differential translation p-value를 계산할 수 없다.

대응:

- 본 분석은 gene별 유의성 검정이 아니라 **gene을 데이터포인트로 한 across-gene association 분석**임을 Methods에 명시한다.
- p-value가 있다면 biological replicate 기반 differential test가 아니라 regression/correlation association의 p-value로 해석한다.
- 결론은 "개별 gene이 유의하게 변했다"가 아니라 "유전자들 사이에서 binding architecture와 Delta ribosome density가 함께 변하는 경향"으로 제한한다.
- 원 논문의 반복 실험과는 데이터 수준이 다르므로, 보고서의 Limitation에 명시한다.

### 6.5 5-prime UTR signal 희박성

논문 Fig 3C 자체가 5-prime UTR에서 LIN28A 결합이 고갈됨을 보였다. 5-prime UTR는 길이도 짧기 때문에 read count가 적고 composition 추정이 noisy할 수 있다.

대응:

- region length와 RNA-control signal을 함께 기록한다.
- low-count region은 필터링하거나 pseudocount sensitivity analysis를 수행한다.
- 5-prime UTR에 CLIP signal이 있는 gene 수를 별도로 보고한다.
- composition 분석에 들어갈 만큼 total CLIP signal이 충분한 gene 수를 별도로 보고한다.
- 5-prime UTR negative result는 "효과 없음"이 아니라 "현재 데이터에서 독립적 기여를 검출하기 어려움"으로 해석한다.

### 6.6 UTR annotation 문제

현재 GENCODE GTF는 `UTR` feature만 제공하고 `five_prime_utr`, `three_prime_utr`로 나뉘어 있지 않다. 따라서 UTR 구분은 Methods의 핵심 작업이다.

대응:

- gene당 대표 protein-coding transcript 1개를 먼저 정한다.
- 대표 transcript의 CDS 시작/끝 좌표를 구한다.
- `+` strand에서는 CDS 시작보다 upstream인 UTR을 5-prime UTR, CDS 끝보다 downstream인 UTR을 3-prime UTR로 분류한다.
- `-` strand에서는 방향을 반대로 적용한다.
- CDS와 겹치는 UTR 또는 ambiguous transcript는 제외하거나 별도 기록한다.
- 5-prime UTR, CDS, 3-prime UTR 중 하나가 없는 gene의 처리 기준을 명시한다.

### 6.7 gene/transcript 단위 정합

Delta ribosome density는 `read-counts.txt`에서 gene-level로 계산되지만, UTR/CDS region은 transcript-level annotation에서 나온다. isoform을 모두 사용하면 한 gene의 여러 transcript region이 gene-level response 하나에 중복 연결될 수 있다.

대응:

- gene당 대표 transcript 1개를 사용한다.
- 우선순위는 protein-coding, transcript_support_level 1, basic/appris/principal tag, longest CDS 또는 longest transcript 순으로 정한다.
- 최종 merge key는 version을 정리한 gene_id로 통일한다.
- 유효 gene 수를 모든 주요 단계에서 보고한다.

### 6.8 CLIP 결합 위치 counting 방식

CLIP read 전체를 region에 overlap시키면 read가 region boundary를 가로질러 신호가 번질 수 있다. 특히 짧은 5-prime UTR에서는 치명적일 수 있다.

대응:

- primary region signal은 full-read overlap보다 **point-based signal**을 사용한다.
- 이 논문 데이터에서는 crosslink marker가 5-prime truncation이 아니라 **CIMS substitution**이므로, 원칙적으로 CIMS site를 primary point signal로 사용한다.
- genome-wide CIMS calling이 일정상 부담되면 CLIP read의 5-prime end 1 bp를 practical primary로 사용하고, CIMS site count를 targeted sensitivity analysis로 수행한다.
- full-read overlap count는 sensitivity analysis로만 사용한다.

### 6.9 composition-eligible gene 필터

total CLIP point가 너무 적은 gene은 region composition이 이산적으로 튀며 의미가 약하다. 예를 들어 total point가 3개인 gene의 fraction은 0, 0.33, 0.67, 1.0 같은 coarse value만 가진다.

대응:

- position composition 분석에는 별도의 minimum total CLIP signal threshold를 둔다.
- 예: total CLIP points >= 10 또는 high-confidence CIMS sites >= 3 같은 기준을 sensitivity analysis로 비교한다.
- 각 region의 RNA-control coverage도 최소 기준을 통과해야 한다. RNA-control signal이 거의 0인 region은 `CLIP/RNA` enrichment가 pseudocount에 의해 허위 중립값처럼 보일 수 있기 때문이다.
- expression/RPF 필터를 통과한 gene 수와 composition-eligible gene 수를 따로 보고한다.
- 5-prime UTR signal 보유 gene 수, CDS signal 보유 gene 수, 3-prime UTR signal 보유 gene 수도 함께 보고한다.

### 6.10 multimapper 처리

CLIP BAM에는 multimapping read가 포함될 수 있다. 논문은 single-best-hit 중심으로 분석했으므로, region counting에서도 처리 방침이 필요하다.

대응:

- primary 분석은 NH=1 read만 사용한다.
- NH tag가 없거나 필터가 너무 강할 경우 MAPQ 또는 unique alignment 기준을 대안으로 쓴다.
- 필요하면 fractional counting을 sensitivity analysis로 둔다.

### 6.11 종속변수의 CDS 중심성

Delta ribosome density는 RPF/RNA 비율이며, 보통 CDS에 매핑된 ribosome footprint로 계산된다. 따라서 종속변수 자체가 CDS translation을 반영한다.

대응:

- CDS composition이 강하게 보이더라도 causality로 과장하지 않는다.
- CDS signal은 elongation과 연결될 수 있지만 상관 기반 분석임을 명시한다.
- 3-prime UTR 결합은 CDS 번역과 간접적으로 연결될 수 있으므로 mRNA localization 또는 translation regulatory complex 관점에서 신중히 해석한다.

---

## 7. 분석 설계

### 입력 데이터

로컬 데이터만 사용한다.

- `week1/work/binfo1-datapack1/CLIP-35L33G.bam`
- `week1/work/binfo1-datapack1/RNA-control.bam`
- `week1/work/binfo1-datapack1/RNA-siLuc.bam`
- `week1/work/binfo1-datapack1/RNA-siLin28a.bam`
- `week1/work/binfo1-datapack1/RPF-siLuc.bam`
- `week1/work/binfo1-datapack1/RPF-siLin28a.bam`
- `week1/work/binfo1-datapack1/gencode.gtf`
- `week1/work/read-counts.txt`

추가 genome FASTA 다운로드 없이 진행하는 것을 기본으로 한다.

### Step 1. gene-level response 계산

week1 방식으로 gene별 Delta ribosome density를 계산한다.

```text
Delta ribosome density =
log2((RPF_siLin28a / RNA_siLin28a) / (RPF_siLuc / RNA_siLuc))
```

동시에 baseline predictor인 total CLIP enrichment도 계산한다.

```text
total CLIP enrichment =
log2(CLIP_35L33G / RNA_control)
```

이 baseline은 W1/Fig 4D 재현에 해당한다. 분석의 새 기여는 total CLIP 위에 position composition과 dose가 설명력을 추가하는지를 보는 것이다.

주의:

- 원 논문의 ribosome density는 RPF와 RNA를 CDS에 매핑해 계산한 값이다.
- `week1/work/read-counts.txt`는 gene/exon-level featureCounts 결과이므로, 이를 사용하면 원 논문 정의의 근사값이다.
- 시간이 허용되면 representative CDS에 대해 RPF/RNA를 재카운트하여 CDS-limited Delta ribosome density를 만들고, 그렇지 않으면 gene-level 근사임을 Methods에 명시한다.

### Step 2. 대표 transcript set 구축

GTF에서 gene당 대표 protein-coding transcript 1개를 선택한다.

우선순위:

1. protein-coding transcript
2. transcript_support_level 1
3. `basic`, `appris_principal`, `CCDS` 등 신뢰도 tag
4. longest CDS 또는 longest transcript

출력:

- `representative_transcripts.tsv`
- 대표 transcript를 가진 gene 수
- 최종 분석에 들어간 gene 수

### Step 3. 5-prime UTR/CDS/3-prime UTR annotation 생성

대표 transcript의 `CDS`와 `UTR`를 추출한다. UTR는 strand-aware하게 5-prime UTR과 3-prime UTR로 직접 분류한다.

생성할 BED:

- `regions.5utr.bed`
- `regions.cds.bed`
- `regions.3utr.bed`
- `regions.exon.bed`

각 BED는 gene_id, transcript_id, region, length 정보를 포함한다.

### Step 4. point-based CLIP signal 생성

CLIP BAM에서 primary signal을 만든다.

기본 방침:

- NH=1 read만 사용한다.
- 원칙적으로 week3 entropy/CIMS logic을 genome-wide 또는 expressed gene subset에 적용해 high-confidence CIMS site를 호출한다.
- CIMS site point가 어느 representative transcript region에 들어가는지 count한다.

일정상 CIMS calling이 부담될 경우:

- CLIP read의 5-prime end를 1 bp BED로 변환해 practical primary signal로 사용한다.
- CIMS site count는 대표 target 또는 strong binder subset에서 sensitivity analysis로 수행한다.

full-read overlap count는 보조 sensitivity analysis로만 사용한다.

### Step 5. amount와 RNA-control-adjusted composition predictor 계산

각 gene에 대해 다음을 계산한다.

Amount predictor:

```text
total_CLIP_amount =
log2((CLIP_points_total + pseudocount) / (RNA_control_signal + pseudocount))
```

또는 week1 방식의 total CLIP enrichment를 baseline으로 사용한다.

Region-level enrichment:

```text
E_5utr = (CLIP_points_5utr + pc) / (RNA_control_signal_5utr + pc)
E_cds  = (CLIP_points_cds  + pc) / (RNA_control_signal_cds  + pc)
E_3utr = (CLIP_points_3utr + pc) / (RNA_control_signal_3utr + pc)
```

Position composition predictor:

```text
f_5utr = E_5utr / (E_5utr + E_cds + E_3utr)
f_cds  = E_cds  / (E_5utr + E_cds + E_3utr)
f_3utr = E_3utr / (E_5utr + E_cds + E_3utr)
```

회귀에는 세 fraction을 그대로 모두 넣지 않는다.

대안 1:

```text
M2: DeltaRD ~ covariates + total_CLIP_amount + f_cds + f_3utr
```

이때 5-prime UTR fraction은 reference 역할을 한다.

대안 2:

```text
alr_cds  = log(f_cds / f_5utr)
alr_3utr = log(f_3utr / f_5utr)
```

5-prime UTR가 매우 sparse하면 reference를 CDS로 바꾸거나 ilr 변환을 사용한다.

Composition eligibility:

- total CLIP point 또는 CIMS site 수가 너무 적은 gene은 position composition 분석에서 제외한다.
- threshold는 main analysis와 sensitivity analysis로 나누어 보고한다.
- excluded gene은 amount/dose 분석에는 포함할 수 있지만, position composition model에는 넣지 않는다.

### Step 6. 공선성 및 유효 n 진단

다음을 보고한다.

- total_CLIP_amount와 composition predictor 간 correlation
- VIF 또는 condition number
- final model에 들어간 gene 수
- composition-eligible gene 수
- 5-prime UTR, CDS, 3-prime UTR에 CLIP signal이 있는 gene 수
- low-count filter 이후 남은 gene 수

### Step 7. nested model로 position effect 검정

primary model:

```text
M0: DeltaRD ~ expression + gene_length
M1: DeltaRD ~ expression + gene_length + total_CLIP_amount
M2: DeltaRD ~ expression + gene_length + total_CLIP_amount + position_composition
```

핵심 비교:

```text
Delta R2_amount   = R2(M1) - R2(M0)
Delta R2_position = R2(M2) - R2(M1)
```

Primary endpoint:

```text
Observed Delta R2_position
```

검정:

- primary test는 nested model F-test로 `M1`과 `M2`를 비교한다.
- permutation test는 robustness check로 둔다. 단순히 composition을 gene 간에 섞으면 composition과 length/expression 같은 covariate의 관계까지 깨져 null이 과낙관적일 수 있다.
- permutation을 수행할 경우, 가능하면 composition predictor를 covariate에 대해 residualize한 뒤 residual을 섞거나, 비슷한 expression/length bin 안에서 permutation한다.
- 관측된 `Delta R2_position`이 permutation null의 상위 몇 %에 해당하는지 empirical rank로 보조 보고한다.

해석:

- `Delta R2_amount`는 원 논문/W1의 transcript-level binding effect를 재현하는 baseline이다.
- `Delta R2_position`이 본 분석의 핵심 novelty이다.
- `Delta R2_position`이 작으면 region-independent 또는 amount-driven model을 지지한다.

### Step 8. binding dose 분석

position composition이 noisy하거나 추가 설명력이 작을 때를 대비해, binding dose 분석을 보조 축으로 수행한다.

가능한 dose 변수:

- total CLIP point count
- normalized CLIP point density
- region별 point count
- CIMS/high-confidence site count

기본 분석:

- total binding dose 분위수별 Delta ribosome density 분포
- region별 binding count 분위수별 Delta ribosome density 분포
- Spearman correlation 또는 Mann-Whitney U test

이 분석은 "어디에 결합하는가"뿐 아니라 "얼마나 많이 결합하는가"가 번역 탈억제와 연결되는지를 확인한다.

---

## 8. 예상 결과 해석 프레임

### 결과 A: amount가 강하고 composition Delta R2가 작음

해석:

- LIN28A 억제 효과는 특정 region보다 전체 결합량 또는 binding burden이 더 중요할 수 있다.
- 원 논문의 transcript-level model을 지지하며, region-specific mechanism은 약하다는 결론이다.

### 결과 B: CDS composition이 추가 설명력 제공

해석:

- 같은 total CLIP amount라도 CDS 쪽에 결합이 몰린 mRNA에서 derepression이 더 크다는 의미이다.
- elongation 또는 CDS 상 ribosome movement와 연결해 논의할 수 있다.
- 단, Delta ribosome density가 CDS translation을 반영하므로 causality를 과장하지 않는다.

### 결과 C: 3-prime UTR composition이 추가 설명력 제공

해석:

- LIN28A 결합이 mRNA localization, translation regulatory complex, closed-loop regulation과 연결될 가능성을 논의한다.
- ER-associated translation model과 연결 가능하다.

### 결과 D: dose가 강함

해석:

- site 위치보다 mRNA당 결합 수 또는 binding burden이 중요할 수 있다.
- 원 논문에서 언급한 "multiple locations across CDS and 3-prime UTR"와 잘 맞는다.

### 결과 E: total amount도 재현되지 않음

해석:

- W1/Fig 4D 결과와 맞지 않으므로 pipeline/filtering 문제를 먼저 점검해야 한다.
- 이 경우 biological conclusion을 내기보다 count/normalization/filtering을 재검토한다.

### 결과 F: raw composition 효과가 adjusted composition에서 사라짐

해석:

- raw point fraction에서 보인 position effect가 region length 또는 RNA coverage 차이에 의해 생긴 confound였을 가능성이 크다.
- 최종 결론은 adjusted composition 결과를 기준으로 낸다.

---

## 9. 기존 후보와 비교한 최종 판단

| 후보 | 장점 | 약점 | 최종 판단 |
|---|---|---|---|
| Fig 2 motif 발굴 | 그림이 좋고 생물학적으로 핵심 | 원 논문 Fig 2 재현에 가까움, FASTA 필요 | 보류 |
| RNA biotype 분류 | 쉽고 안전 | Fig 3A 재현에 가까움, thesis 약함 | 제외 |
| Fig 4A/4E 재현 | 로컬 데이터로 쉬움 | 원 논문 figure 재현 중심 | 제외 |
| RPF periodicity | QC로 깔끔 | 생물학적 질문이 약함 | 제외 |
| membrane topology | ER model과 연결 | W1 Fig 5B와 가까움, annotation 부담 | 보조 후보 |
| **amount + position composition + dose-aware LIN28A binding** | 원 논문 결론 2와 3을 연결, 총량과 위치를 분리, 새 질문, 로컬 데이터로 가능 | 반복 부재와 UTR/coordinate 처리 주의 필요 | **최종 채택** |

---

## 10. 최종 topic statement

### English

This project extends Cho et al. (2012) by asking whether the amount and RNA-normalized positional composition of LIN28A binding within mRNAs explain translational derepression after Lin28a knockdown. While the original paper showed that transcript-level LIN28A CLIP enrichment correlates with ribosome density changes, this analysis separates total binding amount from regional binding preference across the 5-prime UTR, CDS, and 3-prime UTR. Regional preference will be estimated from CLIP point signals, ideally CIMS sites, normalized by region-specific RNA-control signal to reduce length and expression confounding. The central test is whether positional composition adds explanatory power beyond total CLIP amount, using nested model Delta R2, F-test/permutation testing, partial correlation, and binding-dose analysis. Because the available datapack contains one library per condition, conclusions will be framed as across-gene associations rather than per-gene differential tests.

### Korean

본 분석은 Cho et al. (2012)의 LIN28A CLIP-seq 및 ribosome profiling 데이터를 이용해, LIN28A 결합의 총량과 RNA-control로 보정한 mRNA 내부 위치 구성이 번역 억제 효과와 어떻게 연결되는지 검정한다. 원 논문은 전체 transcript 수준에서 LIN28A 결합과 ribosome density 변화의 상관을 보였지만, 본 분석은 total CLIP amount와 5-prime UTR/CDS/3-prime UTR 결합 preference를 분리하여, 같은 결합량에서도 결합 위치 구성이 Lin28a knockdown 후 번역 탈억제를 추가로 설명하는지 평가한다. 위치 신호는 가능하면 CIMS site 같은 point signal을 region별 RNA-control signal로 보정해 정의하며, `M1 = amount`와 `M2 = amount + adjusted composition`의 Delta R2를 nested F-test와 permutation으로 검정한다. 반복 라이브러리가 없는 로컬 datapack의 한계를 고려해, 결론은 gene별 유의성 검정이 아니라 gene 간 association 분석으로 제한하며, partial correlation과 dose-response 분석은 보조 해석으로 사용한다.
