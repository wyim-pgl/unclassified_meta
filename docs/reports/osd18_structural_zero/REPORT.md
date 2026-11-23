# OSD18 RNA 구조적 0 진단: Tobamovirus cacti와 Citrus yellow vein-associated virus

작성 2026-11-19. 데이터: Pronghorn `/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18/RNA`
(18 샘플 = 식물 P257·P262·P353 × healthy·macho × 반복 3, Centrifuger nr `-k 5`, `--min-score 64`, 종 수준).
진단 결과 파일: `OSD18/szdiag/struc_zero_RNA_{S,G}.tsv`, `OSD18/szdiag/reads_RNA_S.tsv` (Slurm 6279815, 8분 22초).
스크립트: 이 폴더의 `struc_zero.R`, `reads_diag.py`, `szdiag.sbatch`.

## 무엇을 했나

- **ANCOM-BC2 재적합 없음.** 저장된 표만 사용했다.
  - 구조적 0 판정: 설치된 ANCOMBC 2.14.0의 `ANCOMBC:::.get_struc_zero`(ancombc2가 `struc_zero = TRUE`일 때 부르는 함수)를
    병합 count 표에 직접 적용했다. group = `condition`, 주 분석 `neg_lb = FALSE`, 민감도 `neg_lb = TRUE`.
    대상: LCA 표(`taxa_S_lca/merged`), EM 기대값 표(`taxa_S_em/merged`), imputation 10개(`taxa_S_em/merged/imputations`).
  - 고유 리드 vs 공유 리드: taxa-table이 읽은 per-read 파일을 같은 규칙(점수 ≥ 64, Viridiplantae·Metazoa 제외, 종 수준
    투영)으로 다시 읽어, 대상 종만 후보인 리드(고유)와 다른 종과 함께 후보인 리드(공유), 공유 상대 종을 셌다.
- ANCOMBC 2.14.0 동작(소스 확인): 판정 규칙은 그룹 안 모든 샘플에서 count == 0(`present_table[present_table != 0] = 1`).
  판정된 특징은 검정 대상 표에서만 빠져 결과(`res`)에 LFC·SE·W·p가 없고, 판정은 `zero_ind`로 따로 반환된다.
  `neg_lb = TRUE`는 검출 비율의 정규근사 하한이 0 이하인 그룹도 구조적 0으로 본다(그룹당 9 샘플이라 주 분석에서 제외).

## 확인된 사실

1. **종 수준 LCA count = 고유 리드 수.** 두 바이러스의 샘플별 LCA count가 per-read 파일에서 센 고유 리드 수와
   전부 같다(예: T. cacti healthy 1, 1, 4, 214, 69, 126, 28, 19, 43). 종 수준 LCA는 후보가 한 종으로 모이는 리드만
   세기 때문이다.
2. **EM 분수 배정이 0 조건을 없앤 두 사례.**

| 분류군 | 구조적 0 조건 | LCA | EM 기대값 | imputation 10개 중 구조적 0 | 그 조건의 고유 리드 | 그 조건의 공유 리드(상대) | `neg_lb = TRUE` |
|---|---|---|---|---|---|---|---|
| Tobamovirus cacti | macho | 0 (구조적 0) | P262_M_rep1에 0.33 → 구조적 0 아님 | 7 / 10 | 0 (9 샘플) | 4 리드 (Tobamovirus muricaudae, Guangxi virgavirus, Cactus tobamovirus 1) | 모든 출처에서 macho 구조적 0 |
| Citrus yellow vein-associated virus | healthy | 0 (구조적 0) | P257_H_rep3에 0.50 → 구조적 0 아님 | 3 / 10 | 0 (9 샘플) | 10 리드, P257 healthy 3 샘플 (Opuntia umbra-like virus, Fig umbra-like virus) | 모든 출처에서 healthy 구조적 0 |

   - DA(구조적 0 검정 꺼짐): LCA에서는 두 분류군 모두 LFC NA(완전 분리). EM 기대값 fit에서는 T. cacti LFC −6.62,
     q 0.0017(robust), CYVaV LFC 7.56, q 0.035. 이 유의성은 각각 2리드·10리드의 공유 리드 분수 배정으로 0이 깨져서
     생긴 것이고, 그 조건의 고유 리드 근거는 0이다.
   - imputation에서는 구조적 0 여부가 fit마다 바뀐다(7/10, 3/10). 구조적 0인 fit에서는 계수를 추정할 수 없으므로
     MI 결합(`537f358` 이후)에서 이 특징은 모든 fit에서 추정되지 않아 p = 1이다. 구조적 0 비율은 안정성 지표로만 쓰고
     p 값과 합치지 않는다.
3. **구조적 0이 곧 생물학적 부재는 아니다.**
   - **Tobamovirus cacti (macho).** healthy에서 T. cacti 고유 리드 505개, 공유 리드 9,307개이고 공유 상대는 거의 모두
     Tobamovirus muricaudae다. 고유 리드 수의 식물·샘플 간 변화가 T. muricaudae를 따라간다(P262에서 최대). T. muricaudae
     자체는 macho에도 고유 리드 7,723개로 있다. healthy의 비율(T. cacti 고유 505 / T. muricaudae 고유 11,154,218 =
     4.5 × 10⁻⁵)을 macho의 T. muricaudae 7,723개에 적용하면 기대값은 0.35리드다. macho에서 0을 관찰한 것은 이
     검출 한계와 구별되지 않는다. 따라서 "macho에 T. cacti가 없다"는 결론은 데이터가 지지하지 않는다.
   - **Citrus yellow vein-associated virus (healthy).** 같은 umbra-like 바이러스 무리가 행 여러 개로 나뉜다.
     식물 P257의 healthy 3 샘플에는 Opuntia umbra-like virus 고유 리드 798–1,026개와 패널 행 "Citrus yellow
     vein-associated virus isolate YV920-like" 86–166리드가 있다. P257 macho에는 거의 없다(1–12). 반대로 P262·P353에서는
     macho에만 있다(Opuntia umbra-like virus 2,502–154,785, YV920-like 428–2,434). nr의 CYVaV 행만 보면 "healthy에 없음"이지만,
     같은 바이러스 무리의 리드는 P257 healthy에 있다. 이 신호는 조건보다 식물(감염 개체)에 따라 갈린다.
4. **전체 규모(종 수준, `neg_lb = FALSE`).** LCA 표에서 구조적 0인 특징 12,580개, EM 기대값 표 22,074개(EM 표가 행이 더 많다).
   LCA에서는 구조적 0인데 EM 기대값에서는 아닌 특징 2,386개. imputation 간 판정이 갈리는 특징(0 < 비율 < 1)
   6,186개, 10개 모두 구조적 0인 특징 18,828개. LCA 표는 `--panel-taxonomy` 없이 병합되어 EM 표와 패널 collapse가
   다르다(ISSUE 2d). 같은 collapse를 적용한 LCA 표(`taxa_{S,G}_lca_ptx`, Slurm 6280255, 각 약 6분)로 다시 계산해도
   LCA에서 구조적 0이지만 EM 기대값에서는 아닌 특징은 종 2,383개(속 193)이고, 반대 방향은 0개다
   (`OSD18/szdiag/struc_zero_RNA_{S,G}_ptx.tsv`). 분수 배정이 0 조건을 없애는 일은 이 두 바이러스에 한정되지 않는다.

## 해석 규칙(이 보고서에서 쓴 것)

- 구조적 0은 ANCOM-BC2 규칙(그룹 안 모든 샘플 count == 0)으로만 판정하고, EM 값을 바꾸거나 문턱을 두지 않았다.
- 출처별 판정(LCA, 고유 리드, EM 기대값, imputation 비율)을 나란히 두고, 그 조건에 고유 리드가 0이고 EM 기대값이
  공유 리드에서만 오면 "고유 리드 근거 없음; EM 값은 공유 리드 배정"으로 적는다.
- 구조적 0 특징에는 LFC·SE·W·p·q를 만들지 않는다. 구조적 0 결과와 Wald 검정 결과는 별도 표로 둔다.
- 구조적 0을 생물학적 부재로 쓰려면 같은 생물의 다른 행(관련 종, 패널 `X-like` 행)과 검출 한계를 함께 확인한다.

## 남은 것

- 위 판정을 패키지 출력으로 만드는 작업(`structural_zeros.tsv` 설명 열, imputation별 구조적 0 보존과 요약,
  `unique_counts.tsv`, `detection` 명령)은 진행 중이다.
- RNA 종 수준 MI(6272794)를 고친 코드로 재결합하고, 구조적 0이 imputation마다 갈리는 특징이 MI에서 p = 1로 남는지 확인한다.
- 같은 생물이 여러 행으로 나뉘는 문제(umbra-like 바이러스, tobamovirus 무리)는 ISSUE 2d에서 다룬다.

## v2 실행 결과 (코드 ace3467, 2026-11-20)

같은 조건(`--struc-zero`, `neg_lb = FALSE`, BH, formula `condition + plant`)으로 LCA(`taxa_*_lca_ptx`, `--panel-taxonomy`
병합), EM 기대값, MI(10 imputation, Barnard–Rubin df)를 다시 돌렸다. 출력 `OSD18/<T>/da_<R>_{lca,em}_sz`, `detection_<R>`.

| | LCA 검정 / q<0.05 / robust | EM 기대값 | MI | MI∩LCA (q<0.05) | MI만 | EM 기대값만 |
|---|---|---|---|---|---|---|
| RNA 종 | 6,786 / 16 / 8 | 9,930 / 21 / 11 | 13,138 / 3 / 3 | 3 | 0 | 7 |
| RNA 속 | 2,465 / 15 / 8 | 2,678 / 11 / 6 | 3,077 / 2 / 2 | 2 | 0 | 1 |
| DNA 종 | 10,094 / 973 / 459 | 11,633 / 1,182 / 621 | 12,573 / 486 / 287 | 475 | 11 | 374 |
| DNA 속 | 2,756 / 643 / 433 | 2,804 / 649 / 447 | 2,888 / 508 / 394 | 508 | 0 | 36 |

- RNA 종에서 EM 기대값 fit에만 유의한 7개는 Tobamovirus cacti, Citrus yellow vein-associated virus, Fig umbra-like
  virus, Sclerotinia nivalis, Alternaria incomplexa, Bottle gourd mottle virus, Watermelon green mottle mosaic virus로,
  모두 공유 리드 배정에 기대는 특징이다. MI에서는 하나도 유의하지 않다. T. cacti는 10개 fit 중 8개, CYVaV는 3개에서
  구조적 0이라 MI 상태 `incomplete`, p = 1, 결합 LFC 없음(`mi_structural_zeros.tsv`에 비율 기록).
- RNA 종 MI 유의 3개: Tobamovirus muricaudae(LFC −7.98, q 4.5 × 10⁻⁸), cyanobacterium TDX16(2.69, q 0.0018),
  Ampelomyces quisqualis(4.51, q 0.032). 셋 다 robust이고 LCA에서도 유의하다.
- MI 상태(RNA 종): pooled 3,238, untestable 5,641, partly_untestable 2,210, not_estimated 1,492, incomplete 557.
  결합된 특징의 fmi 중앙값 0.34, df 중앙값 1.04(DNA 속 0.026, 5.58).
- `detection`(RNA 종): 한 조건은 고유 리드 지지, 다른 조건은 공유 리드 배정만 있는 특징 2,382개.
- DNA의 LCA q<0.05가 이전 61개에서 973개로 늘어난 것은 BH 기본값(이전 holm), 구조적 0 특징의 검정 제외,
  `--panel-taxonomy` 재병합이 함께 바뀐 결과다.
