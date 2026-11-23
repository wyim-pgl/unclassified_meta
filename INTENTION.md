# INTENTION — unclassified_meta (em-imputation 작업)

_최종 갱신: 2026-11-23_

## 왜 이 작업을 하는가
- Centrifuger로 분류하면 unknown(미분류) 리드가 너무 많고, 데이터베이스를 키우자니 메모리가 부족하다. 그래서 미분류 리드가 무엇인지 설명하려고 이 레포(unclassified_meta)를 만들었다.
- 사용자 목표(2026-11-18): 이 도구를 **centrifuger-quant와 합쳐서 차등 존재량(DE/DA)을 보고**, https://github.com/mourisl/centrifuger 와 https://github.com/FrederickHuangLin/ANCOMBC 를 고려해 **함께 쓸 수 있게** 하고, **DE와 중복(여러 후보에 걸린) 리드를 통계적으로 잘 처리하는 코드**를 만든다.

## 결정된 방향
- 애매한 리드(Centrifuger `-k` > 1, 패널 feature 간 동점)는 LCA나 버림 대신 **EM 최대우도로 분배**한다. 기본값(LCA, ledger 정책)은 바꾸지 않고 옵션으로 둔다(`--assign em`, `--ambiguous em`).
- 배정 불확실성은 **다중 대체(`--imputations M`) + Rubin 규칙(Barnard–Rubin df)** 으로 ANCOM-BC2 검정에 반영한다(`da --imputations`, `mi_results.tsv`).
- 같은 생물이 분류된 taxon 행과 패널 feature로 두 번 검정되지 않도록 **panel-taxonomy로 패널 서열을 같은 Centrifuger 인덱스로 분류해 taxon 행에 합친다**(`--panel-taxonomy`; 미분류 리드에 센 패널에서만).
- 구조적 0(2026-11-20 사용자 요구): ANCOM-BC2 자체 규칙(그룹 안 모든 샘플 count == 0, `neg_lb = FALSE`가 주 분석)을 그대로 쓰고, EM 분수 값을 바꾸거나 문턱을 두지 않는다. 대신 고유 리드·LCA·EM 기대값·imputation 근거를 나란히 보고(`detection`), 구조적 0 특징에는 LFC·SE·p를 만들지 않으며, imputation별 구조적 0 비율은 서술 지표로만 두고 Rubin 결합이나 p 값에 넣지 않는다. LCA의 구조적 0도 검출 한계·같은 생물의 다른 행을 확인하기 전에는 생물학적 부재로 쓰지 않는다.
- 다중검정 보정 기본값 BH(2026-11-20 사용자 결정, main `06a4ecd`).
- Barnard–Rubin 자유도: 표준 공식(mice) 유지, imputation이 일치해도 같은 공식(2026-11-20 사용자 결정, `4319f64`).
- 차등 분석 엔진은 ANCOM-BC2 유지(조성 편향 보정, plant 공변량·랜덤효과, Eli OSD 분석과 연속성).
- 저자: Won Yim <wyim@unr.edu>, AI trailer 없음. 커밋·푸시는 사용자 요청 시에만.

## 미결정 사항
- `em-imputation`(`ace3467`)을 main에 병합할지, 두 브랜치를 푸시할지(사용자 결정 대기).
- 시뮬레이션 원 실행 폴더(23 GB) 보존 여부.
- 변이체(strain) 수준 분석을 OSD18에 적용할지(도구는 main `4603112`로 구현됨; 현재 potexvirus-like 4개는 같은 종으로 합침).
- 사용자 정정: OSD 연구는 **식물 3개체(P257, P262, P353)** × healthy/macho × 3 반복. DA 공식은 `condition + plant`.
