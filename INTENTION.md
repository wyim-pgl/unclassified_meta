# INTENTION — unclassified_meta (em-imputation 작업)

_최종 갱신: 2026-11-19_

## 왜 이 작업을 하는가
- Centrifuger로 분류하면 unknown(미분류) 리드가 너무 많고, 데이터베이스를 키우자니 메모리가 부족하다. 그래서 미분류 리드가 무엇인지 설명하려고 이 레포(unclassified_meta)를 만들었다.
- 사용자 목표(2026-11-18): 이 도구를 **centrifuger-quant와 합쳐서 차등 존재량(DE/DA)을 보고**, https://github.com/mourisl/centrifuger 와 https://github.com/FrederickHuangLin/ANCOMBC 를 고려해 **함께 쓸 수 있게** 하고, **DE와 중복(여러 후보에 걸린) 리드를 통계적으로 잘 처리하는 코드**를 만든다.

## 결정된 방향
- 애매한 리드(Centrifuger `-k` > 1, 패널 feature 간 동점)는 LCA나 버림 대신 **EM 최대우도로 분배**한다. 기본값(LCA, ledger 정책)은 바꾸지 않고 옵션으로 둔다(`--assign em`, `--ambiguous em`).
- 배정 불확실성은 **다중 대체(`--imputations M`) + Rubin 규칙(Barnard–Rubin df)** 으로 ANCOM-BC2 검정에 반영한다(`da --imputations`, `mi_results.tsv`).
- 같은 생물이 분류된 taxon 행과 패널 feature로 두 번 검정되지 않도록 **panel-taxonomy로 패널 서열을 같은 Centrifuger 인덱스로 분류해 taxon 행에 합친다**(`--panel-taxonomy`; 미분류 리드에 센 패널에서만).
- 차등 분석 엔진은 ANCOM-BC2 유지(조성 편향 보정, plant 공변량·랜덤효과, Eli OSD 분석과 연속성).
- 저자: Won Yim <wyim@unr.edu>, AI trailer 없음. 커밋·푸시는 사용자 요청 시에만.

## 미결정 사항
- `em-imputation` 브랜치를 main에 병합할지(사용자 결정 대기).
- 시뮬레이션(오류율·coverage) 결과에 따라 방법론 서술·기본값을 조정할지.
- 변이체(strain) 수준 분석이 필요한지(현재 potexvirus-like 4개는 같은 종으로 판정되어 합침).
- 사용자 정정: OSD 연구는 **식물 3개체(P257, P262, P353)** × healthy/macho × 3 반복. DA 공식은 `condition + plant`.
