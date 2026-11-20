# TODO — em-imputation

_최종 갱신: 2026-11-19. 관련 이슈는 [ISSUE.md](ISSUE.md)._

## 다음 할 일
0a. OSD18에 `detection` 적용: EM taxa 폴더를 이 코드로 다시 만들어야 `unique_counts.tsv`가 생김(RNA 종 2h42m·속 1h05m). 다시 만들 때 `da --struc-zero`(LCA·EM)도 함께 실행.
0. 종 수준 S(6272794)·DNA EM DA가 끝나면 `OSD18/repool.py`로 고친 코드 기준 MI 재결합 → `cmp_da.py`로 LCA·EM 기대값·MI 비교 재보고.
0b. OSD18 실데이터에 robustness 확장 적용: `strain-evidence`(potexvirus 계열) → `panel-build --strain-clusters` 97/98/99 → `panel-quant --ambiguous em` → `em-uncertainty` → `da-robustness --formula 'condition + plant' --unit plant` → `da-resolutions`. Pronghorn R 환경에 ALDEx2 1.42.0·MicrobiomeStat 1.4 설치 필요(위키 installs.md 확인 후).
1. OSD18 실데이터 결과 수집·보고 (RNA EM DA 6272794·6272796, DNA stage A 6272798 → DA 6272799–802). EM+MI와 LCA 결과 비교표, fmi 분포, 실행시간 정리. → ISSUE: In Progress #1
2. 시뮬레이션: 결과 보존 완료(`docs/reports/simulation_em_mi_20261119/`). 남은 것: 실데이터형 모호성 시나리오(우세 분류군의 저빈도 자매종, 상위 계급 리드 집중, 조건 전용 특징, dof NA 특징) 추가; `runs/` 23 GB 보존 여부 결정.
3. (완료 2026-11-19, 커밋 전) taxa-table EM 속도 개선 → docs/BENCHMARK_em.md. 남은 일: 고친 코드를 Pronghorn `em_code/`에 반영하고 RNA G·DNA에서도 실행시간 확인.
4. round-11 수정분에 대한 Codex 재검토(round 12) 1회.
5. `em-imputation`(cabe340) PR/병합 여부를 사용자에게 확인. main 병합 시 전체 테스트 재실행.
6. (선택) GenBank의 Opuntia potexvirus A / Opuntia virus X 게놈을 `--extra-virus-fasta`로 추가하는 방안 (변이체 수준 분석 자체는 `strain-evidence`로 구현됨).

8. 패널(unclassified assembly)과 classified 배경 점검 → ISSUE 2d: (a) OSD18 LCA 표를 `--panel-taxonomy`로 재병합해 EM과 collapse 조건 통일 후 비교 재실행, (b) above_rank·no_taxid 상위 패널 특징과 대형 collapse 항목의 contig taxon 판정·chimera 확인, (c) 같은 생물이 배경 행과 패널 행으로 나뉜 경우 목록화(예: Citrus yellow vein-associated virus).

## 완료
- 2026-11-20: 구조적 0 안전장치(structural_zeros 열, imputation별 구조적 0, unique_counts.tsv, detection 명령), 테스트 255 통과.
- 2026-11-19: MI 결합에서 ANCOM-BC2 p=1(dof NA) fit을 실패로 처리 + 회귀 테스트, RNA 속 수준 재결합(q<0.05 70 → 2).
- 2026-11-19 4603112 (worktree-additional): robustness 확장 — `strain-evidence`, `panel-build --strain-clusters`, panel-quant의 `equivalence_classes.tsv`/`unique_counts.tsv`, `em-uncertainty`, `da-robustness`(ANCOM-BC2·ALDEx2·LinDA), `da-resolutions`. 설계 `docs/DESIGN_robustness.md`, README 10절. 이전 TODO 6(변이체 분석)·7(ALDEx2/LinDA 교차 확인) 포함.
- 2026-11-18 776a488: staged 출력 권한(umask)과 da `results:` 경로 수정 (p257-fixes, 푸시).
- 2026-11-19 cabe340: EM 분배, 다중 대체, Rubin 결합, panel-taxonomy·collapse (em-imputation, 푸시). 테스트 246 통과.
