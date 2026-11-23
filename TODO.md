# TODO — em-imputation

_최종 갱신: 2026-11-23. 관련 이슈는 [ISSUE.md](ISSUE.md), 현재 상태는 [resume.md](resume.md)._

## 다음 할 일
1. 사용자 결정: 미커밋 문서(REPORT.md 'v2 실행 결과', resume·TODO·ISSUE·INTENTION·QUARANTINE·AGENTS) 커밋, `em-imputation`의 main 병합, 푸시. → ISSUE Waiting 5
2. 새 EM 표와 이전 EM 표의 차이(칸당 최대 1.65 리드, 61행) 원인 확인. → ISSUE 3b
3. MI 결합 df가 작은 이유 확인(RNA 종 중앙값 1.04): `combine_imputations`가 복원한 완전데이터 df와 ANCOM-BC2 특징별 dof 비교. → ISSUE 3c
4. OSD18 실데이터에 robustness 확장 적용: `strain-evidence`(potexvirus 계열) → `panel-build --strain-clusters` 97/98/99 → `panel-quant --ambiguous em` → `em-uncertainty` → `da-robustness --formula 'condition + plant' --unit plant` → `da-resolutions`. Pronghorn R 환경에 ALDEx2 1.42.0·MicrobiomeStat 1.4 설치 필요(위키 installs.md 먼저 확인).
5. 패널(unclassified assembly)과 classified 배경 점검: (b) above_rank·no_taxid 상위 패널 특징과 대형 collapse 항목의 contig taxon 판정·chimera 확인, (c) 같은 생물이 배경 행과 패널 행으로 나뉜 경우 목록화(예: Citrus yellow vein-associated virus, Opuntia umbra-like virus). (a) LCA 재병합은 완료. → ISSUE 2d
6. 시뮬레이션: 실데이터형 모호성 시나리오(우세 분류군의 저빈도 자매종, 상위 계급 리드 집중, 조건 전용 특징, ANCOM-BC2 dof NA 특징) 추가. `runs/` 23 GB 보존 여부 결정. → ISSUE Waiting 6
7. `da --imputations` fit 병렬화(종 수준 8.5–9시간). → ISSUE 3a
8. Codex 재검토(round 12): round-11 수정분과 이후 변경(MI 결합, 구조적 0, detection). → ISSUE 4
9. (선택) GenBank의 Opuntia potexvirus A / Opuntia virus X 게놈을 `--extra-virus-fasta`로 추가.

## 완료
- 2026-11-20 `ace3467`: main(`integrate`)을 em-imputation에 병합, 전체 테스트 375 통과. OSD18 v2 실행(EM 표 재작성, detection, LCA·EM `--struc-zero` DA, MI) 완료 — resume.md.
- 2026-11-20 `e27d383`: 구조적 0 근거(structural_zeros 열, imputation별 구조적 0, unique_counts.tsv, detection 명령)와 보고서 2개.
- 2026-11-20 `4319f64`: Barnard–Rubin df를 FMI = 0에서도 같은 공식으로(사용자 결정).
- 2026-11-20: 종 수준 MI 재결합(고친 코드) RNA 242→2, DNA 197→33. LCA `--panel-taxonomy` 재병합.
- 2026-11-19 main `3a3c95e`: taxa-table EM 속도 개선(docs/BENCHMARK_em.md). OSD18 RNA 종 2:41:41 → 28:54.
- 2026-11-19 main `4603112`: robustness 확장(`strain-evidence`, `em-uncertainty`, `da-robustness`, `da-resolutions`).
- 2026-11-19 `537f358`: MI 결합에서 ANCOM-BC2 p=1(dof NA) fit을 실패로 처리.
- 2026-11-19 `cabe340`: EM 분배, 다중 대체, Rubin 결합, panel-taxonomy·collapse.
- 2026-11-18 `776a488`: staged 출력 권한(umask)과 da `results:` 경로 수정.
