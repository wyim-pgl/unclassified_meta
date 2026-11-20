# TODO — em-imputation

_최종 갱신: 2026-11-19. 관련 이슈는 [ISSUE.md](ISSUE.md)._

## 다음 할 일
0. 종 수준 S(6272794)·DNA EM DA가 끝나면 `OSD18/repool.py`로 고친 코드 기준 MI 재결합 → `cmp_da.py`로 LCA·EM 기대값·MI 비교 재보고.
1. OSD18 실데이터 결과 수집·보고 (RNA EM DA 6272794·6272796, DNA stage A 6272798 → DA 6272799–802). EM+MI와 LCA 결과 비교표, fmi 분포, 실행시간 정리. → ISSUE: In Progress #1
2. 시뮬레이션: 연장 실행 결과 확인. 실데이터형 모호성 시나리오(우세 분류군의 저빈도 자매종, 상위 계급 리드 집중, 조건 전용 특징, dof NA 특징) 추가. Codex 시뮬레이션 결과 검토 (`/home/wyim/.claude/jobs/d4415b8c/tmp/sim/REPORT.md`, `results.tsv`). 잡 임시폴더라 **보존 필요하면 영구 위치로 복사**. → ISSUE #2
3. (완료 2026-11-19, 커밋 전) taxa-table EM 속도 개선 → docs/BENCHMARK_em.md. 남은 일: 고친 코드를 Pronghorn `em_code/`에 반영하고 RNA G·DNA에서도 실행시간 확인.
4. round-11 수정분에 대한 Codex 재검토(round 12) 1회.
5. `em-imputation`(cabe340) PR/병합 여부를 사용자에게 확인. main 병합 시 전체 테스트 재실행.
6. (선택) potexvirus 변이체 수준 분석: contig를 97–99 % 동일성으로 재묶음 후 EM으로 정량. GenBank의 Opuntia potexvirus A / Opuntia virus X 게놈을 `--extra-virus-fasta`로 추가하는 방안.
7. (선택) ALDEx2 또는 LinDA로 DA 교차 확인.

## 완료
- 2026-11-19: MI 결합에서 ANCOM-BC2 p=1(dof NA) fit을 실패로 처리 + 회귀 테스트, RNA 속 수준 재결합(q<0.05 70 → 2).
- 2026-11-18 776a488: staged 출력 권한(umask)과 da `results:` 경로 수정 (p257-fixes, 푸시).
- 2026-11-19 cabe340: EM 분배, 다중 대체, Rubin 결합, panel-taxonomy·collapse (em-imputation, 푸시). 테스트 246 통과.
