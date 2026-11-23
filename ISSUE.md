# ISSUE — em-imputation

_최종 갱신: 2026-11-23. 작업 항목은 [TODO.md](TODO.md). 현재 상태는 [resume.md](resume.md)._

## Blocker
- 없음.

## In Progress
- 없음(2026-11-23 확인: 이 프로젝트의 실행 중 잡·에이전트 없음).

## To Go
2d. **unclassified assembly 패널과 classified 배경의 중복·분할 점검** — 질문: unclassified 리드로 만든 패널 contig에 이미 분류된 생물이 들어가 리드 카운트가 틀어지는가.
   - 확인된 것(OSD18 RNA): `panel-quant` reads_scope = unclassified(mode A), `merged/reconciliation.tsv` 18 샘플 모두 ok → 한 리드가 두 번 세어지지는 않음(분류 리드와 패널 리드는 서로 겹치지 않는 집합).
   - 남은 위험: (1) 같은 생물이 배경 taxon 행과 `X-like` 패널 행으로 **나뉨** — EM 종 수준 collapse 결과 collapsed 40개(7.76 M 리드), above_rank 3개(41.3만 리드, 1개가 40.1만), no_taxid 39개(5.2만), dropped_lineage 12개(1.5만). 예: Citrus yellow vein-associated virus는 nr 행(4,180) 외에 YV920-like(9,038)·protein-only(1,189) 패널 행이 따로 남음. (2) 패널 contig의 taxon 판정이 틀리거나 chimeric contig면 다른 생물 행에 합쳐짐(collapse 7.76 M 리드의 대형 항목 확인 필요). (3) OSD18 LCA 표는 `--panel-taxonomy` 없이 병합되어 Rattail cactus necrosis virus 패널 행이 따로 남고 EM 표는 Tobamovirus muricaudae로 합쳐짐 → LCA vs EM 비교에 collapse 차이가 섞임.
   - 2026-11-20: LCA 표를 `--panel-taxonomy`로 다시 병합(`taxa_{S,G}_lca_ptx`)해 v2 비교에 사용함(조건 통일 완료).
   - 다음: above_rank·no_taxid 상위 특징과 대형 collapse 항목의 contig 판정(동일성·커버리지·chimera) 점검; 구조적 0 분석(2b)에 패널 행 포함. 출처: 이 세션(`RNA/taxa_S_em/merged/collapsed.tsv`, `quant*/provenance.json`).
3a. **da --imputations가 fit을 순차 실행** — 종 수준 EM DA가 8.5–9시간(fit 10개 × 약 45분). fit 병렬화 검토.
3b. **새 EM 표와 이전 EM 표가 다름** — `ace3467`(main의 속도 개선 포함)으로 만든 EM 기대값이 round-10 코드 표와 칸당 최대 1.65 리드 다름(RNA 종 0.01 넘는 칸 921개, 행 31,881 vs 31,942). main의 BENCHMARK_em.md는 새 해가 최대우도에 더 가깝다고 기록. 원인(정지 규칙, 평평한 방향의 비식별 분할) 확인 필요. 출처: 2026-11-20 `OSD18/szdiag/v2check.py`.
3c. **MI 결합 df가 매우 작음** — RNA 종 결합 특징의 df 중앙값 1.04(DNA 속 5.58). 완전데이터 df 복원값(ANCOM-BC2 특징별 dof)이 희소 특징에서 작은 것으로 보이나 확인하지 않음.
4. **round-11 수정분 미검토** — Codex round 11은 "READY AFTER FIXES", 그 수정(링크된 imputations 폴더 거부, dropped_lineages 기록, round_preserving, 정수 추출)은 재검토 안 받음.

## Waiting
5. **main 병합·푸시 여부** — 사용자 결정 대기. `em-imputation`(`ace3467`)은 main(`a547229`)을 포함하고 `4319f64`·`e27d383`·`ace3467`이 더 있다. 둘 다 원격보다 앞서 있음(푸시 안 함).
6. **시뮬레이션 `runs/` 23 GB 보존 여부** — 잡 임시폴더에만 있어 잡 정리 시 사라짐. 요약·스크립트·결과표는 `docs/reports/simulation_em_mi_20261119/`에 보존됨.

## Resolved
- 2026-11-20: **OSD18 실데이터 실행(구 In Progress #1)** — v2 실행 완료(resume.md, REPORT.md 'v2 실행 결과'). 이전 배치와 고친 코드 재결합 결과는 `OSD18/<T>/repool_fix/`.
- 2026-11-20: **EM 분수 배정이 완전 분리를 깸(구 2b)** — 안전장치 적용 결과: RNA 종에서 EM 기대값에만 유의한 7개가 MI에서 모두 유의하지 않음; T. cacti는 imputation 8/10, CYVaV 3/10에서 구조적 0 → MI `incomplete`, p = 1.
- 2026-11-20: **EM의 한 종 집중(구 2c)** — Sclerotinia nivalis는 EM 기대값에서만 유의, MI에서 유의하지 않음.
- 2026-11-20: **시뮬레이션(구 2)** — 완료, `docs/reports/simulation_em_mi_20261119/`. 결과 요약은 resume.md.
- 2026-11-20: **Barnard–Rubin 자유도(구 2e)** — 사용자 결정: 표준 공식(mice `barnard.rubin`) 유지, docstring 정정, FMI = 0에서도 같은 공식을 써서 불연속 제거(`em.rubin`). 일치하는 imputation도 ν_com(ν_com+1)/(ν_com+3) 사용. 테스트 `tests/test_em.py` 3개(수정 전 실패 확인). README에 보수성과 시뮬레이션 검정력(LCA 12 % → MI 5 %) 기록.
- 2026-11-20: **구조적 0 안전장치** — ANCOM-BC2 규칙(count != 0)은 그대로 두고 문턱·값 변경 없음. `structural_zeros.tsv`에 그룹별 samples·present_samples·reads·max_count 열; `da --imputations --struc-zero`가 fit별 `imputations/structural_zeros_NNN.tsv`와 서술용 `mi_structural_zeros.tsv`(n_fits_absent, frac, inconsistent; p 값과 결합 안 함)를 봉인; `taxa-table --assign em`이 `unique_counts.tsv`(고유 리드, 정수) 봉인; 새 명령 `detection`(LCA·고유·EM 기대값·imputation 근거와 규칙 기반 interpretation). 테스트 `tests/test_struc_zero_em.py` 8개. 전체 테스트 246 통과 + DIAMOND 정지로 오류난 `test_panel_e2e.py` 9개 재실행 통과(255 전부 통과).
- 2026-11-20: 전체 테스트가 12시간 정지 — 제3자 DIAMOND 2.1.10 `blastx` 교착(futex 대기, CPU 0 %; 상류 bbuchfink/diamond #943·#849와 같은 증상, 랩 위키 항목 없음). 해당 프로세스만 종료하고 영향받은 9개 테스트 재실행 통과. 재발 시 같은 조치.
- 2026-11-19: **taxa-table EM이 느림**(RNA 종 수준 2시간 42분, LCA 5분) — 원인: (1) 접평면 gap ≤ 0.5 정지 규칙이 작은 몫의 target 하나 때문에 최대 147,634회 반복(로그우도는 1e-4 수준에서 이미 정지), (2) imputation의 SIR 시작(트리 재귀 256회/체인)과 Gibbs 배정(np.add.at, 체인·성분별 순차). 해결: 라그랑주 쌍대 인증(_dual_gap, 그룹당 1e-4 nat), 배치 Dirichlet-tree 추출, 모든 crossing 성분·M 체인 동시 진행, 클래스 단위 다항 추출. 18 샘플 EM 3,299→432 s, imputation 4,300→551 s; 기대값은 기준보다 최대우도에 더 가까움(+0.17–1.88 nat), imputation 분포 동등(M=200, 3 샘플). 보고서 docs/BENCHMARK_em.md.
- 2026-11-20: **다중검정 보정을 BH로** — 사용자 결정. `da --p-adjust` 기본값을 holm → BH로 바꿈(CLI, `ancombc2.R`). 속 수준 재결합 BH 결과 `OSD18/RNA/repool_fix/mi_results_G_BH.tsv`: q<0.05 5개(Tobamovirus·Ampelomyces robust, Ophiobolus, Potexvirus·Carlavirus 0.04999). 저장된 fit의 diff_robust는 holm 기준이라 robust 판정은 보수적.
- 2026-11-19: **MI 결합 감사(audit)** — 537f358의 판정은 W가 NA이거나 W 열이 없으면 p=1(검정 불가) fit을 유효로 보고 다시 검정해 유의하게 만들었다(합성 예: p 2e-9). 해결: W가 없으면 lfc/se로 판단, |W| ≤ 1e-12만 진짜 p=1, df는 집계된 fit에서만 복원, fit 2개 미만·계수 열 불일치는 거부, `status_`/`n_untestable_` 열 추가. 저장된 OSD18 RNA 속 fit 재결합은 기준(`repool_fix/mi_results_G.tsv`)과 4e-16 이내로 같음. 기준의 q는 **Holm**이다(BH 아님; BH면 q<0.05 5개). `fmi`는 λ다. 상세 `docs/MI_POOLING_AUDIT.md`, 테스트 `tests/test_mi_pooling.py`.
- 2026-11-19: ccloop 런 `unclassified-meta-p257-20261118`(round 4에서 blocked)을 사용자 지시로 종료(`stop`, phase stopped). P257 목표는 이 worktree에서 수동으로 이어 진행됨. 기록 `~/.local/state/ccloop/unclassified-meta-p257-20261118/ccloop.md`.
- 2026-11-19: **MI 결합이 ANCOM-BC2의 p=1(검정 불가)을 무시함** — ANCOM-BC2는 특징의 dof가 NA면 W가 커도 p=1을 준다(`.ancombc2_core`: `p_hat[is.na(p_hat)] = 1`). `combine_imputations`는 df를 복원하지 못하면 ∞로 두고 p를 다시 계산해, 각 fit에서 검정 불가였던 특징이 유의해졌다(OSD18 RNA 속 수준 q<0.05 70개; Aedoeadaptatus는 모든 fit p=1인데 MI q 9e-57). 해결: p ≥ 1이고 W가 0이 아닌 유한값인 fit은 실패로 처리(`da.py` `combine_imputations`), 회귀 테스트 `test_a_fit_ancombc2_could_not_test_stays_untested_after_pooling`. 재결합 결과 q<0.05 70 → 2(Tobamovirus robust, Ophiobolus).
- 2026-11-19: EM 실데이터 비수렴(샘플 20, 10만 반복) — 원인: 수치 잡음 step, 평평한 방향, `dominated()` 반복 호출. 해결: 성분별 닫힌 해(Dirichlet tree), 지배 후보 제거, 인증된 likelihood gap(0.5 nat) 정지, projected SQUAREM. 샘플 20: 95초 수렴.
- 2026-11-19: 이전 OSD18 배치가 큐에 남아 있었음 — `scancel -n a,b,...`(쉼표 목록)이 아무것도 취소하지 않음. ID로 취소해 해결.
- 2026-11-19: centrifuger 단일 입력 옵션을 `-U`로 써서 실패 — 올바른 옵션은 `-u`.
- 2026-11-19: 서브에이전트(Agent) 호출이 auto mode에서 거부 — context-mode 훅이 입력을 바꿈. 대신 `codex exec`로 위임.
