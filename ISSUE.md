# ISSUE — em-imputation

_최종 갱신: 2026-11-19. 작업 항목은 [TODO.md](TODO.md)._

## Blocker
- 없음.

## In Progress
1. **OSD18 실데이터(18 샘플, RNA·DNA) 실행** — Pronghorn `/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18`.
   - stage A: RNA 6272793 완료, DNA 6272798 실행 중. stage B(DA): RNA S/G em 6272794·6272796 실행 중, RNA S/G lca 6272795·6272797 완료, DNA 6272799–802 대기.
   - 이 배치는 **round-10 코드**로 실행 중. round-11 차이는 반올림 수준(≤1e-6 리드)과 안전 검사뿐.
   - 11:00 PST: RNA G em DA 6272796 완료. S em DA 6272794는 imputation fit 10개 순차 실행(fit당 약 45분), 17–18시 완료 예상. DNA stage A 6272798 실행 중.
   - MI 결합 버그(아래 Resolved)를 고친 코드로 RNA 속 수준 MI를 저장된 fit에서 재결합: `OSD18/RNA/repool_fix/mi_results_G.tsv`(고친 코드 사본 `unclassified_meta_test/em_code_fix/src`). 종 수준 S와 DNA는 잡이 끝난 뒤 같은 방식(`OSD18/repool.py <RNA|DNA> <S|G>`)으로 재결합해야 한다. 잡이 봉인한 `da_*_em/mi_results.tsv`는 고치기 전 코드의 결과다. 출처: 이 세션.
2. **시뮬레이션(오류율·CI coverage)** — 완료(alt 86·null 87 데이터셋, cabe340). 결과 `docs/reports/simulation_em_mi_20261119/`(REPORT.md, 표, results.tsv.gz; 원 실행 폴더 `runs/` 23 GB는 잡 임시폴더 `/home/wyim/.claude/jobs/d4415b8c/tmp/sim/runs`에만 있음).
   - 전역 귀무 FDR 세 방법 모두 0/87(95 % CI 0–4.2 %), raw FPR LCA 1.01 %, EM 0.98 %, MI 0.78 %. 변화 종 coverage LCA 98.6 %, EM 98.1 %, MI 99.2 %. 편향 차이는 MCSE 안.
   - 검정력 LCA 12.0 %, EM 11.4 %, **MI 5.2 %**(그룹당 6 샘플). 시뮬레이션 FMI는 0.002 수준이라 MI의 검정력 손실은 대부분 Barnard–Rubin 자유도(Resolved 2026-11-20 참조)에서 온다.
   - 고친 코드(537f358)로 저장된 fit 171개 런을 재결합: 판정·p 변화 0(per-fit p=1·W≠0 사례 없음). 시뮬레이션은 실데이터형 모호성(우세 분류군 옆 저빈도 자매종, 상위 계급 리드 집중, 조건 전용 특징, dof NA 특징)을 담지 않아 실데이터 문제를 잡지 못했다.
   - alt_085는 연장 실행 마지막에 EM 단계가 없어 실패(집계에서 빠짐).
2b. **EM 분수 배정이 완전 분리를 깸** (코드 안전장치 완료, Resolved 참조; OSD18 적용 남음) — 한 조건이 모두 0인 특징(LCA에서는 ANCOM-BC2 NA)이 EM의 0.3–0.5 리드 배정으로 유한 LFC와 q<0.05를 얻음(Tobamovirus cacti, Citrus yellow vein-associated virus). 두 방법이 같은 기준을 쓰도록 structural-zero 검정을 켜거나 조건 전용 특징을 따로 보고하는 방안 검토.
2c. **EM의 한 종 집중** — Sclerotinia nivalis 고유 리드 3개가 EM에서 329개(상위 계급 리드가 한 종으로 몰림). imputation 간 변동은 있으나 신호는 남음. S MI 결과로 판단.

## To Go
2d. **unclassified assembly 패널과 classified 배경의 중복·분할 점검** — 질문: unclassified 리드로 만든 패널 contig에 이미 분류된 생물이 들어가 리드 카운트가 틀어지는가.
   - 확인된 것(OSD18 RNA): `panel-quant` reads_scope = unclassified(mode A), `merged/reconciliation.tsv` 18 샘플 모두 ok → 한 리드가 두 번 세어지지는 않음(분류 리드와 패널 리드는 서로 겹치지 않는 집합).
   - 남은 위험: (1) 같은 생물이 배경 taxon 행과 `X-like` 패널 행으로 **나뉨** — EM 종 수준 collapse 결과 collapsed 40개(7.76 M 리드), above_rank 3개(41.3만 리드, 1개가 40.1만), no_taxid 39개(5.2만), dropped_lineage 12개(1.5만). 예: Citrus yellow vein-associated virus는 nr 행(4,180) 외에 YV920-like(9,038)·protein-only(1,189) 패널 행이 따로 남음. (2) 패널 contig의 taxon 판정이 틀리거나 chimeric contig면 다른 생물 행에 합쳐짐(collapse 7.76 M 리드의 대형 항목 확인 필요). (3) OSD18 LCA 표는 `--panel-taxonomy` 없이 병합되어 Rattail cactus necrosis virus 패널 행이 따로 남고 EM 표는 Tobamovirus muricaudae로 합쳐짐 → LCA vs EM 비교에 collapse 차이가 섞임.
   - 다음: LCA 표를 `--panel-taxonomy`로 다시 병합(taxa-table LCA ~5분)해 비교 조건 통일; above_rank·no_taxid 상위 특징과 대형 collapse 항목의 contig 판정(동일성·커버리지·chimera) 점검; 구조적 0 분석(2b)에 패널 행 포함. 출처: 이 세션(`RNA/taxa_S_em/merged/collapsed.tsv`, `quant*/provenance.json`).
3a. **da --imputations가 fit을 순차 실행** — S 수준 fit 10개 × 45분. fit 병렬화(스레드 분할) 검토.
3. **taxa-table EM이 느림** — RNA 종 수준 2시간 42분(LCA 5분). 다음: crossing 성분 SQUAREM 반복수·imputation 체인 비용 프로파일링.
4. **round-11 수정분 미검토** — Codex round 11은 "READY AFTER FIXES", 그 수정(링크된 imputations 폴더 거부, dropped_lineages 기록, round_preserving, 정수 추출)은 재검토 안 받음.

## Waiting
5. **main 병합 여부** — 사용자 결정 대기. `em-imputation`은 `p257-fixes` 위에 있음.

## Resolved
- 2026-11-20: **Barnard–Rubin 자유도(구 2e)** — 사용자 결정: 표준 공식(mice `barnard.rubin`) 유지, docstring 정정, FMI = 0에서도 같은 공식을 써서 불연속 제거(`em.rubin`). 일치하는 imputation도 ν_com(ν_com+1)/(ν_com+3) 사용. 테스트 `tests/test_em.py` 3개(수정 전 실패 확인). README에 보수성과 시뮬레이션 검정력(LCA 12 % → MI 5 %) 기록.
- 2026-11-20: **구조적 0 안전장치** — ANCOM-BC2 규칙(count != 0)은 그대로 두고 문턱·값 변경 없음. `structural_zeros.tsv`에 그룹별 samples·present_samples·reads·max_count 열; `da --imputations --struc-zero`가 fit별 `imputations/structural_zeros_NNN.tsv`와 서술용 `mi_structural_zeros.tsv`(n_fits_absent, frac, inconsistent; p 값과 결합 안 함)를 봉인; `taxa-table --assign em`이 `unique_counts.tsv`(고유 리드, 정수) 봉인; 새 명령 `detection`(LCA·고유·EM 기대값·imputation 근거와 규칙 기반 interpretation). 테스트 `tests/test_struc_zero_em.py` 8개. 전체 테스트 246 통과 + DIAMOND 정지로 오류난 `test_panel_e2e.py` 9개 재실행 통과(255 전부 통과).
- 2026-11-20: 전체 테스트가 12시간 정지 — 제3자 DIAMOND 2.1.10 `blastx` 교착(futex 대기, CPU 0 %; 상류 bbuchfink/diamond #943·#849와 같은 증상, 랩 위키 항목 없음). 해당 프로세스만 종료하고 영향받은 9개 테스트 재실행 통과. 재발 시 같은 조치.
- 2026-11-19: ccloop 런 `unclassified-meta-p257-20261118`(round 4에서 blocked)을 사용자 지시로 종료(`stop`, phase stopped). P257 목표는 이 worktree에서 수동으로 이어 진행됨. 기록 `~/.local/state/ccloop/unclassified-meta-p257-20261118/ccloop.md`.
- 2026-11-19: **MI 결합이 ANCOM-BC2의 p=1(검정 불가)을 무시함** — ANCOM-BC2는 특징의 dof가 NA면 W가 커도 p=1을 준다(`.ancombc2_core`: `p_hat[is.na(p_hat)] = 1`). `combine_imputations`는 df를 복원하지 못하면 ∞로 두고 p를 다시 계산해, 각 fit에서 검정 불가였던 특징이 유의해졌다(OSD18 RNA 속 수준 q<0.05 70개; Aedoeadaptatus는 모든 fit p=1인데 MI q 9e-57). 해결: p ≥ 1이고 W가 0이 아닌 유한값인 fit은 실패로 처리(`da.py` `combine_imputations`), 회귀 테스트 `test_a_fit_ancombc2_could_not_test_stays_untested_after_pooling`. 재결합 결과 q<0.05 70 → 2(Tobamovirus robust, Ophiobolus).
- 2026-11-19: EM 실데이터 비수렴(샘플 20, 10만 반복) — 원인: 수치 잡음 step, 평평한 방향, `dominated()` 반복 호출. 해결: 성분별 닫힌 해(Dirichlet tree), 지배 후보 제거, 인증된 likelihood gap(0.5 nat) 정지, projected SQUAREM. 샘플 20: 95초 수렴.
- 2026-11-19: 이전 OSD18 배치가 큐에 남아 있었음 — `scancel -n a,b,...`(쉼표 목록)이 아무것도 취소하지 않음. ID로 취소해 해결.
- 2026-11-19: centrifuger 단일 입력 옵션을 `-U`로 써서 실패 — 올바른 옵션은 `-u`.
- 2026-11-19: 서브에이전트(Agent) 호출이 auto mode에서 거부 — context-mode 훅이 입력을 바꿈. 대신 `codex exec`로 위임.
