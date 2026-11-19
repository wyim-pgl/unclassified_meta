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
2. **시뮬레이션(오류율·CI coverage)** — alt 5·null 3 결과는 REPORT.md에 있음(resume.md 요약). 연장 실행 중. 시뮬레이션의 FMI(중앙값 0.001)가 실데이터(0.24)보다 훨씬 작아, 실데이터형 모호성(우세 분류군 옆 저빈도 자매종, 상위 계급 리드의 한 종 집중)을 넣은 시나리오가 필요. 시뮬레이션에는 ANCOM-BC2 dof NA 특징이 없어 MI 결합 버그(Resolved, 2026-11-19)를 잡지 못했다.
2b. **EM 분수 배정이 완전 분리를 깸** — 한 조건이 모두 0인 특징(LCA에서는 ANCOM-BC2 NA)이 EM의 0.3–0.5 리드 배정으로 유한 LFC와 q<0.05를 얻음(Tobamovirus cacti, Citrus yellow vein-associated virus). 두 방법이 같은 기준을 쓰도록 structural-zero 검정을 켜거나 조건 전용 특징을 따로 보고하는 방안 검토.
2c. **EM의 한 종 집중** — Sclerotinia nivalis 고유 리드 3개가 EM에서 329개(상위 계급 리드가 한 종으로 몰림). imputation 간 변동은 있으나 신호는 남음. S MI 결과로 판단.

## To Go
3a. **da --imputations가 fit을 순차 실행** — S 수준 fit 10개 × 45분. fit 병렬화(스레드 분할) 검토.
3. **taxa-table EM이 느림** — RNA 종 수준 2시간 42분(LCA 5분). 다음: crossing 성분 SQUAREM 반복수·imputation 체인 비용 프로파일링.
4. **round-11 수정분 미검토** — Codex round 11은 "READY AFTER FIXES", 그 수정(링크된 imputations 폴더 거부, dropped_lineages 기록, round_preserving, 정수 추출)은 재검토 안 받음.

## Waiting
5. **main 병합 여부** — 사용자 결정 대기. `em-imputation`은 `p257-fixes` 위에 있음.
6. **ccloop 런 `unclassified-meta-p257-20261118`이 blocked** (round 4, "Native process exited unsuccessfully"). 원래 목표(P257 테스트 보고)는 이후 수동으로 진행됨. 재개·종료 결정 필요.

## Resolved
- 2026-11-19: **MI 결합이 ANCOM-BC2의 p=1(검정 불가)을 무시함** — ANCOM-BC2는 특징의 dof가 NA면 W가 커도 p=1을 준다(`.ancombc2_core`: `p_hat[is.na(p_hat)] = 1`). `combine_imputations`는 df를 복원하지 못하면 ∞로 두고 p를 다시 계산해, 각 fit에서 검정 불가였던 특징이 유의해졌다(OSD18 RNA 속 수준 q<0.05 70개; Aedoeadaptatus는 모든 fit p=1인데 MI q 9e-57). 해결: p ≥ 1이고 W가 0이 아닌 유한값인 fit은 실패로 처리(`da.py` `combine_imputations`), 회귀 테스트 `test_a_fit_ancombc2_could_not_test_stays_untested_after_pooling`. 재결합 결과 q<0.05 70 → 2(Tobamovirus robust, Ophiobolus).
- 2026-11-19: EM 실데이터 비수렴(샘플 20, 10만 반복) — 원인: 수치 잡음 step, 평평한 방향, `dominated()` 반복 호출. 해결: 성분별 닫힌 해(Dirichlet tree), 지배 후보 제거, 인증된 likelihood gap(0.5 nat) 정지, projected SQUAREM. 샘플 20: 95초 수렴.
- 2026-11-19: 이전 OSD18 배치가 큐에 남아 있었음 — `scancel -n a,b,...`(쉼표 목록)이 아무것도 취소하지 않음. ID로 취소해 해결.
- 2026-11-19: centrifuger 단일 입력 옵션을 `-U`로 써서 실패 — 올바른 옵션은 `-u`.
- 2026-11-19: 서브에이전트(Agent) 호출이 auto mode에서 거부 — context-mode 훅이 입력을 바꿈. 대신 `codex exec`로 위임.
