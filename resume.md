# resume — 현재 상태 (2026-11-19, 약 11:00 PST 기준)

관련 기록: [INTENTION.md](INTENTION.md) · [TODO.md](TODO.md) · [ISSUE.md](ISSUE.md) · [QUARANTINE.md](QUARANTINE.md) · [AGENTS.md](AGENTS.md)

## 위치
- 작업 worktree: `/mnt/f/scratch/develop/unclassified_meta-p257fix` (브랜치 `em-imputation`, HEAD `cabe340`, 푸시됨, main 미병합). 이 핸드오프 파일 7개는 **미커밋**.
- 브랜치: main `57448db` ← `p257-fixes` `776a488`(푸시) ← `em-imputation` `cabe340`(푸시).
- 메인 checkout `/mnt/f/scratch/develop/unclassified_meta`: main, 미커밋 `wiki/`·`.ignore`(이 작업과 무관). ccloop 런이 이 폴더를 감시하다 blocked 상태.

## 한 일 (cabe340)
- `src/unclassified_meta/em.py`: 성분별 EM(laminar = Dirichlet tree 닫힌 해, crossing = projected SQUAREM, 인증 gap ≤ 0.5 nat), 지배 후보 제거, imputation(laminar 정확 추출, crossing = SIR 시작 Gibbs 500 sweep), Rubin/Barnard–Rubin 결합, t 꼬리 계산, round_preserving.
- taxa-table `--assign em --imputations --seed`, panel-quant `--ambiguous em --imputations --seed`, merge의 imputation 결합·검증·collapse, `panel-taxonomy` 명령, `da --imputations`(mi_results.tsv).
- README 섹션 추가, 테스트 파일 9개(test_em*.py, test_panel_taxonomy.py, test_review_em1/2.py).

## 검증
- 전체 테스트: 246 passed (최종 코드, 2026-11-19, R 포함 실행; 명령은 AGENTS.md의 환경 설정).
- Codex 독립 검토 11회. r10·r11 "READY AFTER FIXES" → 수정 반영. r11 수정분은 재검토 안 함.
- 실데이터 RNA 샘플 20(15.9 M 리드): 성분 1,586개 중 laminar 1,541개, 최대 crossing 성분 5,156 target. EM 95초 수렴(36,678 반복), imputation 10개 47.6초.

## 실데이터 OSD18 (Pronghorn `/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18`, round-10 코드 `em_code/`)
- RNA stage A 6272793 완료: taxa-S-em 2:41:41, taxa-G-em 1:05:31, LCA 각 ~5분, centrifuger(nr) 패널 분류 3.5분 202 GB.
- 결과: 18 샘플 reconciliation 전부 ok(차이 0). 종 수준 collapse 40개(potexvirus-like 4개 → Opuntia potexvirus A). 종 수준 kept: EM 58.38 M vs LCA 53.13 M, above-rank 0.43 M vs 5.75 M.
- LCA DA(종): 8,806 검정, q<0.05 7, robust 5 (Tobamovirus muricaudae, Rattail cactus necrosis virus 관련 2개 macho↓; cyanobacterium TDX16, Ampelomyces quisqualis macho↑).
- 11:00 PST: RNA G em DA 6272796 완료(2:00:13, maxRSS 0.5 GB). RNA S em DA 6272794 실행 중 — imputation 10개 fit을 **순차** 실행, S 수준 fit 1개 약 45분(기대값 fit 50분) → 완료 예상 17–18시 PST. DNA stage A 6272798 실행 중, DNA DA 6272799–802 대기.

## LCA vs EM 비교 (RNA, 2026-11-19 11:00; 스크립트 Pronghorn `OSD18/cmp_da.py`, `cmp3.py`, `cmp4.py`)
- 종(S): LCA 검정 8,806 / q<0.05 7 / robust 5; EM 기대값 fit 12,900 / 10 / 4 (staging에서 읽음).
  - 일치: Tobamovirus muricaudae(macho↓, LFC −8.0), cyanobacterium TDX16, Ampelomyces quisqualis(robust), Ophiobolus disseminans.
  - LCA의 Rattail cactus necrosis associated virus 패널 행 2개(robust)는 EM에서 Tobamovirus muricaudae로 collapse됨(같은 종, ICTV 이명).
  - EM에만 나온 것: Tobamovirus cacti·Citrus yellow vein-associated virus — LCA에서는 한 조건 전체가 0이라 ANCOM-BC2가 NA, EM은 0.33·0.50 리드의 분수 배정이 완전 분리를 깨서 유한 LFC가 나옴. Sclerotinia nivalis(LCA 3 리드 → EM 329, P353 macho rep2·3에 집중), Alternaria incomplexa(LFC 부호 반대), Penicillium malachiteum, Alternaria arbusti.
- 속(G): LCA 2,853 / 7 / 4; EM 기대값 3,081 / 5 / 2; 양쪽 추정 가능 2,295개 LFC Pearson 0.893.
- 속(G) EM-MI(잡이 봉인한 결과, 고치기 전 코드): 3,486 / q<0.05 70 / robust 1 — MI 결합 버그(ISSUE Resolved 2026-11-19).
- 속(G) EM-MI 재결합(고친 코드, `OSD18/RNA/repool_fix/mi_results_G.tsv`): 결합된 특징 1,656, q<0.05 **2**(Tobamovirus robust, LFC −8.09, fmi 0.07; Ophiobolus LFC 3.91, q 0.038), fmi 중앙값 0.21, 90분위 0.70.
  - LCA 유의 항목의 MI 값: Ampelomyces q 0.059(fmi 0.08, robust_frac 1.0), Carlavirus q 0.25, Parastagonospora q 1(fmi 0.66, df 2.7), Neocucurbitaria q 1(fmi 0.79). 배정에 따라 결과가 바뀌는 속은 MI에서 유의하지 않다.
- potexvirus-like 4개 점검(`RNA/potex_check/`): 특징 간 동일성 92–98 %, Spearman 0.88–0.99 → 같은 종.

## 시뮬레이션
- 폴더 `/home/wyim/.claude/jobs/d4415b8c/tmp/sim/`(잡 임시폴더 — 끝나면 영구 위치로 복사). REPORT.md는 alt 5·null 3 데이터셋 기준. 연장 실행(alt 12·null 8 목표, 4시간 예산) 진행 중.
- 결과(alt 5·null 3): 세 방법 편향 유사(changed LFC 편향 LCA 0.091, EM 0.081, MI 0.083), 변화 없는 종 coverage LCA 97.1 %, EM 96.5 %, MI 99.4 %, 검정력 LCA 10 %, EM 6.7 %, MI 0 %(그룹당 6 샘플), global-null FDR 모두 0. 시뮬레이션 FMI 중앙값 0.001로 실데이터(0.24)보다 훨씬 작음 → 시뮬레이션의 모호성 구조가 실데이터를 대표하지 못함.

## 다음 단계
1. S em DA(6272794)·DNA DA 완료 후 `OSD18/repool.py <RNA|DNA> <S|G>`로 고친 코드 기준 MI 재결합, `cmp_da.py`로 재비교.
2. MI 수정 커밋 후 Pronghorn `em_code/`를 갱신해 이후 잡이 고친 코드를 쓰게 한다(실행 중 잡이 끝난 뒤).
3. 시뮬레이션 연장 결과 확인·복사, TODO.md 나머지.
