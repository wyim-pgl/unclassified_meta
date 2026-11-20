# resume — 현재 상태 (2026-11-19, 약 11:00 PST 기준)

관련 기록: [INTENTION.md](INTENTION.md) · [TODO.md](TODO.md) · [ISSUE.md](ISSUE.md) · [QUARANTINE.md](QUARANTINE.md) · [AGENTS.md](AGENTS.md)

## 위치
- 작업 worktree: `/mnt/f/scratch/develop/unclassified_meta-p257fix` (브랜치 `em-imputation`, HEAD `cabe340`, 푸시됨, main 미병합). 이 핸드오프 파일 7개는 **미커밋**.
- 브랜치: main `57448db` ← `p257-fixes` `776a488`(푸시) ← `em-imputation` `cabe340`(푸시).
- 메인 checkout `/mnt/f/scratch/develop/unclassified_meta`: main, 미커밋 `wiki/`·`.ignore`(이 작업과 무관). ccloop 런 `unclassified-meta-p257-20261118`은 2026-11-19 사용자 지시로 종료(stopped).

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

## DNA 비교 (2026-11-19 14:30 PST; stage A 6272798 완료 2:13:30, DA LCA S·G와 EM G 완료, EM S 6272799 실행 중)
- 종(S): LCA 12,612 / q<0.05 61 / robust 26; EM 기대값 14,054 / 54 / 32; LFC Pearson 0.901(양쪽 추정 9,781); robust 공통 24.
- 속(G): LCA 3,255 / 71 / 24; EM 기대값 3,243 / 105 / 27; LFC Pearson 0.983; robust 공통 22.
- 속(G) MI 재결합(고친 코드, `OSD18/DNA/repool_fix/mi_results_G.tsv`): q<0.05 40(수정 전 75), 40개 모두 LCA에서도 q<0.05; robust 25(LCA와 공통 21, MI만 4, LCA만 3). fmi 중앙값 0.029, 90분위 0.25 — RNA(0.21)보다 배정 불확실성이 작다.

## 종 수준 MI 재결합 (2026-11-20 02:40 PST; RNA EM DA 6272794 7:43:16, DNA EM DA 6272799 9:07:48 완료, 고친 코드로 `OSD18/<T>/repool_fix/mi_results_S.tsv`)
- RNA 종: 수정 전 MI q<0.05 242 → **2**(Tobamovirus muricaudae LFC −7.99, cyanobacterium TDX16 LFC 2.68; 둘 다 robust, 둘 다 LCA에서도 유의). LCA 7/robust 5, EM 기대값 10/4. fmi 중앙값 0.34, 90분위 0.82.
- DNA 종: 수정 전 197 → **33**(LCA와 공통 32, MI만 1); robust 17(전부 LCA robust와 공통, LCA만 9). LCA 61/26, EM 기대값 54/32. fmi 중앙값 0.15, 90분위 0.61.
- MI는 LCA보다 적게 부르고, MI가 부른 것은 거의 모두 LCA도 부른다. EM 기대값 fit에만 있는 유의 항목은 MI에서 남지 않는다.

## 구조적 0 진단 (RNA; 보고서 `docs/reports/osd18_structural_zero/REPORT.md`)
- ANCOMBC 2.14.0 `.get_struc_zero`(group condition)를 저장된 표에 직접 적용(재적합 없음), per-read 파일로 고유·공유 리드 계산(Slurm 6279815).
- 종 수준 LCA count = 고유 리드 수. Tobamovirus cacti(macho)·Citrus yellow vein-associated virus(healthy)는 그 조건의 고유 리드 0, 공유 리드 4·10개의 분수 배정이 EM 기대값의 0을 없앰; imputation 중 구조적 0 7/10·3/10.
- 두 경우 모두 구조적 0을 부재로 읽을 근거가 없다: T. cacti macho 기대 리드 0.35(검출 한계 안), CYVaV는 같은 umbra-like 무리가 P257 healthy에 있음(식물별 감염).
- 전체(종, `--panel-taxonomy`로 재병합한 LCA `taxa_S_lca_ptx` 기준): LCA 구조적 0이지만 EM 기대값은 아닌 특징 2,383개(속 193), 반대 방향 0개, imputation 간 판정이 갈리는 특징 6,186개(속 617). 재병합 전 LCA와 차이는 3개뿐.

## 시뮬레이션
- 폴더 `/home/wyim/.claude/jobs/d4415b8c/tmp/sim/`(잡 임시폴더 — 끝나면 영구 위치로 복사). REPORT.md는 alt 5·null 3 데이터셋 기준. 연장 실행(alt 12·null 8 목표, 4시간 예산) 진행 중.
- 결과(alt 5·null 3): 세 방법 편향 유사(changed LFC 편향 LCA 0.091, EM 0.081, MI 0.083), 변화 없는 종 coverage LCA 97.1 %, EM 96.5 %, MI 99.4 %, 검정력 LCA 10 %, EM 6.7 %, MI 0 %(그룹당 6 샘플), global-null FDR 모두 0. 시뮬레이션 FMI 중앙값 0.001로 실데이터(0.24)보다 훨씬 작음 → 시뮬레이션의 모호성 구조가 실데이터를 대표하지 못함.

## 다음 단계
1. S em DA(6272794)·DNA DA 완료 후 `OSD18/repool.py <RNA|DNA> <S|G>`로 고친 코드 기준 MI 재결합, `cmp_da.py`로 재비교.
2. MI 수정 커밋 후 Pronghorn `em_code/`를 갱신해 이후 잡이 고친 코드를 쓰게 한다(실행 중 잡이 끝난 뒤).
3. 시뮬레이션 연장 결과 확인·복사, TODO.md 나머지.
