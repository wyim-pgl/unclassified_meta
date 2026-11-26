# resume — 현재 상태 (2026-11-25 기준; 아래 v2·시뮬레이션 절은 2026-11-23 확인분)

관련 기록: [INTENTION.md](INTENTION.md) · [TODO.md](TODO.md) · [ISSUE.md](ISSUE.md) · [QUARANTINE.md](QUARANTINE.md) · [AGENTS.md](AGENTS.md) · [README.md](README.md)

## 2026-11-25 갱신
- `origin/main` = `734c77e`(em-imputation `ace3467` 포함, 푸시됨). worktree `additional`(브랜치 `integrate`)도 이 위치로 fast-forward.
- 전체 테스트 366 passed(`7c6ae6d`, 2026-11-24, 29분, 건너뛴 것 없음).
- OSD18 robustness 확장 실행 중: ISSUE In Progress 1(host genome 경로 변경으로 첫 제출 실패 → 재제출 6331075–82).
- 관련 별도 저장소: Centrifuger lab fork `plantgenomicslab/centrifuger`(tag `v1.2.0-pgl.1`, 상류 PR mourisl/centrifuger#114) — `centrifuger-quant`의 EM read count(emNumReads, `--em-reads`), protein index weight, kraken-report 합계 수정. 위키 `guide/centrifuger.md`.
- Eli의 Opuntia microbiome pilot 데이터(`/data/gpfs/assoc/pgl/Eli/opuntia_microbiome/pilot_data/`, 원본 NL4CGS 40 샘플 + rerun KYZMNB 4 샘플) QC 진행 중: 출력 `unclassified_meta_test/eli_pilot/`.

## 위치와 git (2026-11-23 기록)
- 작업 worktree: `/mnt/f/scratch/develop/unclassified_meta-p257fix`, 브랜치 `em-imputation`, HEAD `ace3467`.
  - `origin/em-imputation`은 `cabe340`에 있다(로컬이 18 커밋 앞섬, **푸시 안 함**).
  - 미커밋: `TODO.md`, `resume.md`, `docs/reports/osd18_structural_zero/REPORT.md`('v2 실행 결과' 절), 이 핸드오프 갱신분.
- 메인 checkout `/mnt/f/scratch/develop/unclassified_meta`: `main` `a547229`(사용자가 `integrate`를 병합; `origin/main`보다 17 커밋 앞섬, 푸시 안 함). 미추적 `wiki/`·`.ignore`는 이 작업과 무관.
- 커밋 흐름: `cabe340`(EM·MI·panel-taxonomy) → `aea3ae7`(핸드오프) → `537f358`(MI p=1 수정) → `e190ef0`(ccloop 종료 기록) → main의 `integrate`(EM 속도 개선, BH 기본값, MI 결합 감사, strains, uncertainty, robustness) → `4319f64`(Barnard–Rubin df) → `e27d383`(구조적 0 근거·`detection`·보고서) → `ace3467`(main을 em-imputation에 병합).
- `em-imputation`은 main의 모든 커밋을 포함한다. main에는 `4319f64`·`e27d383`·`ace3467`이 아직 없다.

## 검증
- 전체 테스트(`ace3467`, 로컬): **375 passed**, 35분 36초(2026-11-20). 명령은 AGENTS.md.
- 병합 전 대상 테스트 96 passed; Barnard–Rubin 테스트는 수정 전 실패·수정 후 통과 확인.

## 실데이터 OSD18 v2 (완료, 2026-11-20; Pronghorn `/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18`, 코드 `../em_code_v2/` = `ace3467`)
- 잡(모두 COMPLETED): 1단계 RNA 6302073(1:31:01)·DNA 6302076(1:22:49); EM DA RNA 종 6302074(8:29:16)·속 6302075(1:56:00), DNA 종 6302077(9:02:31)·속 6302078(2:02:40). 스크립트 `OSD18/v2A.sbatch`, `v2B.sbatch`, 로그 `OSD18/logs/v2*_*.log`.
- 출력(각 `OSD18/<RNA|DNA>/`): `taxa_{S,G}_em_v2`(unique_counts.tsv 포함), `taxa_{S,G}_lca_ptx`(`--panel-taxonomy`), `detection_{S,G}`, `da_{S,G}_lca_sz`, `da_{S,G}_em_sz`(mi_results.tsv, mi_structural_zeros.tsv).
- 결과표와 해석: `docs/reports/osd18_structural_zero/REPORT.md` 'v2 실행 결과'. 요약(q<0.05 / robust, BH, `--struc-zero`):

| | LCA | EM 기대값 | MI | MI∩LCA |
|---|---|---|---|---|
| RNA 종 | 16 / 8 | 21 / 11 | 3 / 3 | 3 |
| RNA 속 | 15 / 8 | 11 / 6 | 2 / 2 | 2 |
| DNA 종 | 973 / 459 | 1,182 / 621 | 486 / 287 | 475 |
| DNA 속 | 643 / 433 | 649 / 447 | 508 / 394 | 508 |

  - RNA 종 MI 3개: Tobamovirus muricaudae(macho↓), cyanobacterium TDX16, Ampelomyces quisqualis(macho↑). EM 기대값에만 유의한 7개(T. cacti, Citrus yellow vein-associated virus, Fig umbra-like virus, Sclerotinia nivalis, Alternaria incomplexa, Bottle gourd mottle virus, Watermelon green mottle mosaic virus)는 공유 리드 배정에 의존하며 MI에서 모두 유의하지 않다.
  - EM 표 작성 시간: RNA 종 2:41:41 → 28:54(main의 속도 개선). `unique_counts.tsv` = LCA count(네 표 모든 칸). 새 EM 표는 이전 EM 표와 칸당 최대 1.65 리드 다름(TODO 9).
- 이전 배치(round-10 코드 `em_code/`, holm, 구조적 0 검정 꺼짐)의 `da_*_em/mi_results.tsv`는 MI 결합 버그 이전 결과라 쓰지 않는다. 고친 코드로 재결합한 것은 `OSD18/<T>/repool_fix/`.

## 시뮬레이션 (완료)
- 결과 `docs/reports/simulation_em_mi_20261119/`(alt 86·null 87, cabe340): 전역 귀무 FDR 0/87(세 방법), 검정력 LCA 12.0 %, EM 11.4 %, MI 5.2 %(그룹당 6 샘플; MI 손실은 Barnard–Rubin df 때문).
- 원 실행 폴더 `runs/`(23 GB)는 잡 임시폴더 `/home/wyim/.claude/jobs/d4415b8c/tmp/sim/runs`에만 있고, 잡이 정리되면 사라진다(보존 여부 미결정, TODO 2).

## 서브에이전트·백그라운드 작업 (2026-11-23 확인)
- 실행 중인 것 없음. 이 프로젝트의 Pronghorn 잡 없음(`squeue -u wyim`에는 다른 프로젝트 잡만 있음).
- 구조적 0 구현 작성 에이전트: 완료, 결과는 `e27d383`에 반영됨.
- Codex 시뮬레이션(`codex exec`): 완료, 결과 복사됨(위).
- ccloop 런 `unclassified-meta-p257-20261118`: 2026-11-19 사용자 지시로 stopped.

## 다음 단계
1. 미커밋 문서(REPORT.md v2 절, resume·TODO 등) 커밋 여부와 `em-imputation`의 main 병합·푸시 여부를 사용자에게 확인.
2. TODO 9: 새 EM 표와 이전 EM 표 차이의 원인 확인.
3. TODO 0b·8: robustness 확장(OSD18)과 패널/배경 분할 점검.
