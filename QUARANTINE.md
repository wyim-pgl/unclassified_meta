# QUARANTINE — 폐기된 가정·접근 (재사용 금지)

_최종 갱신: 2026-11-23. 현재 정본은 [resume.md](resume.md)와 `src/unclassified_meta/em.py`._

| 날짜 | 폐기된 것 | 이유 | 현재 근거 |
|---|---|---|---|
| 2026-11-19 | "OSD P257은 식물 한 그루" | 사용자 정정: 식물 3개체(P257/P262/P353) × 조건 × 3 반복 | INTENTION.md, `meta_*.tsv`의 plant 열 |
| 2026-11-19 | "centrifuger-quant EM이 리드 수에도 적용된다" | 소스(1.2.0-r371 Quantifier.hpp) 확인: EM은 abundance 열만, numReads는 균등 분할 | README "Reads that fit several taxa" |
| 2026-11-19 | 고정 theta에서 뽑는 imputation | 모수 불확실성 누락(Codex r1) | em.imputations (Dirichlet tree / SIR 시작 체인) |
| 2026-11-19 | Gibbs 체인만 쓰는 imputation (prior 시작, 60·500 sweep) | 대규모 공유 리드에서 혼합이 매우 느려 사후분포보다 좁고 치우침(Codex r2–r4) | 위와 같음 |
| 2026-11-19 | Bayesian-bootstrap EM imputation | 희소 증거·겹치는 비식별 비율에서 오차(Codex r3) | 위와 같음 |
| 2026-11-19 | EM 정지 규칙: 전역 tol, per-target 상대 tol, Aitken count 외삽, 로그우도 증가 외삽 | 모두 MLE에서 먼 곳에서 거짓 수렴 가능(Codex r1–r8) | 인증된 tangent-plane gap ≤ 0.5 nat (em._gap) |
| 2026-11-19 | SQUAREM 제안을 0으로 clip / 음수면 거부 | clip은 target 영구 제거, 거부는 느린 기어감 | 현재 값의 1e-3 바닥으로 projection |
| 2026-11-19 | dropped·unresolved 후보를 하나의 토큰으로 묶기 | 서로 무관한 생물이 분배에 영향(Codex r1) | 토큰에 taxID 유지 (`#dropped|taxid`) |
| 2026-11-19 | panel_id 앞 12자 FASTA 태그 | 다른 패널 세대가 통과 가능(Codex r4) | 전체 panel_id 태그 |
| 2026-11-19 | "potexvirus-like 4개(CVX/PitVX/ZyVX/SchVX)는 서로 다른 바이러스" | 특징 간 염기 동일성 92–98 %, 샘플 간 상관 0.88–0.99 → 같은 종의 조각·변이체 | RNA/potex_check (Pronghorn OSD18), panel-taxonomy → Opuntia potexvirus A |
| 2026-11-19 | centrifuger 단일 입력 `-U` | 옵션은 `-u` | stageA.sbatch, README |
| 2026-11-20 | EM 기대값 0.5 리드 미만을 구조적 0 판정에서 "없음"으로 보기(구조적 0 문턱) | 사용자 요구: EM 값을 바꾸거나 임의 문턱을 두지 않는다 | ANCOM-BC2 규칙 그대로 + `detection` 근거표 (INTENTION.md) |
| 2026-11-20 | "imputation이 모두 같으면 MI는 완전데이터 t 검정을 재현한다"(df = ν_com) | 사용자 결정으로 Barnard–Rubin 공식을 FMI = 0에도 적용(ν_com(ν_com+1)/(ν_com+3)) | `em.rubin`, `4319f64`, tests/test_mi_pooling.py |
| 2026-11-20 | LCA의 구조적 0 = 그 조건에 생물이 없음 | T. cacti macho 기대 리드 0.35(검출 한계 안), CYVaV는 같은 umbra-like 무리가 P257 healthy에 있음 | docs/reports/osd18_structural_zero/REPORT.md |
| 2026-11-20 | OSD18 첫 배치의 `da_*_em/mi_results.tsv`(round-10 코드) | MI 결합 버그(p=1 fit 재검정) 이전 결과: RNA 종 q<0.05 242개 | `OSD18/<T>/repool_fix/`, v2 `da_*_em_sz/` |
| 2026-11-20 | OSD18 첫 배치 LCA 표(`taxa_*_lca`, `--panel-taxonomy` 없음)로 EM과 비교 | EM 표와 패널 collapse 조건이 달라 비교에 섞임 | `taxa_*_lca_ptx`, v2 결과 |
| 2026-11-20 | "새 EM 코드(main 속도 개선)는 이전 EM과 같은 표를 낸다" | OSD18에서 칸당 최대 1.65 리드 차이 관측 | ISSUE 3b |
| 2026-11-19 | "Validated by three independent statistical methods" 문구 | 세 방법은 같은 데이터·같은 조성 가정을 공유해 독립이 아님(사용자 명세) | "Robust across three complementary differential abundance frameworks", `docs/DESIGN_robustness.md` |
| 2026-11-20 | group term을 이름 접두어로 찾기 | group 이름으로 시작하는 다른 공변량의 term도 잡음 | `robustness.is_group_term`(model-matrix 열 이름 매칭, 7c6ae6d) |
| 2026-11-20 | alpha와 무관한 mde_80 | `--alpha`를 바꿔도 최소 검출 효과가 그대로 | `robustness.z_mde80(alpha)` (7c6ae6d) |
| 2026-11-25 | host genome 경로 `/data/gpfs/assoc/pgl/Eli/OSD/kraken2/input/Ococh_genome/Ococh_noHighGC.fa` | Eli가 `Eli/OSD/input/`으로 옮김(sha256 같음); robustness 잡 6328729–31·6328734–36이 이 때문에 실패 | `panel-build --host ococh=/data/gpfs/assoc/pgl/Eli/OSD/input/Ococh_noHighGC.fa` (`OSD18/rob_res.sbatch`) |
