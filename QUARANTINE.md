# QUARANTINE — 폐기된 가정·접근 (재사용 금지)

_최종 갱신: 2026-11-19. 현재 정본은 [resume.md](resume.md)와 `src/unclassified_meta/em.py`._

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
