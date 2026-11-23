# AGENTS.md — unclassified_meta

## Read first
- `INTENTION.md` (goal and decisions), `resume.md` (current state), `TODO.md`, `ISSUE.md`, `QUARANTINE.md` (settled mistakes; do not reuse).
- User-facing documentation: `README.md` (sections 10, "Reads that fit several taxa", "One organism, one feature").

## Rules
- Authorship: Won Yim <wyim@unr.edu> only; no AI co-author trailers. Commit or push only when the user asks.
- The editable install in the `unclassified_meta` micromamba env points at the main checkout `/mnt/f/scratch/develop/unclassified_meta`. To test this worktree run with `PYTHONPATH=src`, `PATH=<env>/bin:$PATH` and `UNCLASSIFIED_META_RSCRIPT=<unclassified_meta_r env>/bin/Rscript`; otherwise tests skip or test the wrong code.
- Every producer seals its outputs (`MANIFEST.json` / `provenance.json`) last; every consumer verifies before reading. Keep that pattern for new outputs; remove stale files only through `manifest.remove_stale`.
- Defaults must stay backward compatible: `--assign lca`, `--ambiguous ledger`.
- Centrifuger single-end input is `-u` (lowercase). On Pronghorn cancel jobs by ID (`scancel -n a,b` with a comma list cancels nothing).
- Slurm: partitions `cpu-s1-pgl-0` or `cpu-s2-core-0`; 13-day limit; never `cpu-s3-sponsored-0`.
- Test command (from the worktree root): `E=/home/wyim/micromamba/envs; PYTHONPATH=src PATH=$E/unclassified_meta/bin:$PATH UNCLASSIFIED_META_RSCRIPT=$E/unclassified_meta_r/bin/Rscript python -m pytest -q -p no:warnings tests` (about 35-40 min; run in the background).
- DIAMOND 2.1.10 `blastx` can deadlock in `tests/test_panel_e2e.py` (futex wait, 0 % CPU; upstream bbuchfink/diamond #943, #849). Kill only that `diamond` process and rerun the affected tests.
- Structural zeros: keep ANCOM-BC2's own rule (count != 0); never threshold or alter EM values; never Rubin-pool structural-zero indicators or fabricate LFC/SE/p for them (see INTENTION.md).
- Real-data runs on Pronghorn use a code copy, not the checkout: `unclassified_meta_test/em_code_v2/` (with a `COMMIT` file) for OSD18 v2; sync with rsync from a clean committed tree.
