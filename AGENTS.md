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
