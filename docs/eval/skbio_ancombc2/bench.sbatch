#!/bin/bash
#SBATCH --job-name=skbio_vs_r
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/skbio_eval/logs/%x_%j.log
# R ANCOM-BC2 (project r/ancombc2.R, ANCOMBC 2.14.0) vs scikit-bio 0.7.4 ancombc2 on identical OSD18 RNA matrices.
# Fits run one at a time, single-threaded, on the same node. Each fit includes the pseudo-count sensitivity
# analysis (4 core fits), as in production.
set -uo pipefail
umask 022
E=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/skbio_eval
W=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18
RS=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta_r/bin/Rscript
RSCRIPT=$E/ancombc2.R                    # copy of src/unclassified_meta/r/ancombc2.R at the evaluated commit
SIF=$HOME/scratch/singularity/micromamba.sif
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=$E/bench_$SLURM_JOB_ID; mkdir -p "$O"
echo "node $(hostname) cpus $SLURM_CPUS_PER_TASK $(date)"
lscpu | grep -E "Model name|^CPU\(s\)" ; free -g | head -2

counts_of() {   # level matrix -> path
  if [ "$2" = exp ]; then echo "$W/RNA/taxa_$1_em/merged/counts.tsv"; else echo "$W/RNA/taxa_$1_em/merged/imputations/counts_$2.tsv"; fi
}
run_r() {
  local lv=$1 m=$2 d=$O/r_$1_$2; mkdir -p "$d"
  /usr/bin/time -v -o "$d/time.txt" "$RS" "$RSCRIPT" --counts "$(counts_of $lv $m)" \
    --features "$W/RNA/taxa_${lv}_em/merged/features.tsv" --metadata "$W/meta_RNA.tsv" --group condition \
    --reference-level healthy --formula "condition + plant" --threads 1 --outdir "$d" > "$d/log.txt" 2>&1
  echo "R  $lv $m exit=$? $(grep -E 'Elapsed|Maximum resident' "$d/time.txt" | tr -s ' ' | tr '\n' ' ')"
}
run_py() {
  local lv=$1 m=$2 rep=$3 d=$O/py_$1_$2_$3; mkdir -p "$d"
  /usr/bin/time -v -o "$d/time.txt" singularity exec -B /data/gpfs "$SIF" "$E/venv/bin/python" -I "$E/skbio_fit.py" \
    "$(counts_of $lv $m)" "$W/meta_RNA.tsv" "$d/results.tsv" --r-na --timing "$d/timing.json" > "$d/log.txt" 2>&1
  echo "Py $lv $m rep$rep exit=$? $(grep -E 'Elapsed|Maximum resident' "$d/time.txt" | tr -s ' ' | tr '\n' ' ')"
}
for rep in 1 2 3; do run_py G exp $rep; done
for m in 001 002 003; do run_py G $m 1; done
run_py S exp 1
for m in exp 001 002 003; do run_r G $m; done
run_r S exp
echo "done $(date)"
