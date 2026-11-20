#!/bin/bash
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=emopt_run
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt/logs/run_%x_%A_%a.log
# usage: sbatch --array=1-18 run_em.sbatch <code dir: baseline_src|opt_src> <T_R prefix e.g. RNA_S> <tag>
set -euo pipefail
B=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta/bin:$PATH
export PYTHONPATH=$B/$1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
f=$(ls $B/classes/$2_*.pkl | sort -V | sed -n "${SLURM_ARRAY_TASK_ID}p")   # reports-file order
mkdir -p $B/results/$3
/usr/bin/time -f "TIME elapsed=%E user=%U sys=%S cpu=%P maxrss_kb=%M" \
  python $B/run_em.py $f $B/results/$3/$(basename $f) 10 11 $((SLURM_ARRAY_TASK_ID - 1))
