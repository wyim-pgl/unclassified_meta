#!/bin/bash
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=emopt_imp
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt/logs/imp_%j.log
# usage: sbatch imp_draws.sbatch <code dir> <sample stem e.g. RNA_S_23_RNA_P262_H_rep2> <M> <seed>
set -euo pipefail
B=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta/bin:$PATH
export PYTHONPATH=$B/$1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $B/impval
python $B/imp_draws.py $B/classes/$2.pkl $B/impval/$1_$2_M$3_s$4.npz $3 $4
