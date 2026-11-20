#!/bin/bash
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=emopt_dump
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt/logs/dump_%j.log
set -euo pipefail
B=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta/bin:$PATH
export PYTHONPATH=$B/baseline_src
/usr/bin/time -f "TIME dump $1 $2 elapsed=%E maxrss_kb=%M" python $B/dump_classes.py $1 $2 $B/classes
