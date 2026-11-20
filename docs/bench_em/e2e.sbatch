#!/bin/bash
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=emopt_e2e
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt/logs/e2e_%x_%j.log
# usage: sbatch e2e.sbatch <code dir: baseline_src|opt3_src> <set: small|medium|complete> [rank S|G]
# The production taxa-table command of OSD18 stage A (RNA), on a subset of samples for small/medium.
set -euo pipefail
umask 022
B=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/em_opt
W=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18
TX=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/osd_taxonomy
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta/bin:$PATH
export PYTHONPATH=$B/$1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R=${3:-S}
O=$B/e2e/$1_$2_$R
rm -rf $O
hostname; lscpu | grep "Model name"; python -c "import unclassified_meta, numpy; print('CODE', unclassified_meta.__file__, 'numpy', numpy.__version__)"
EXTRA=""
if [ $2 = complete ]; then
  REP=$W/reports_RNA.tsv
  EXTRA="--merge-panel $W/RNA/quant_em/counts.tsv --panel-taxonomy $W/RNA/ptx/panel_taxonomy.tsv"
else
  REP=$B/reports_RNA_$2.tsv
fi
/usr/bin/time -f "TIME $1 $2 $R elapsed=%e s user=%U sys=%S cpu=%P maxrss_kb=%M" \
  python -m unclassified_meta.cli taxa-table --reports $REP --per-read --nodes $TX/cfr_nr.nodes.dmp \
  --names $TX/cfr_nr.names.tsv --taxonomy-name cfr_nr_20260908 --min-score 64 --rank $R \
  --drop-lineage Viridiplantae,Metazoa --assign em --imputations 10 --seed 11 $EXTRA -o $O
