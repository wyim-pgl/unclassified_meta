#!/bin/bash
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --job-name=sz_diag
set -euo pipefail
W=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18; D=$W/szdiag; TX=$(grep -o "TX=[^ ;]*" $W/stageA.sbatch | head -1 | cut -d= -f2)
export PATH=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta/bin:$PATH
RS=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta_r/bin/Rscript
for R in S G; do $RS $D/struc_zero.R $W/RNA $R $D/struc_zero_RNA_$R.tsv; done
python $D/reads_diag.py $TX/cfr_nr.nodes.dmp $TX/cfr_nr.names.tsv $W/reports_RNA.tsv $D/reads_RNA_S.tsv
echo DONE
