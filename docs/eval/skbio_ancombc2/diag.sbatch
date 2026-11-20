#!/bin/bash
#SBATCH --job-name=skbio_diag
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --time=13-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/skbio_eval/logs/%x_%j.log
# Root-cause diagnostics: initial (pre-bias-correction) estimates and bias EM of R and scikit-bio on the OSD18 RNA
# genus expected-count table; R fit on the subset of features with at most two empty condition x plant cells.
set -uo pipefail
umask 022
E=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/skbio_eval
W=/data/gpfs/assoc/pgl/develop/unclassified_meta_test/OSD18
RS=/data/gpfs/assoc/pgl/bin/conda/conda_envs/unclassified_meta_r/bin/Rscript
SIF=$HOME/scratch/singularity/micromamba.sif
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=$E/diag_$SLURM_JOB_ID; mkdir -p "$O/sub"
C=$W/RNA/taxa_G_em/merged/counts.tsv
singularity exec -B /data/gpfs "$SIF" "$E/venv/bin/python" "$E/init_py.py" "$C" "$W/meta_RNA.tsv" "$O/py"
/usr/bin/time -f "init_r %e s %M KB" "$RS" "$E/init_r.R" "$C" "$W/meta_RNA.tsv" "$O/r"
cat "$O/r_bias.tsv"
singularity exec -B /data/gpfs "$SIF" "$E/venv/bin/python" "$E/subsets.py" "$C" "$W/meta_RNA.tsv" "$O/sub"
for s in le2empty; do
  mkdir -p "$O/sub/r_$s"
  "$RS" "$E/ancombc2.R" --counts "$O/sub/sub_$s.tsv" --features "$W/RNA/taxa_G_em/merged/features.tsv" \
    --metadata "$W/meta_RNA.tsv" --group condition --reference-level healthy --formula "condition + plant" \
    --prv-cut 0 --outdir "$O/sub/r_$s" > "$O/sub/r_$s/log.txt" 2>&1; echo "R $s exit=$?"
  singularity exec -B /data/gpfs "$SIF" "$E/venv/bin/python" -I "$E/skbio_fit.py" "$O/sub/sub_$s.tsv" \
    "$W/meta_RNA.tsv" "$O/sub/py_$s.tsv" --prv-cut 0 --r-na; echo "Py $s exit=$?"
done
echo "done $(date)"
