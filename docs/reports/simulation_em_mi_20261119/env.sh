# Source from sim/. Keep package caches, temporary files and plots in this directory.
export PYTHONPATH=/mnt/f/scratch/develop/unclassified_meta-p257fix/src
export PATH=/home/wyim/micromamba/envs/unclassified_meta/bin:$PATH
export UNCLASSIFIED_META_RSCRIPT=/home/wyim/micromamba/envs/unclassified_meta_r/bin/Rscript
export MPLCONFIGDIR="$PWD/.mpl" TMPDIR="$PWD/.tmp" XDG_CACHE_HOME="$PWD/.cache"
export R_USER_CACHE_DIR="$PWD/.cache/R" R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
