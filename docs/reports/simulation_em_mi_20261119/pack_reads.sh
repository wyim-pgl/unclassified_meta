#!/usr/bin/env bash
# Archive only the raw hit listings produced by this simulation and validate every byte.
set -euo pipefail
cd -- "$(dirname -- "$0")"
test ! -e per_read_inputs.tar.zst
test ! -e per_read_inputs.partial.tar.zst
shopt -s nullglob
dirs=(runs/alt_[0-9][0-9][0-9] runs/null_[0-9][0-9][0-9])
((${#dirs[@]} > 0))
inputs=()
for d in "${dirs[@]}"; do
    test -f "$d/singleton.tsv"
    for i in {00..11}; do
        file="$d/s$i.cfr.tsv"
        test -s "$file"
        inputs+=("$file")
    done
done
tar -I 'zstd -1 -T1' -cf per_read_inputs.partial.tar.zst "${inputs[@]}"
tar -I zstd -df per_read_inputs.partial.tar.zst "${inputs[@]}"
mv per_read_inputs.partial.tar.zst per_read_inputs.tar.zst
# Every member matched its live input. The archive restores their original paths.
for d in "${dirs[@]}"; do
    while [[ ! -s "$d/run.json" || ( -d "$d/pinned_fit" && ! -s "$d/pinned_fit/run.json" ) ]]; do
        sleep 5
    done
done
rm -- "${inputs[@]}"
