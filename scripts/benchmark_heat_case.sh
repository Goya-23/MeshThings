#!/usr/bin/env bash
set -euo pipefail

process_counts="${1:-1 2 4}"
case_file="${2:-examples/large_composite_heat_multimaterial.plc}"
nx="${3:-110}"
ny="${4:-90}"
nz="${5:-8}"
method="${6:-i}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${repo_root}/build/heat_conduction_mpi"
log_dir="${repo_root}/logs/benchmarks"
mkdir -p "${log_dir}"

if [[ ! -x "${binary}" ]]; then
    echo "Missing ${binary}; build first with: CXX=g++ CC=gcc cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j" >&2
    exit 1
fi

printf 'processes,wedge_cells,total_seconds_max,rss_max_kib,rss_sum_kib,log_file\n'
for np in ${process_counts}; do
    log_file="${log_dir}/heat_np${np}_nx${nx}_ny${ny}_nz${nz}.log"
    mpirun --allow-run-as-root -np "${np}" "${binary}" "${repo_root}/${case_file}" "${method}" "${nx}" "${ny}" "${nz}" \
        > "${log_file}" 2>&1

    wedge_cells="$(awk '/Wedge mesh:/ {print $(NF-2)}' "${log_file}" | tail -n 1)"
    total_metric="$(awk '/METRIC phase=total/ {line=$0} END {print line}' "${log_file}")"
    total_seconds="$(awk -v line="${total_metric}" 'BEGIN {n=split(line, fields, " "); for (i=1; i<=n; ++i) if (fields[i] ~ /^seconds_max=/) {split(fields[i], value, "="); print value[2]}}')"
    rss_max="$(awk -v line="${total_metric}" 'BEGIN {n=split(line, fields, " "); for (i=1; i<=n; ++i) if (fields[i] ~ /^rss_max_kib=/) {split(fields[i], value, "="); print value[2]}}')"
    rss_sum="$(awk -v line="${total_metric}" 'BEGIN {n=split(line, fields, " "); for (i=1; i<=n; ++i) if (fields[i] ~ /^rss_sum_kib=/) {split(fields[i], value, "="); print value[2]}}')"
    printf '%s,%s,%s,%s,%s,%s\n' "${np}" "${wedge_cells:-unknown}" "${total_seconds:-unknown}" "${rss_max:-unknown}" "${rss_sum:-unknown}" "${log_file}"
done
