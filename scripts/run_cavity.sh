#!/usr/bin/env bash
# Run the lid-driven cavity cases and report which ones reached steady state.
#
# Usage:  scripts/run_cavity.sh [config.json ...]
#         (no arguments: configs/cavity_re100.json, cavity_re400.json, cavity_re1000.json)
# Env:    CAVITY=path/to/cavity   (default: build/Release/cavity)
#
# Each run's console output goes to out/logs/<config name>.log.
# Exit code: 0 if every case reached steady state, 1 otherwise.

set -u

# Make config arguments absolute, then work from the repository root
configs=()
for arg in "$@"; do
    configs+=("$(realpath "$arg")")
done
cd "$(dirname "$0")/.." || exit 1

if [[ ${#configs[@]} -eq 0 ]]; then
    configs=(configs/cavity_re100.json configs/cavity_re400.json configs/cavity_re1000.json)
fi

cavity="${CAVITY:-build/Release/cavity}"
if [[ ! -x "$cavity" ]]; then
    echo "error: '$cavity' not found; build it first:" >&2
    echo "  cmake --build build/Release -j4" >&2
    exit 1
fi

logDir=out/logs
mkdir -p "$logDir"

rows=()
failed=0
for cfg in "${configs[@]}"; do
    name=$(basename "$cfg" .json)
    log="$logDir/$name.log"
    echo "running $name  (log: $log)"

    start=$SECONDS
    "$cavity" "$cfg" > "$log" 2>&1
    status=$?                               # capture before any other command overwrites $?
    elapsed=$((SECONDS - start))

    case $status in
        0) result="steady" ;;
        2) result="NOT CONVERGED" ;;
        *) result="ERROR" ;;
    esac
    if [[ $status -ne 0 ]]; then
        failed=$((failed + 1))
    fi

    # The app's last status line: "Steady state at ...", "Not converged ..." or "error: ..."
    detail=$(grep -E '^(Steady state|Not converged|error|usage)' "$log" | tail -n 1)
    rows+=("$(printf '%-22s %-14s %4d %8d   %s' "$name" "$result" "$status" "$elapsed" "$detail")")
done

echo
printf '%-22s %-14s %4s %8s   %s\n' "case" "result" "exit" "seconds" "detail"
printf '%s\n' "${rows[@]}"
echo

if [[ $failed -gt 0 ]]; then
    echo "$failed of ${#configs[@]} case(s) failed"
    exit 1
fi
echo "all ${#configs[@]} case(s) reached steady state"
