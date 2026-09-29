#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="judge-cpp:13"
OUT="$ROOT/extreme_test.txt"
TMP_OUT="$ROOT/.extreme_test.stdout"
TMP_TIME="$ROOT/.extreme_test.time"

n=300000000

USE_DOCKER=0
if [ ! -f /.dockerenv ] && command -v docker >/dev/null 2>&1; then
	USE_DOCKER=1
fi

run_in_env() {
	if [ "$USE_DOCKER" = 1 ]; then
		docker run --rm \
			--cpuset-cpus=0-7 \
			--network=none \
			-e OMP_NUM_THREADS=8 \
			-v "$ROOT:/work:rw" \
			-w /work \
			"$IMAGE" "$@"
	else
		( cd "$ROOT" && "$@" )
	fi
}

echo "Building solution ..."
run_in_env g++ src/solution.cpp -o solution -Ofast -fopenmp -march=native

: > "$OUT"
trap 'rm -f "$TMP_OUT" "$TMP_TIME"' EXIT

echo "Running n=$n ..."
run_in_env bash -c "/usr/bin/time -f 'Time: %es\nMemory: %M' ./solution $n > $TMP_OUT 2> $TMP_TIME"
{
	printf 'n=%s\n' "$n"
	cat "$TMP_OUT"
	awk '{ if ($1 == "Memory:") printf "Memory: %.2f MB\n", $2 / 1024; else print }' "$TMP_TIME"
	printf '\n'
} >> "$OUT"

echo "Wrote $OUT"
