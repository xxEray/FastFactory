#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="judge-cpp:13"
OUT="$ROOT/fingerprint.txt"
TMP_OUT="$ROOT/.fingerprint.stdout"
TMP_TIME="$ROOT/.fingerprint.time"

USE_DOCKER=0
if [ ! -f /.dockerenv ] && command -v docker >/dev/null 2>&1; then
	USE_DOCKER=1
fi

run_in_env() {
	if [ "$USE_DOCKER" = 1 ]; then
		docker run --rm \
			--memory=3g \
			--memory-swap=3g \
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

for k in 4 5 6 7 8; do
	echo "Running k=$k ..."
	run_in_env bash -c "/usr/bin/time -f 'Time: %es\nMemory: %M' ./solution $((10 ** k + 5000)) > .fingerprint.stdout 2> .fingerprint.time"
	{
		printf 'k=%s, n=%s\n' "$k" "$((10 ** k + 5000))"
		cat "$TMP_OUT"
		awk '{ if ($1 == "Memory:") printf "Memory: %.2f MB\n", $2 / 1024; else print }' "$TMP_TIME"
		printf '\n'
	} >> "$OUT"
done

echo "Wrote $OUT"
