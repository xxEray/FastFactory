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

find_python() {
	for candidate in python3 python python.exe py; do
		if command -v "$candidate" >/dev/null 2>&1 &&
			"$candidate" -c 'import sys; raise SystemExit(0 if sys.version_info[0] >= 3 else 1)' >/dev/null 2>&1; then
			PYTHON="$candidate"
			return 0
		fi
	done
	return 1
}

PYTHON=""
if find_python; then
	echo "Updating README.md from $OUT (using $PYTHON) ..."
	if "$PYTHON" "$ROOT/update.py"; then
		echo "README.md is up to date"
	else
		echo "warning: update.py failed; README.md was not updated" >&2
	fi
else
	echo "warning: python 3 not found (tried python3/python/python.exe/py); skipping README.md update" >&2
fi
