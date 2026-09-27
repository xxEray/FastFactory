#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="judge-cpp:13"

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

echo "Building verify.exe ..."
run_in_env g++ src/verify.cpp -o verify.exe -Ofast -fopenmp

for k in 4 5 6 7 8; do
	echo "Running verify.exe $k ..."
	run_in_env ./verify.exe "$k"
done
