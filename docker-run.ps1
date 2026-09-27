param(
    [int]$k,
    [string]$Argument = ""
)


docker run --rm `
  --memory=3g `
  --memory-swap=3g `
  --cpuset-cpus=0-7 `
  --network=none `
  -e OMP_NUM_THREADS=8 `
  -v "${PWD}:/work:rw" `
  -w /work `
  judge-cpp `
  g++ src/solution.cpp -o solution.exe -Ofast -fopenmp $Argument

echo Built

docker run --rm `
  --memory=3g `
  --memory-swap=3g `
  --cpuset-cpus=0-7 `
  --network=none `
  -e OMP_NUM_THREADS=8 `
  -v "${PWD}:/work:rw" `
  -w /work `
  judge-cpp `
  time timeout 600 ./solution.exe $k