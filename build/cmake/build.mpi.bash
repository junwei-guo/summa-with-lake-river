#!/usr/bin/env bash
# Build MPI SUMMA + coupled mizuRoute without deleting existing build trees.

export FC=/opt/homebrew/bin/gfortran
export LIBRARY_LINKS=-llapack 

set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export FC="${FC:-$(command -v gfortran)}"
options=(-DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_SUNDIALS=OFF
         "-DCMAKE_BUILD_TYPE=${BUILD_TYPE:-Release}")
if [[ -n "${LIBRARY_LINKS:-}" ]]; then
  options+=(-DSPECIFY_LAPACK_LINKS=ON)
else
  options+=(-DSPECIFY_LAPACK_LINKS=OFF)
fi
cmake -S "$script_dir/.." -B "$script_dir/../cmake_mpi" "${options[@]}"
cmake --build "$script_dir/../cmake_mpi" -j "${BUILD_JOBS:-4}"
