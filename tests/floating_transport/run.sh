#!/usr/bin/env bash
# Serial (no MPI) checks for the floating-potential boundary condition and the weak Bohm row.
# Compiles the PRODUCTION model600 boundary assembler against small fixture modules.
set -euo pipefail
repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/jorek-floating-tests.XXXXXX")
fc=${FC:-gfortran}
flags="-ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -Wall -Wno-compare-reals -Wno-unused-variable -Wno-unused-dummy-argument"
cd "$build_dir"

python3 "$repo_dir/tests/floating_transport/check_imports.py"

$fc $flags \
  "$repo_dir/tests/floating_transport/fixtures.f90" \
  "$repo_dir/models/model600/mod_floating_u.f90" \
  "$repo_dir/tests/floating_transport/test_floating_u.f90" -o test_floating_u
./test_floating_u

$fc $flags \
  "$repo_dir/tests/floating_transport/fixtures.f90" \
  "$repo_dir/models/model600/mod_boundary_edges.f90" \
  "$repo_dir/tests/floating_transport/test_boundary_edges.f90" -o test_boundary_edges
./test_boundary_edges

$fc $flags \
  "$repo_dir/tests/floating_transport/fixtures.f90" \
  "$repo_dir/models/model600/mod_boundary_matrix_open.f90" \
  "$repo_dir/tests/floating_transport/test_weak_mach.f90" -o test_weak_mach
./test_weak_mach

echo "ALL TESTS PASSED (build retained in $build_dir)"
