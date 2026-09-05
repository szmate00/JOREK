#!/usr/bin/env bash
set -euo pipefail
repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/jorek-floating-tests.XXXXXX")
compiler=${FC:-gfortran}
cd "$build_dir"
"$compiler" -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -fimplicit-none -Wall -Wextra -Wno-compare-reals \
  "$repo_dir/models/model600/mod_floating_transport.f90" \
  "$repo_dir/tests/floating_transport/test_kernels.f90" -o test_kernels
./test_kernels
"$compiler" -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan \
  "$repo_dir/tests/floating_transport/fixtures.f90" \
  "$repo_dir/models/model600/mod_floating_u.f90" \
  "$repo_dir/models/model600/mod_floating_boundary_edges.f90" \
  "$repo_dir/models/model600/mod_floating_transport_diag.f90" \
  "$repo_dir/models/model600/mod_floating_transport.f90" \
  "$repo_dir/models/model600/mod_boundary_matrix_open.f90" \
  "$repo_dir/tests/floating_transport/test_boundary.f90" -o test_boundary
./test_boundary
printf 'Test build retained in %s\n' "$build_dir"
