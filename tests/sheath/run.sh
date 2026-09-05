#!/usr/bin/env bash
set -euo pipefail
sheath_repo=$(cd "$(dirname "$0")/../.." && pwd)
sheath_build=$(mktemp -d "${TMPDIR:-/tmp}/jorek-sheath-tests.XXXXXX")
cd "$sheath_build"
sheath_fc=${FC:-gfortran}
sheath_flags=(-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none)
"$sheath_fc" "${sheath_flags[@]}" -c "$sheath_repo/tests/sheath/stubs.f90" \
  "$sheath_repo/models/model600/mod_sheath_bc.f90" \
  "$sheath_repo/matrix/mod_locate_irn_jcn.f90" \
  "$sheath_repo/models/mod_assembly.f90" \
  "$sheath_repo/models/model600/mod_sheath_trace.f90" \
  "$sheath_repo/models/model600/mod_sheath_boundary_edges.f90"
for sheath_test in test_wall_law test_trace_edges; do
  "$sheath_fc" "${sheath_flags[@]}" "$sheath_repo/tests/sheath/$sheath_test.f90" ./*.o -o "$sheath_test"
  "./$sheath_test"
done
printf 'Standalone test build: %s\n' "$sheath_build"
