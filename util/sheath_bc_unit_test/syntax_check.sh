#!/bin/bash
# Syntax check (V0b) for the files carrying the charge-conserving sheath boundary condition,
# without a full JOREK build. Compiles them against the stub modules in this directory with
# -fsyntax-only, for the four combinations of with_TiTe and with_vpar. Catches undeclared
# variables, name collisions (Fortran is case insensitive!) and expression typos in seconds
# instead of a cluster build cycle.
#
#   ./util/sheath_bc_unit_test/syntax_check.sh
set -e
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
FC=${FC:-gfortran}
sed 's|#include "version.h"||' "$root/models/constants.f90" > "$tmp/constants.f90"
cp "$here/stub_mpi_mod.f90" "$tmp/mpi_mod.f90"

files="$root/models/model600/mod_sheath_bc.f90 $root/models/model600/mod_sheath_diag.f90 \
       $root/models/model600/mod_boundary_matrix_open.f90"

status=0
for cfg in "single_T,with_vpar:0,0:.false.:.true." \
           "TiTe,with_vpar:9,10:.true.:.true."     \
           "single_T,no_vpar:0,0:.false.:.false."  \
           "TiTe,no_vpar:9,10:.true.:.false."; do
  name=${cfg%%:*}; rest=${cfg#*:}; idx=${rest%%:*}; rest=${rest#*:}
  titel=${rest%%:*}; vpar=${rest#*:}; ti=${idx%,*}; te=${idx#*,}
  sed -e "s/var_Ti=0, var_Te=0/var_Ti=$ti, var_Te=$te/" -e "s/n_var=8/n_var=10/" \
      -e "s/with_TiTe=\.false\./with_TiTe=$titel/" -e "s/with_vpar=\.true\./with_vpar=$vpar/" \
      "$here/stub_modules.f90" > "$tmp/stubs.f90"
  rm -f "$tmp"/*.mod
  err=$($FC -fsyntax-only -ffree-line-length-none -J"$tmp" \
        "$tmp/constants.f90" "$tmp/stubs.f90" "$tmp/mpi_mod.f90" $files 2>&1 | grep -i "^Error\|: Error" || true)
  if [ -z "$err" ]; then
    echo "PASS  $name"
  else
    echo "FAIL  $name"; echo "$err"; status=1
  fi
done
exit $status
