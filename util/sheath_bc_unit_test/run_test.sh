#!/bin/bash
# Unit test (V0) for the sheath j-V characteristic used by the electric-potential
# boundary condition in model600. Compiles the REAL models/model600/mod_sheath_bc.f90
# against a stub phys_module, so it needs nothing but a Fortran compiler.
#
#   ./util/sheath_bc_unit_test/run_test.sh
#
# Checks: Lambda(Ti/Te), j = 0 exactly at the floating potential, Phi_float = Lambda*Te/e,
# j -> j_sat and dj/du -> 0 in ion saturation, electron saturation at Phi = 0,
# all analytic derivatives against central differences in three regimes, and that the
# limiter and the characteristic stay C1, monotone and finite for X in [-200,200].
set -e
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
FC=${FC:-gfortran}
sed 's|#include "version.h"||' "$root/models/constants.f90" > "$tmp/constants.f90"
$FC -O1 -ffree-line-length-none -J"$tmp" -o "$tmp/test_sheath_bc" \
    "$tmp/constants.f90" \
    "$here/stub_phys_module.f90" \
    "$root/models/model600/mod_sheath_bc.f90" \
    "$here/test_sheath_bc.f90"
"$tmp/test_sheath_bc"
