#!/bin/bash
# Serial harness for the sheath current row (bcs%sheath_j): compiles the PRODUCTION model600 boundary assembler and
# mod_floating_u against stub modules (fixtures.f90), runs the FD checks of the zj rows, and compares the assembled
# edge with sheath_j off against develop's routine (needs git; DEV_REF = the develop commit, default the merge-base).
# mod_boundary_conditions.f90 (the DOF release) needs the full build and is not compiled here.
set -e
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../.." && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/sheath_j_harness.XXXXXX"); mkdir -p "$work/old"; cd "$work"
F="-ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -Wno-unused-variable -Wno-unused-dummy-argument -Wno-maybe-uninitialized"
M=$root/models/model600
gfortran -c $F "$here/fixtures.f90"
gfortran -c $F -Wall $M/mod_floating_u.f90 $M/mod_boundary_matrix_open.f90
objs="fixtures.o mod_floating_u.o mod_boundary_matrix_open.o"
gfortran $F -o test_sheath_j $objs "$here/test_sheath_j.f90"
./test_sheath_j
gfortran $F -o dump_new $objs "$here/dump_edge.f90"
./dump_new
mv dump.bin dump_new.bin
ref=${DEV_REF:-$(git -C "$root" merge-base origin/develop HEAD)}
( cd old && git -C "$root" show "$ref":models/model600/mod_boundary_matrix_open.f90 > bmo_old.f90 && cp ../*.mod . \
  && gfortran -c $F bmo_old.f90 && gfortran $F -o dump_old ../fixtures.o bmo_old.o "$here/dump_edge.f90" && ./dump_old )
cmp -s old/dump.bin dump_new.bin && echo ' PASS: sheath_j off: edge matrix and residual bitwise identical to develop' \
  || { echo ' FAIL: sheath_j off differs from develop'; exit 1; }
echo "ALL PASS (build in $work)"
