#!/bin/bash
# Serial harness for the floating_u wall BCs: compiles the production model600 boundary assembler, mod_floating_u and
# mod_wall_diag against stub modules (fixtures.f90), runs the FD checks, the diagnostics smoke test and the develop
# identity check (needs git; DEV_REF = the develop routine to compare with, default the merge-base with develop).
set -e
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../.." && pwd)
work=${TMPDIR:-/tmp}/floating_u_nodal_harness; rm -rf "$work"; mkdir -p "$work/old"; cd "$work"
F="-fopenmp -ffree-line-length-none -Wno-unused-variable -Wno-unused-dummy-argument -Wno-maybe-uninitialized"
M=$root/models/model600
gfortran -c $F "$here/fixtures.f90"
gfortran -c $F -Wall $M/mod_floating_u.f90 $M/mod_wall_diag.f90 $M/mod_boundary_matrix_open.f90
objs="fixtures.o mod_floating_u.o mod_wall_diag.o mod_boundary_matrix_open.o"
gfortran $F -o test_nodal_flux $objs "$here/test_nodal_flux.f90" && ./test_nodal_flux
gfortran $F -o test_wall_diag $objs "$here/test_wall_diag.f90" && ./test_wall_diag | head -4
gfortran $F -o dump_new $objs "$here/dump_edge.f90" && ./dump_new && mv dump.bin dump_new.bin
ref=${DEV_REF:-$(git -C "$root" merge-base develop HEAD)}
( cd old && git -C "$root" show "$ref":models/model600/mod_boundary_matrix_open.f90 > bmo_old.f90 && cp ../*.mod . \
  && gfortran -c $F bmo_old.f90 && gfortran $F -o dump_old ../fixtures.o bmo_old.o "$here/dump_edge.f90" && ./dump_old )
python3 - <<PY
import numpy as np
o=np.fromfile('old/dump.bin','<f8'); w=np.fromfile('dump_new.bin','<f8'); d=np.abs(o-w); sc=np.maximum(np.abs(o),np.abs(w))
rel=(d/np.maximum(sc,1e-300)).max(); same=((o==0)==(w==0)).all()
print('develop identity (floating_u off): max rel diff %.1e, sparsity equal: %s' % (rel, same))
assert same and rel < 1e-12, 'FAIL: differs from develop with floating_u off'
PY
echo ALL PASS
