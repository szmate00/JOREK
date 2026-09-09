"""Check production diagnostic expressions against an analytic model600 velocity.

Run: python3 tests/postproc/test_boundary_velocity.py
Checks outward-normal signs, both unit systems and normal-flow decomposition.
"""
from pathlib import Path
import os
import re
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "diagnostics/new_diag/mod_expression.f90").read_text()
names = ["bn_unit", "vexb_norm", "vpar_norm", "vtot_norm", "vpar_speed"]
assignments = []
for index, name in enumerate(names, 1):
    block = re.search(
        r"case\s*\(\s*'" + name + r"'\s*\)(.*?)(?=\n\s*case\s*\()",
        source, re.S,
    ).group(1)
    expression = re.search(r"^\s*res\s*=\s*(.+)", block, re.M).group(1)
    assignments.append(f"out({index}) = {expression}")

program = """
program check
implicit none
real*8 :: R,u0_Z,u0_R,nmlR,nmlZ,vpar0,Bnorm,Btot,fact_time,VR,VZ,out(5),expected(5)
integer :: k, idir
R=2d0; u0_Z=4d0; u0_R=3d0; vpar0=2d0
Btot=sqrt(3d0**2+2.5d0**2+7d0**2)
VR=-R*u0_Z+vpar0*3d0; VZ=R*u0_R-vpar0*2.5d0
do idir=-1,1,2
  nmlR=idir*0.8d0; nmlZ=idir*0.6d0
  Bnorm=3d0*nmlR-2.5d0*nmlZ
  do k=1,2
    fact_time=1d0
    if(k==2) fact_time=3.7d-7
""" + "\n".join(assignments) + """
    expected=(/idir*0.9d0/Btot, -idir*2.8d0/fact_time, idir*1.8d0/fact_time, &
               -idir*1d0/fact_time, 2d0*Btot/fact_time/)
    if(any(abs(out-expected)>1d-12*max(1d0,abs(expected)))) error stop 'velocity/sign/unit mismatch'
    if(abs(out(4)-out(2)-out(3))>1d-12*max(1d0,abs(out(4)))) error stop 'normal velocity sum'
  enddo
enddo
print *, 'PASS: production velocity expressions, normals, units and normal-flow sum'
end program
"""
with tempfile.TemporaryDirectory(prefix="jorek-bnd-velocity-") as directory:
    path = Path(directory)
    (path / "test.f90").write_text(program)
    subprocess.run(
        [os.environ.get("FC", "gfortran"), "-ffree-line-length-none", "-fcheck=all",
         "-ffpe-trap=invalid,zero,overflow", "test.f90", "-o", "test"],
        cwd=path, check=True,
    )
    subprocess.run([str(path / "test")], cwd=path, check=True)
