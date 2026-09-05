"""Optional integration syntax check: requires gfortran and fparser (tested 0.2.5).

This does not resolve module interfaces or replace a full MPI production build.
"""
import pathlib
import subprocess

from fparser.common.readfortran import FortranStringReader
from fparser.two.parser import ParserFactory

root = pathlib.Path(__file__).resolve().parents[2]
parse = ParserFactory().create(std="f2008")
files = [
    "models/model600/mod_boundary_conditions.f90",
    "models/model600/mod_boundary_matrix_open.f90",
    "models/model600/initialise_parameters.f90",
    "models/model600/mod_sheath_diag.f90",
    "matrix/construct_matrix_mod.f90",
]
for filename in files:
    source = subprocess.check_output(
        ["gfortran", "-E", "-cpp", "-DJOREK_MODEL=600", "-Itools",
         "-ffree-line-length-none", filename],
        cwd=root, text=True,
    )
    parse(FortranStringReader(source))
    print("PARSE PASS:", filename)
