"""Optional fparser integration check; does not replace a full MPI build."""
import pathlib
import subprocess
import re
from fparser.common.readfortran import FortranStringReader
from fparser.two.parser import ParserFactory

root = pathlib.Path(__file__).resolve().parents[2]
parse = ParserFactory().create(std="f2008")
files = [
    "models/model600/mod_boundary_conditions.f90",
    "models/model600/mod_boundary_matrix_open.f90",
    "models/model600/mod_elt_matrix_fft.f90",
    "models/model600/initialise_parameters.f90",
    "models/model600/mod_floating_transport.f90",
    "models/model600/mod_floating_transport_diag.f90",
    "models/model600/mod_floating_boundary_edges.f90",
    "matrix/construct_matrix_mod.f90",
    "core/mod_jorek_timestepping.f90",
]
for filename in files:
    source = subprocess.check_output(
        ["gfortran", "-E", "-cpp", "-DJOREK_MODEL=600", "-Itools", "-Imodels",
         "-ffree-line-length-none", filename], cwd=root, text=True)
    # Existing volume code contains binary-plus followed by unary-plus, a GNU
    # extension. Normalise this notation ONLY for the standards-based parser.
    source = re.sub(r"([+-])\s*&\s*\n(?:[ \t]*!.*\n|[ \t]*\n)*[ \t]*([+-])",
                    lambda m: ("+" if m[1] == m[2] else "-") + " &\n ", source)
    parse(FortranStringReader(source))
    print("PARSE PASS:", filename)
