# Plotting model600 boundary velocities

For the Mach1 cutoff A/B comparison, use the same input with only
`mach1_grazing_release` changed:

```fortran
mach1_grazing_release = .true.  ! release at |B.n/|B|| <= sin(min_sheath_angle)
min_sheath_angle = 1.d0        ! degrees
```

The default `.false.` retains Mach enforcement at nonzero incidence. Both runs
use the restored signed, unclipped ExB correction; the flag controls release,
not restoration of the former SOLPS clip. Exact tangency is released in both
cases to avoid division by zero.

Use the input namelist matching each run and postprocess each run in its own
directory. For example, in `jorek2_postproc`:

```text
namelist input
si-units
set nsub_bnd 8
for step 100 do
  expressions length bnd_contour bnd_type R Z bn_unit vexb_norm vpar_norm vtot_norm vpar_speed
  boundary_quantities
done
```

Change the step, or use `for step 20 to 100 do` for a sequence. The output is
`boundary_quantities*.dat`, with the selected expressions as columns.

| Expression | Meaning with `si-units` |
|---|---|
| `length` | Curved-boundary arc length in metres, starting at zero for each contour |
| `bnd_contour` | Connected contour number, starting at one |
| `bnd_type` | Boundary node type; `-1` at interpolated points on mixed-type edges |
| `bn_unit` | Signed magnetic incidence `B.n/|B|`, dimensionless |
| `vexb_norm` | Normal component of the model600 u-driven ExB velocity, m/s |
| `vpar_norm` | Normal component of the parallel flow `Vpar B.n`, m/s |
| `vtot_norm` | Sum of those normal velocities, m/s |
| `vpar_speed` | Signed parallel flow speed `Vpar |B|`, m/s; positive along B |

The three normal velocities are positive **outward from the plasma domain**.
They obey `vtot_norm = vexb_norm + vpar_norm` in model600, including inward
ExB flow and reversed parallel flow. `vexb_norm` uses
`R*(-u_Z*n_R + u_R*n_Z)` directly, without the Mach BC's correction or any
clipping. These are the actual reconstructed velocities, not prescribed targets.

`vu_norm` is the existing equivalent u-velocity diagnostic in model600.
**`ExB_norm` is an electromagnetic energy flux, not a velocity.** The existing
`vpar` expression gives a speed in SI units but the evolved `Vpar` coefficient in
JOREK units. `vpar_speed` includes `|B|` in both unit systems.

To locate the grazing band, compare `abs(bn_unit)` against
`sin(min_sheath_angle*pi/180)`. Equality belongs to the released side when the
release flag is enabled. This sampled incidence shows geometry, not an exact
record of constrained matrix rows: nodal Mach conditions use nodal geometry,
quadrature enforcement uses quadrature geometry, and a shared corner may still
be constrained by its other incident edge. In 3D the sampled toroidal plane may
also differ from the axisymmetric geometry used by the nodal BC.

Boundary points follow the connected boundary-element list and each edge's
stored traversal direction. `length` integrates the curved geometry using
eight-point Gauss quadrature on each sampling interval, including the last
interval to an edge's endpoint. It resets at each disconnected contour; plot
each `bnd_contour` separately rather than connecting distinct walls. The closing
point is not duplicated. Increasing `nsub_bnd` interpolates the existing
solution more finely; it does not increase simulation resolution.

At every sampled vertex, `bnd_type` is the actual stored node boundary number.
At interior samples, it is the common endpoint type, or `-1` if the endpoints
differ. Add `bnd_type1 bnd_type2` to the expressions to see both endpoint types
in traversal order. With `set nsub_bnd 1`, all samples are vertices.

Local regression checks:

```sh
python3 tests/postproc/test_boundary_position.py
python3 tests/postproc/test_boundary_velocity.py
```

These compile the production sampling routine and expression formulas against
analytic fixtures. They do not replace a full postprocessor build or comparison
against simulation restart data.
