# Sheath current boundary condition (model600)

Imposes the full sheath current-voltage characteristic on the electrostatic potential at a
material wall, replacing the Dirichlet `u` rows on the boundary types you flag.

Physics, derivation and linearisation: Artola, *Sheath boundary conditions for the electric
potential in JOREK* (2026-07-30), eqs. (1)-(18). Implementation and the full reasoning:
`models/model600/mod_sheath_current.f90` header. Nothing is duplicated here.

## Switching it on

```fortran
bcs(1)%sheath_j = .true.
bcs(4)%sheath_j = .true.
bcs(5)%sheath_j = .true.
bcs(9)%sheath_j = .true.

bcs(1)%dirichlet%u = .false.   ! optional: the sheath row claims these rows anyway
sheath_Lambda      = 3.d0      ! default; 2.84 is the deuterium Ti=Te value
```

`sheath_Lambda` is the floating sheath potential drop in units of `Te/e`,
`-0.5*log(2*pi*(me/mi)*(1+Ti/Te))`. It is the only number the condition takes, and it is a
physical constant, not a tuning knob. There is no gate, threshold, clip, slope limiter or
relaxation gain anywhere in this boundary condition.

At the type-3/type-2 corner the existing exception is preserved: the sheath row is not assembled
there and the Dirichlet `u` rows are left in place, so `u` is never left without an equation.

## Reading the log

Two lines per matrix construction, whenever at least one sheath row was assembled:

```
 [sheath_j] rows=   1248 nsat=     3  ePhi/Te-Lam= -2.104E+00  4.551E+00  j/jsat= -8.12E-01  9.98E-01
 [sheath_j] rowscale=  1.207E-03  4.512E+01  max|res|[V]=  3.418E+00  at R,Z=   1.6012  -1.1104  bnd type  4
```

| field | meaning | what is wrong if it moves |
|---|---|---|
| `rows` | sheath rows assembled | a drop means boundary nodes stopped qualifying |
| `nsat` | nodes at or past ion saturation, `j/jsat >= 1` | a few is expected; a whole target is not |
| `ePhi/Te-Lam` | normalised potential measured from floating | the two targets on opposite signs is the thermoelectric current, and is physical |
| `j/jsat` | current relative to ion saturation | sustained values far past 1 mean the plasma is demanding current the sheath cannot pass |
| `rowscale` | the exact factor divided out of each row | grows like `exp(ePhi/Te)`; a spread of decades is the two branches, a runaway is not |
| `max|res|` | worst node's distance from the characteristic, in volts, with its position and boundary type | **this is the number that grows before a crash**, and the boundary type says which node family loses the condition first |

The self-test of the normalisation runs once, on the first step where the condition is active, and
aborts the run if it fails. It checks the sign convention, the floating level in volts, and every
Jacobian column against a finite difference of the characteristic coded independently from eq. (6).

## What has been checked, and what has not

Verified locally, with `gfortran`:

- `mod_sheath_current.f90` compiles clean with `-Wall`, and its self-test passes for both signs of
  `F0`, for `central_mass` 2 and 2.5, and for `sheath_Lambda` 3 and 2.84.
- Every Jacobian column matches a central finite difference of Artola eq. (6) to 1e-5 relative.
- The row stays finite 300 e-foldings into both the ion-saturation and the electron branch.
- The assembly block compiles against stub dependencies with the real `mod_assembly` signatures.

**Not done:** MPI build, cluster run, regression comparison, or any statement about long-time
stability. No run of any kind has exercised this code.
