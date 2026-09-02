# The sheath `a_n` sign: evidence for review

Self-contained. Every claim below is a code snippet you can check. Nothing here relies on the
module header of `mod_sheath_bc.f90`, on `doc/sheath_bc_whitepaper.tex`, or on `doc/HANDOFF.md` -
those three all state `Phi = -F0*u` and are the artefacts under suspicion.

---

## 1. What the code does today

`models/model600/mod_sheath_bc.f90`, `sheath_norm`:

```fortran
  ! --- NOTE the minus sign: Phi = -F0*u in the code's variables (see the module header)
  a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
  c_sat = - 0.5d0 * a_n
  vw    =   EL_CHG * V_wall_use * MU_ZERO * central_density * 1.d20
```

so today `a_n < 0` and `c_sat > 0`. These feed

```fortran
  x      = ( 0.5d0 * a_n * u - vw ) / Te_l - lam      ! X = e*Phi/(k*Te) - Lambda
  zj_sat = c_sat * rho_l * g_eff * cs / Btot
```

---

## 2. The original implementation had BOTH signs the other way

`bb5d7fc31`, "Sheath j-V boundary condition for the electric potential (model600)",
17 Aug 11:00, `models/model600/mod_boundary_conditions.f90:808-810`:

```fortran
            sh_a_n   = 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * central_density * 1.d20 * ...) &
                     / (central_mass * ATOMIC_MASS_UNIT)
            sh_c_sat = - 0.5d0 * sh_a_n
```

`a_n > 0`, `c_sat < 0`.

`0f6c8b881`, "j bc first go", 17 Aug 16:44 - five hours later - changed only `a_n`:

```diff
-            sh_a_n   = 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * ...)
+            ! --- NOTE the leading minus: the electrostatic potential is Phi = -F0*u in the code's
+            sh_a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * ...)
```

`c_sat = -0.5*a_n` was left untouched, **so `c_sat` flipped as a side effect**. One deliberate
change, two coefficients moved.

---

## 3. The original derivation is STILL IN THE CODE, and it says `Phi = +F0*u`

`models/model600/mod_boundary_conditions.f90:842-854` - never deleted by the 17 Aug change:

```fortran
          ! --- J. Artola, "Sheath boundary conditions for the electric potential in JOREK".
          ! --- Stangeby 2.68 with the sheath factor:
          !
          ! ---     j = j_sat * f ,   f = 1 - exp(-X) ,   X = e*Phi/(k*Te) - Lambda
          !
          ! --- Phi = V_sheath_entrance - V_wall is referenced to the WALL, so f = 0 (zero current)
          ! --- sits at Phi = Lambda*Te/e, the floating potential, and Phi = 0 is electron
          ! --- saturation. In JOREK units, with Phi = F0*u,
          !
          ! ---     e*Phi/(k*Te) = ( a_n*u/2 - e*V_wall*mu0*n0 ) / Te
          ! ---     a_n   = 2*e*F0*sqrt(mu0*rho0)/(m_c*m_amu)
          ! ---     j_sat = c_sat*rho*v_par ,  c_sat = -e*F0*n0*sqrt(mu0/rho0) = -a_n/2
```

Three things to note:

1. **`with Phi = F0*u`** - positive, stated explicitly.
2. `a_n = 2*e*F0*sqrt(mu0*rho0)/m_i` - **positive**.
3. `c_sat = -e*F0*n0*sqrt(mu0/rho0)` is a **closed form with its own explicit minus**, given
   BEFORE the `= -a_n/2` identity. And `-e*F0*n0*sqrt(mu0/rho0) == -e*F0*sqrt(mu0*rho0)/m_i`
   since `rho0 = n0*m_i`. So `c_sat = -a_n/2` is a **signed** identity in the `Phi = +F0*u`
   convention, not a magnitude relation - `c_sat` is negative and its sign does not depend on
   `a_n`'s.

---

## 4. Why the 17 Aug change was made, and why the reason is wrong

The justification, in `mod_sheath_bc.f90:26-32`:

> the JOREK reference paper (Hoelzl et al 2021, eq. 26) defines `u = Phi/F0` together with a
> velocity ansatz `v_pol = -R grad(u) x e_phi`, whereas the model600 element matrix implements
> `v_pol = +R grad(u) x e_phi` ... The code's u is therefore minus the paper's.

The code's velocity, `particles/mod_fields.f90:581-582`:

```fortran
v_R   = -R * U_Z + vpar * R_inv *psi_Z
v_Z   =  R * U_R - vpar * R_inv *psi_R
```

i.e. `v_pol = (-R u_Z, +R u_R)`. Now evaluate `grad(u) x e_phi` in each basis:

| basis | cross products | `grad(u) x e_phi` | matches the code as |
|---|---|---|---|
| JOREK `(e_R, e_Z, e_phi)` RH | `e_R x e_phi = -e_Z`, `e_Z x e_phi = +e_R` | `( u_Z, -u_R )` | **`-R grad(u) x e_phi`** |
| standard `(e_R, e_phi, e_Z)` RH | `e_R x e_phi = +e_Z`, `e_Z x e_phi = -e_R` | `( -u_Z, +u_R )` | **`+R grad(u) x e_phi`** |

**The same components read as `-R grad(u) x e_phi` in JOREK's basis and `+R grad(u) x e_phi` in
the standard one.** The 17 Aug reasoning compared the code's components against the paper while
implicitly using the standard basis, found a sign difference, and attributed it to `u` rather
than to the basis. In JOREK's own basis the code implements the paper's ansatz exactly, so
`u_code = u_paper = Phi/F0`.

**Which basis is JOREK's?** Three independent checks, none of them a comment:

- `diagnostics/jorek2_povray.f90:102-103`, `vacuum/mod_vacuum_fields.f90:1172-1173`:
  `x = R cos(phi)`, **`y = -R sin(phi)`**. Hence `e_R x e_Z = e_phi`: `(R,Z,phi)` is right-handed
  and JOREK's `e_phi` is minus the standard one. (Same as the JOREK wiki's "phi goes CLOCKWISE if
  looked at from above".)
- `particles/mod_fields.f90:590-596` `rot_tmp` is the right-handed curl for `(R,Z,phi)` with
  scale factors `(1,1,R)`. Applied to the code's own `A` (`:775`) it returns the code's own `B`
  (`:791`) exactly. In the standard basis the same input returns `-B` in all three components, so
  this test discriminates.
- `diagnostics/new_diag/mod_expression.f90:1354`: `Jtor = -zj0/BigR`, i.e. `mu0*j_phi = -zj/R`,
  which is the third curl component in that basis.

---

## 5. Independent check of the magnitude

From `Phi_SI = F0*u/sqrt(mu0*rho0)` and `T_J = mu0*n0*kT_SI`:

```
0.5*a_n = e*mu0*n0*F0/sqrt(mu0*rho0)  =>  a_n = 2*e*F0*sqrt(mu0*rho0)/m_i
```

Numerically at your parameters: **1.856385e+02 both ways, 12 significant digits.** So the
magnitude in the code is right; only the sign moved.

---

## 6. A second, compensating error in the diagnostic

`models/model600/mod_sheath_diag.f90:159-161`:

```fortran
  ! --- J.n = J_par*(b.n) = zj*(B.n)/(F0*mu0)   [A/m^2]
  jn_sheath = zj_sh * Bdotn / (F0 * MU_ZERO)
  jn_amp    = zj0   * Bdotn / (F0 * MU_ZERO)
```

From `mu0*j_phi = -zj/R` and `B_phi = F0/R`, a field-aligned current is
`j = (j_phi/B_phi)*B = -(zj/(mu0*F0))*B`, so `j.n = -zj*(B.n)/(mu0*F0)`. **The minus is missing.**

With today's `c_sat > 0` the two sign errors cancel and the reported `I_sheath` looks physically
correct - which is exactly why "the currents look sane" and why this survived. **They must be
fixed together**, or the report inverts once `a_n` is corrected.

---

## 7. What is NOT evidence

- **`I_sheath == I_Ampere` to four significant figures.** Both are computed with the same
  formula at `:160-161`, differing only in `zj_sh` vs `zj0`. The closure tests that the row drove
  `zj_sh -> zj0`. It is invariant under flipping `c_sat` and says nothing about the sign.
- **`ePhi/kTe ~ +Lambda` at float.** The diagnostic uses the same `a_n` as the boundary
  condition, so it reports `+Lambda` under either sign convention.
- **The impurity trapping.** It shows the code produced `Phi < 0` in the divertor - the symptom,
  not a diagnosis. Under the corrected sign a floating wall gives `Phi = +Lambda*Te/e`, a hill
  that expels positive ions, so that well is an artefact of the bug rather than support for it.

---

## 8. The fix

```fortran
! models/model600/mod_sheath_bc.f90, sheath_norm - revert 0f6c8b881
a_n   = + 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
c_sat = - 0.5d0 * a_n            ! UNCHANGED - c_sat becomes negative again, as originally

! models/model600/mod_sheath_diag.f90:159-160 - the compensating error
jn_sheath = - zj_sh * Bdotn / (F0 * MU_ZERO)
jn_amp    = - zj0   * Bdotn / (F0 * MU_ZERO)
```

Plus the stale comments: `mod_sheath_bc.f90:26-32`, `doc/sheath_bc_whitepaper.tex`,
`doc/HANDOFF.md`.

---

## 9. The measurement that can still falsify it

Zero code change; both quantities already exist and neither goes through `mod_sheath_bc`:

```fortran
diagnostics/new_diag/mod_expression.f90:1356-1357
  Jpar        = (JpolR*BR + JpolZ*BZ + Jtor*Btor) / Btot * sign(1.d0, F0)
  Jpar_ionsat = r0 * vpar0 * Btot
```

`Jpar` is built from `Jtor = -zj0/R`; `Jpar_ionsat` is `e*n*v_par`. At a strike point in ion
saturation, with `F0 > 0`, **they must have the same sign** if the sheath is passing physical ion
saturation current.

Prediction: they are OPPOSITE today, and agree after the fix. **If they already agree, the fix is
wrong** - do not apply it.

Secondary: new_diag `Phi` (`mod_expression.f90:1757-1758`) contains no `a_n`. At a near-floating
strike point it must read `~ +Lambda*Te[eV]` volts. It should read `~ -Lambda*Te` today.
