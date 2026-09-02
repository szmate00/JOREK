# The sign of `a_n` in the sheath characteristic

Companion to `doc/jorek_solps_parallel_ohm.md`, which establishes the convention this rests on.
Written 2026-09-01 after the observation that `u` — and therefore `Phi` — is negative through
essentially the whole SOL, which is backwards for a sheath.

**Conclusion: `a_n` has the wrong sign, `c_sat` does not, and the fix is two characters.**

---

## 1. The observation that started it

`u` at cycle 1000: `min = -1264`, `max = +5.4`, negative through the SOL and the core.
`Phi = +F0*u` with `F0 = +2.972`, so `Phi < 0` almost everywhere.

That is backwards. Electrons are more mobile than ions, so a floating or grounded wall charges
NEGATIVE relative to the plasma: the SOL plasma should sit at `Phi = +Lambda*Te/e` above the
wall, i.e. POSITIVE, order tens of volts.

Meanwhile `sheath_diag_report` prints `ePhi/kTe = +2.97`, which looks correct. Both cannot be
right, and the diagnostic uses the same coefficient as the boundary condition, so it cannot be
the independent check it appears to be.

---

## 2. Where the sign error entered

`mod_sheath_bc.f90:26-32` states the reasoning:

> the JOREK reference paper (Hoelzl et al 2021, eq. 26) defines `u = Phi/F0` together with a
> velocity ansatz `v_pol = -R grad(u) x e_phi`, whereas the model600 element matrix implements
> `v_pol = +R grad(u) x e_phi` ... The code's u is therefore minus the paper's, which is where the
> minus sign in `a_n` comes from.

The code's velocity components are `v_R = -R u_Z`, `v_Z = +R u_R`
(`particles/mod_fields.f90:582-583`). Now evaluate `grad(u) x e_phi` in each basis:

**JOREK's basis, (e_R, e_Z, e_phi) right-handed** (wiki: `x = R cos(phi)`, `y = -R sin(phi)`,
"phi goes clockwise if looked at from above"):

    e_R x e_phi = -e_Z ,   e_Z x e_phi = +e_R
    grad(u) x e_phi = u_R(-e_Z) + u_Z(+e_R) = ( u_Z , -u_R )
    -R grad(u) x e_phi = ( -R u_Z , +R u_R )        <-- matches the code

**Standard basis, (e_R, e_phi, e_Z) right-handed:**

    e_R x e_phi = +e_Z ,   e_Z x e_phi = -e_R
    grad(u) x e_phi = ( -u_Z , +u_R )
    +R grad(u) x e_phi = ( -R u_Z , +R u_R )        <-- ALSO matches the code

**The same component expression is `-R grad(u) x e_phi` in JOREK's basis and
`+R grad(u) x e_phi` in the standard one.** The header compared the code's components against the
paper's formula while implicitly using the standard basis, found a sign discrepancy, and
attributed it to `u` rather than to the basis. In JOREK's own basis the code implements the
PAPER's ansatz exactly, so `u_code = u_paper = Phi/F0`, i.e. **`Phi = +F0*u`**.

This is the same trap recorded in `doc/jorek_solps_parallel_ohm.md` section 8.1: component
identities are handedness-independent, cross products are not.

---

## 3. Deriving `a_n` from first principles

Unit relations, both verified in the code:

    Phi_SI [V]  = F0 * u / sqrt(mu0*rho0)          (mod_expression.f90, res = u0*F0/fact_time)
    T_J         = mu0 * n0 * kT_SI                  (=> T_J = 2.036e-5 * T_eV here; the code's
                                                     documented 2.01e-5*central_density = 2.032e-5)

with `n0 = central_density*1e20`, `m_i = central_mass*ATOMIC_MASS_UNIT`, `rho0 = n0*m_i`.

Then, with `V_wall = 0`:

    e*Phi/kTe = e * Phi_SI / kTe_SI
              = e * [ F0*u/sqrt(mu0*rho0) ] / [ T_J/(mu0*n0) ]
              = [ e*mu0*n0*F0 / sqrt(mu0*rho0) ] * u / T_J

The code computes `e*Phi/kTe = ( 0.5*a_n*u - vw ) / Te_l`, so

    0.5*a_n = e*mu0*n0*F0 / sqrt(mu0*rho0)

and since `mu0*n0 / sqrt(mu0*n0*m_i) = sqrt(mu0*rho0)/m_i`,

    a_n = + 2*e*F0*sqrt(mu0*rho0) / m_i        <-- POSITIVE

The code has `a_n = - 2*e*F0*sqrt(mu0*rho0)/m_i`. **Magnitude verified identical to 12
significant digits** (1.856385e+02 both ways). Only the sign differs.

---

## 4. Why flipping `a_n` alone would be wrong

`c_sat = -a_n/2`, and `zj_sat = c_sat * rho * g_bn * cs / Btot`. Flipping `a_n` with that formula
intact flips `c_sat`, hence `zj_sat`, hence the direction of the sheath current.

**We do not want that.** The current side has been validated independently of the potential
convention:
- `I_sheath == I_Ampere` to four significant figures on both targets and the wall total;
- the `pl_sh` relative-sign error was found and fixed in July against the halo-current
  diagnostic as an external anchor, not against this module.

The two signs answer different questions. `a_n`'s sign encodes the **potential** convention
(`Phi` vs `u`). `c_sat`'s sign encodes the **current** convention (`j` vs `zj = Delta*psi`). The
identity `|c_sat| = |a_n|/2` is a magnitude relation from the shared normalisation; the signs are
independent, and writing `c_sat = -a_n/2` silently tied them together. With `a_n` wrong, that
coupling delivered the RIGHT `c_sat` — two errors cancelling, which is why the currents look sane
and only the potential is visibly wrong.

---

## 5. The fix

In `sheath_norm` (`models/model600/mod_sheath_bc.f90`):

```fortran
! --- CURRENT
a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
c_sat = - 0.5d0 * a_n

! --- CORRECTED
a_n   = + 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
c_sat = + 0.5d0 * a_n
```

`c_sat` is numerically UNCHANGED. `a_n` flips. Consequences:

| quantity | change |
|---|---|
| `x = (0.5*a_n*u - vw)/Te - lam` | correct sign in `u`: the BC now floats the wall at `u > 0`, i.e. `Phi > 0` |
| `zj_sat` | **unchanged** — current magnitude and direction untouched |
| `dzj_du = zj_sat*fp*dx_du`, `dx_du = 0.5*a_n/Te` | flips sign, consistently with the residual |
| `sh_u_float = 2*(Te*lam + vw)/a_n` | becomes positive — the floating `u` is now positive |
| `sheath_init_potential` | uses `sheath_norm`, so it initialises to positive `u` automatically |
| `vw = e*V_wall*mu0*n0` | unchanged and correct; `X` is the drop `V_sheath - V_wall` |

Also fix the stale header text at `:26-32` and the same claim in `doc/sheath_bc_whitepaper.tex`
and `doc/HANDOFF.md`, or the error regenerates - a reviewer reproduced it from those comments
this week.

---

## 6. How to verify, in order

1. **The current must NOT change.** Same restart, one step, compare per-target `I_sheath` and
   `I_Ampere` before and after. If they flip sign, the `c_sat` reasoning in section 4 is wrong and
   the patch must be reconsidered - this is the test that can falsify it.
2. **`u` must become positive through the SOL.** Re-plot the field that started this.
3. **`ePhi/kTe` should still read ~ +Lambda at float.** The diagnostic shares `a_n`, so the
   NUMBER is unchanged; what changes is that it now agrees with the physical `Phi`.
4. **Expect a violent transient.** The restart sits at `u ~ -1264`; the corrected equilibrium is
   `u > 0`. That is a large perturbation, not a small correction. Use `sheath_ramp_time`, or
   restart from a state initialised with `sheath_init_u` under the corrected sign.

---

## 7. What this does and does not invalidate

**Unchanged:** all current magnitudes and directions; `I_sheath == I_Ampere` closure; the weak
Galerkin row-replacement result; the magnitudes of `ePhi/kTe`; every conclusion that depended only
on `|Phi|` or on currents.

**Sign-flipped:** the physical sign of `Phi` everywhere, hence `E_r`, hence the direction of the
`ExB` drift, hence the direction of any PFR transfer inferred from it, and the polarity of the
imposed bias in the `sheath_V_wall_asym` scan. `hfshd-bias-scan-negative` should be re-read with
that in mind: a scan that moved the ratio the "wrong" way may have been applying the bias in the
opposite direction to what was intended.
