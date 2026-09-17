# Audit of the sheath current BC, current-slot form (2026-09-17)

Branch `sheath-j-clean` @ baffa84ee. Form under audit: `bcs(1)%sheath_j` with `sheath_j_current_row`:
the characteristic zj = j_sat*(1 - exp(x)) is the zj row at wall Gauss points where the field is not
grazing, u is left to the vorticity equation, psi and w Dirichlet, marginal weak Bohm row, inflow closure
on. It is the structure of the old `sheath_zj_weak` route that ran ~3900 steps on type 1 with an ion tail
and `use_sc = 10`.

## Evidence so far (type 1 alone, no use_sc, no slopes, from t = 0)

| step | j/j_sat range | |j|>j_sat | Phi [V] | Inet/Isat | state |
|---|---|---|---|---|---|
| 10  | -1.9 .. +0.96 | 0.52 | 74 .. 234   | -0.98 | current inside the characteristic from the first solve |
| 100 | -2.7 .. +1.00 | 0.54 | 66 .. 683   | -1.33 | ion rail reached at the outer target, Phi rising there |
| 178 | -7.2 .. +1.00 | 0.49 | 38 .. 2600  | -0.51 | transient peak; min rho 2e-4 at the 1/3 junction |
| 280 | -6.8 .. +1.00 | 0.06 | 73 .. 410   | +0.17 | recovered: interior inside the characteristic |
| ~400 | -26 .. +47   | 0.06 | -468 .. 344 | +0.10 | electron rail at the outer target, Phi below the wall, rho < 0 there |
| 438 | blow-up | | | | |

Every earlier run of this structure had types 4 and 9 on and died in 4-19 steps at their corner; the
old branch found those two types fail alone, whatever the formulation. The structure itself, on type 1,
anchors the potential and drives the current to the characteristic; it fails only where the plasma
demands a current outside [-(e^Lambda - 1), +1] j_sat, at both rails.

## Code audit (what is assembled, and what is wrong or fragile)

1. **Row form.** Zbig*dl*(zj - j_sat f) added to the volume current-definition row at boundary test
   functions (volume entries ~1e3 against ~3e9: 1e-6 relative pollution). The old route cleared the row
   and replaced it by a mass-matrix projection. Equivalent for type 1; on ill-framed nodes (4, 9) the
   leftover volume row is the ill-conditioned Delta*psi and should be cleared. **Do the exact
   replacement** (clear zj row across adjacent blocks, project, equilibrate) when going beyond type 1.
2. **Signs and normalisation.** j_sat = c_sat*rho*sign(B.n)*cs/|B| reproduces Artola eq. 5 with the
   marginal Bohm outflow; -zj*(B_pol.n)/F0 = e*n*cs*|b.n|*f > 0 into the wall for both F0 signs (tested);
   x = Lambda - a_n*(u - C_V*V_wall)/(2Te) is the floating row's normalisation (root test). Sound.
3. **Raw Te and rho in the row.** Te0 and r0 enter without the corr_neg map that every other natural
   row uses (cs0 already uses Ti0_corr + Te0_corr). Te0 <= 0 flips the sign of x; rho <= 0 flips j_sat.
   Use Te0_corr/r0_corr with dcorr_neg_temp_dT1 in the Te column. **Robustness fix, small.**
4. **Rails.** f -> 1 as x -> -inf and f = 1 - e^Lambda for x >= Lambda; on both rails df/du = 0, the row
   has no u column, and wherever the plasma delivers more than the rail the vorticity equation drives Phi
   without limit. Measured both ways. `sheath_j_ion_slope` / `sheath_j_e_slope` give a finite conductance
   there (default 0). The ion slope is physical (sheath expansion); the electron one is a regularisation of
   the model's inability to let the circuit reduce the demand fast enough.
5. **Coupling path is complete.** The wall zj enters the vorticity row through v*[psi, zj] (not integrated
   by parts, no dropped surface term), the potential enters the induction row through [u, psi] in the
   last element, and the current at the wall relaxes on the resistive time: the 178 -> 280 recovery is
   this loop working. The dropped surface terms of the retained vorticity row (polarisation, viscosity)
   are the physical statement "no perpendicular current into the wall".
6. **Junction node type 1/3.** Its zj is released (per-node decision) and gets the current row from the
   type-1 edge only (Zbig-dominant, so determined). Its u is from the vorticity row while the type-3
   neighbour is Dirichlet floating and type 2 beyond is u = 0: a potential step over one element at the
   end of the target. It is also the coldest, most rarefied node (rho 2e-4, Ti 8 eV) and hence the
   smallest j_sat: the rails are reached there first. The old branch had the same junction.
7. **Inflow closure under a sheath potential.** In D (floating) type 1 had 0% inflow, so the closure was
   inert on the targets. With the sheath potential the wall ExB has inward regions on 10-50% of type 1
   (0.13 at 100, 0.48 at 178, 0.23 at 280), and the closure relaxes rho toward 0 there at the rate |vn|.
   That drains the last element exactly where j_sat is already small and pushes the plasma toward the
   rails. The closure was built for the floating-u grazing wall; on a current-carrying target its datum
   rho_in = 0 is questionable (the wall is a neutral source, not a plasma sink). **Test
   `mach1_weak_inflow = .f.` with the sheath.** Candidate driver of the 438 crash.
8. **Bohm row.** Marginal form with u free is consistent (no u column). The drift form's Zbig u column
   on a free wall potential was implicated at the 4/9 corners; on type 1 alone it is untested. Stay
   marginal until the current row is stable; the drift can be added as a flux on rho/T (plan R5).
9. **Ramp.** `sheath_j_ramp` acts only on the potential-form row. The current-slot form starts at full
   gain and its transient (2600 V) recovered; a ramp zj = zj_pinned + alpha*(j_sat f - zj_pinned) is
   cheap if the start proves marginal.
10. **Energy channel.** gamma_sheath_e is fixed; with current the electron transmission depends on the
    drop (SOLPS B.52, q_e ~ Gamma_e*(2Te + e*DeltaPhi)). At an electron-saturated spot the model dumps
    no extra heat, so Te there does not respond to the current it draws. Physics item, second order for
    stability, first order for the target temperature under current.
11. **Lambda fixed at 3.** The old route had Lambda(Ti/Te); the electron amplitude scaling with cs
    (sqrt(Ti+Te)) instead of sqrt(Te) is inherited. Second order.
12. **n_tor > 1.** The row is evaluated per plane on the real state; fine. Not run.
13. **Diagnostics.** The table reports the max |j/j_sat| location; it does not yet report WHERE Phi
    min/max sit or the current at the strike points vs the junction. Add (R,Z) of Phi extremes and a
    split by R (inner/outer target) as the old route had.

## Model audit (what the model can and cannot do)

- The wall current at t = 0 is the equilibrium's residual current at 4e17 m^-3, i.e. +-1.5 j_sat with
  sign changes; a sheath cannot carry it and the induction equation relaxes it on the resistive time.
  The relaxation works (178 -> 280) but passes through the ion rail. Higher target density (realistic
  divertor, 1e19) makes j_sat 25x larger and the whole transient benign.
- Where the plasma demands beyond a rail, the physical answer is the flux-tube circuit reducing the
  demand through the parallel field; the model has this path (induction row with [u,psi]) but it is
  slower than the algebraic row's response. The slopes buy the time; the honest alternative is a
  smaller demand (density) or a current ramp.
- The two-target thermoelectric loop closes through the vessel at V_wall = 0; Inet/Isat is the net
  current to the vessel per type, +0.17 at 280 with a 25/49 eV target contrast. A per-tile floating
  vessel is a later scalar constraint.

## Plan for boundary type 1

1. Robustness: corr_neg on Te and rho in the row (item 3); exact row replacement optional here.
2. Runs from t = 0, type 1 only, marginal Bohm, no use_sc:
   - C: ion slope 0.03, D: ion + electron slope 0.03 (queued);
   - F: no slopes, `mach1_weak_inflow = .f.` (item 7);
   - E: no slopes, use_sc = 10 (the old route's control).
   Read: Phi min/max and |j|>j_sat at 300/600/1200, Inet/Isat sign vs the target Te contrast, min rho at
   the junction.
3. If F holds without slopes, the inflow closure is the sheath-run driver and the slopes stay at 0;
   if only C/D hold, the rails are the limit and the slopes are the model statement, documented as such.
4. Then the physics items: drop-dependent electron transmission (10), Lambda(Ti/Te) (11), the
   (R,Z)/target-split diagnostics (13).

## Plan for the other node types

- **Type 5 (upper HFS wall)**: the old route ran 1+5 for 305 steps converged; expect it to work as type 1
  does once type 1 holds. Add it second.
- **Type 3 (target end / outer-boundary corner)** and **type 2**: make the potential continuous across
  the junction. Type 3 as a sheath type gets its current row from the 1-3 edge (the 3-2 edge stays
  excluded); type 2 should carry the floating potential of its own Te instead of u = 0, with the absolute
  gauge then fixed by V_wall through the sheath. Removes the 100 V step at the end of the target.
- **Types 4 and 9 (ray-cast wall extension)**: fail alone on both branches; their node frames are
  degenerate (det 0.05-0.38) so Delta*psi at the node is ill-conditioned and releasing zj amplifies it.
  Two routes: (a) keep them floating with pinned current (1% of the wetted area, the old campaign's
  choice); (b) fix the frames in the grid builder (`grid_xpoint_wall.f90`, the ray-cast extension's
  x(1,2,:)/x(1,3,:)), which is the root cause and benefits every BC. Exact row replacement (item 1) is
  a prerequisite for even trying (b) on them.
- Sort by incidence stays as it is: the row is never imposed below sin(min_sheath_angle); a wall the
  field does not reach carries no sheath current in this model.
