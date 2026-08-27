#!/usr/bin/env python3
"""
analyse_hill.py -- the potential hill: Delta_Phi from the X-point to each strike point.

Reads pol_line output from postproc_hill.in and reports, per leg,

    Delta_Phi   = Phi(X-point) - Phi(plate)          measured
    int eta*j_par dl                                 predicted by Ohm's law

If those agree, the hill IS the resistive drop needed to close the divergent grad-B
current through the plates -- the Senichenkov (CPP 2022) HFSHD mechanism -- confirmed in
your own run rather than by analogy. If Delta_Phi is much smaller than the integral, the
current is not closing where you think it is.

    util/analyse_hill.py pol_line_*.dat

Because Phi_hill is LINEAR in eta and eta ~ T_e^-3/2, the script also reports what the
hill would become at a colder leg -- the single most useful extrapolation here, since the
whole question is whether a detached leg produces a SOLPS-scale hill.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyse_targets import read_boundary_file  # noqa: E402


def leg(d):
    """Everything for one traverse."""
    # SIGN FIX: postproc's `Phi` returns +F0*u, but model600 implements
    # v_pol = +R grad(u) x e_phi, i.e. Phi = -F0*u. Same fix as analyse_targets.
    phi = -d["Phi"]
    R, Z = d["R"], d["Z"]
    dl = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(R), np.diff(Z)))])

    B = d["B_abs"]
    jpar = np.where(B > 0,
                    (d["JR"] * d["BR"] + d["JZ"] * d["BZ"] + d["Jtor"] * d["Btor"]) / B,
                    np.nan)
    eta = d.get("eta_T", np.full_like(B, np.nan))

    # The field line is not the traverse: only the component of the path along b_pol
    # contributes to a PARALLEL integral. Project it.
    Bpol = np.hypot(d["BR"], d["BZ"])
    frac = np.where(B > 0, Bpol / B, np.nan)          # dl_pol -> dl_par is /frac
    integrand = eta * jpar / np.where(frac > 1e-6, frac, np.nan)
    good = np.isfinite(integrand)
    ohm = np.trapz(integrand[good], dl[good]) if good.sum() > 2 else np.nan

    return dict(dl=dl, phi=phi, Te=d.get("T_e"), eta=eta, jpar=jpar, ohm=ohm,
                R=R, Z=Z, span=float(np.nanmax(phi) - np.nanmin(phi)),
                dphi=float(phi[0] - phi[-1]))


def fmt(x, p=3):
    if x is None or not np.isfinite(x):
        return "n/a"
    a = abs(x)
    return "%.*e" % (p, x) if (a >= 1e4 or (0 < a < 1e-2)) else "%.*f" % (p, x)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    print("=" * 78)
    print("POTENTIAL HILL   Delta_Phi = Phi(X-point) - Phi(plate)")
    print("=" * 78)
    print("  Phi sign FLIPPED on read: postproc returns +F0*u, model600 needs -F0*u.")

    legs = []
    for path in sys.argv[1:]:
        for step, t_now, d in read_boundary_file(path):
            if "Phi" not in d:
                continue
            L = leg(d)
            legs.append((os.path.basename(path), step, t_now, L))

    for name, step, t_now, L in legs:
        print("")
        print("  %s   step %s   t_now %s" % (name, step, fmt(t_now)))
        print("    from (%.3f, %.3f) to (%.3f, %.3f),  length %.4f m"
              % (L["R"][0], L["Z"][0], L["R"][-1], L["Z"][-1], L["dl"][-1]))
        print("    Delta_Phi (start - end)   = %s V" % fmt(L["dphi"]))
        print("    Phi range along the line  = %s V" % fmt(L["span"]))
        print("    int eta*j_par dl          = %s V   (Ohm's law prediction)" % fmt(L["ohm"]))
        if np.isfinite(L["ohm"]) and abs(L["ohm"]) > 1e-12:
            print("    ratio measured/predicted  = %s" % fmt(L["dphi"] / L["ohm"]))
        if L["Te"] is not None:
            te = L["Te"][np.isfinite(L["Te"])]
            if te.size:
                temin, temax = float(np.min(te)), float(np.max(te))
                print("    T_e along the line        = %s .. %s eV" % (fmt(temin), fmt(temax)))
                # Phi_hill is LINEAR in eta and eta ~ T^-3/2, so this extrapolation is
                # exact under the assumption that j_par is unchanged -- which is the
                # point: j_par is set by div.j = 0, not by Ohm's law.
                if temin > 0:
                    for tcold in (5.0, 1.0):
                        if temin > tcold:
                            f = (temin / tcold) ** 1.5
                            print("      -> at a %.0f eV leg, hill would be ~%s V  (x%.0f)"
                                  % (tcold, fmt(L["dphi"] * f), f))

    print("")
    print("  SOLPS (Senichenkov CPP 2022) gets a hill of order 100 V. That is the number")
    print("  to beat. Phi_hill is linear in eta and eta ~ T_e^-3/2, so the leg temperature")
    print("  is the whole amplifier: 25 eV -> 1 eV is a factor 126.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
