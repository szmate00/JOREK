#!/usr/bin/env python3
"""Plot grazing-incidence diagnostics from JOREK boundary_quantities output.

The input file must contain these quantities (the standard length/R/Z columns
may also be present):

  theta_geo Bnorm B_abs vpar vparB_norm vExB_norm vflow_norm V_sound

Several files can be supplied to compare configurations at the same timestep.
All velocity quantities produced by jorek2_postproc are in SI units (m/s),
while Bnorm/B_abs and the derived Mach number are dimensionless.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


REQUIRED = (
    "theta_geo",
    "Bnorm",
    "B_abs",
    "vpar",
    "vparB_norm",
    "vExB_norm",
    "vflow_norm",
    "V_sound",
)


def read_boundary_file(path: Path) -> tuple[dict[str, np.ndarray], str]:
    """Read a whitespace boundary_quantities table using its comment header."""
    header = None
    time_label = ""
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped.startswith("#"):
                continue
            text = stripped[1:].strip()
            if text.startswith("time step"):
                time_label = text
            elif header is None and "theta_geo" in text:
                header = text.split()

    if header is None:
        raise ValueError(f"{path}: could not find the column header")

    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(header):
        raise ValueError(
            f"{path}: header has {len(header)} columns but data has "
            f"{values.shape[1]}"
        )

    columns = {name.lower(): values[:, i] for i, name in enumerate(header)}
    missing = [name for name in REQUIRED if name.lower() not in columns]
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
    return columns, time_label


def col(data: dict[str, np.ndarray], name: str) -> np.ndarray:
    return data[name.lower()]


def finite_log_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(arrays[0].shape, dtype=bool)
    for array in arrays:
        mask &= np.isfinite(array) & (array > 0.0)
    return mask


def print_summary(label: str, data: dict[str, np.ndarray]) -> None:
    bn_ratio = col(data, "Bnorm") / col(data, "B_abs")
    vpar = col(data, "vpar")
    vpar_n = col(data, "vparB_norm")
    vexb_n = col(data, "vExB_norm")
    vflow_n = col(data, "vflow_norm")
    cs = col(data, "V_sound")

    normal_identity = vpar - np.divide(
        vpar_n,
        bn_ratio,
        out=np.full_like(vpar_n, np.nan),
        where=np.abs(bn_ratio) > 1.0e-12,
    )
    closure = vflow_n - vpar_n - vexb_n
    good_identity = np.isfinite(normal_identity)

    print(f"\n{label}")
    print(f"  samples: {len(vpar)}")
    print(f"  min |Bn/B|:              {np.min(np.abs(bn_ratio)):.6e}")
    print(f"  max |Vpar|:              {np.max(np.abs(vpar)):.6e} m/s")
    print(f"  max |Vpar|/cs:           {np.nanmax(np.abs(vpar/cs)):.6e}")
    print(f"  max |Vpar Bn/B|:         {np.max(np.abs(vpar_n)):.6e} m/s")
    print(f"  max |vExB,n|:            {np.max(np.abs(vexb_n)):.6e} m/s")
    print(f"  max |total normal flow|: {np.max(np.abs(vflow_n)):.6e} m/s")
    print(f"  max normal-flow closure: {np.max(np.abs(closure)):.6e} m/s")
    if np.any(good_identity):
        print(
            "  max Vpar identity error:  "
            f"{np.max(np.abs(normal_identity[good_identity])):.6e} m/s"
        )


def plot_profiles(runs: list[tuple[str, dict[str, np.ndarray]]]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

    for label, data in runs:
        theta = np.degrees(col(data, "theta_geo"))
        bn_ratio = col(data, "Bnorm") / col(data, "B_abs")
        vpar = col(data, "vpar")
        cs = col(data, "V_sound")
        vpar_n = col(data, "vparB_norm")
        vexb_n = col(data, "vExB_norm")
        vflow_n = col(data, "vflow_norm")

        # Scatter is intentional: theta_geo is not necessarily a one-to-one
        # boundary coordinate on a re-entrant or artificial-PFR boundary.
        axes[0].scatter(theta, bn_ratio, s=4, label=label)
        axes[1].scatter(theta, vpar / 1.0e3, s=4, label=f"{label}: Vpar")
        axes[1].scatter(theta, cs / 1.0e3, s=4, alpha=0.5,
                        label=f"{label}: cs")
        axes[2].scatter(theta, vpar_n / 1.0e3, s=4,
                        label=f"{label}: Vpar Bn/B")
        axes[2].scatter(theta, vexb_n / 1.0e3, s=4,
                        label=f"{label}: vExB,n")
        axes[2].scatter(theta, vflow_n / 1.0e3, s=4,
                        label=f"{label}: total")
        axes[3].scatter(theta, vpar / cs, s=4, label=label)

    axes[0].axhline(0.0, color="black", linewidth=0.6)
    axes[0].set_ylabel(r"$B_n/B$")
    axes[1].set_ylabel("velocity [km/s]")
    axes[2].set_ylabel("normal velocity [km/s]")
    axes[3].set_ylabel(r"$V_\parallel/c_s$")
    axes[3].set_xlabel(r"geometric angle $\theta$ [deg]")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(loc="best", markerscale=2)
    fig.suptitle("Floating-u boundary profiles")
    fig.tight_layout()


def plot_grazing_correlations(
    runs: list[tuple[str, dict[str, np.ndarray]]]
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for label, data in runs:
        bn = np.abs(col(data, "Bnorm") / col(data, "B_abs"))
        vpar = np.abs(col(data, "vpar"))
        mach = np.abs(col(data, "vpar") / col(data, "V_sound"))
        vpar_n = np.abs(col(data, "vparB_norm"))
        vexb_n = np.abs(col(data, "vExB_norm"))
        vflow_n = np.abs(col(data, "vflow_norm"))

        mask = finite_log_mask(bn, vpar)
        axes[0, 0].scatter(bn[mask], vpar[mask], s=6, alpha=0.55,
                           label=label)
        mask = finite_log_mask(bn, mach)
        axes[0, 1].scatter(bn[mask], mach[mask], s=6, alpha=0.55,
                           label=label)
        mask = finite_log_mask(vexb_n, vpar_n)
        axes[1, 0].scatter(vexb_n[mask], vpar_n[mask], s=6, alpha=0.55,
                           label=label)
        mask = finite_log_mask(bn, vflow_n)
        axes[1, 1].scatter(bn[mask], vflow_n[mask], s=6, alpha=0.55,
                           label=label)

    axes[0, 0].set_xlabel(r"$|B_n|/B$")
    axes[0, 0].set_ylabel(r"$|V_\parallel|$ [m/s]")
    axes[0, 0].set_title("Raw parallel-flow amplification")
    axes[0, 1].set_xlabel(r"$|B_n|/B$")
    axes[0, 1].set_ylabel(r"$|V_\parallel|/c_s$")
    axes[0, 1].set_title("Parallel Mach number near grazing incidence")
    axes[1, 0].set_xlabel(r"$|v_{E,n}|$ [m/s]")
    axes[1, 0].set_ylabel(r"$|V_\parallel B_n/B|$ [m/s]")
    axes[1, 0].set_title("Normal-flow compensation")
    axes[1, 1].set_xlabel(r"$|B_n|/B$")
    axes[1, 1].set_ylabel(r"$|v_{\mathrm{flow},n}|$ [m/s]")
    axes[1, 1].set_title("Residual normal flow")

    # Equality means the parallel contribution has the same magnitude as ExB.
    all_positive = []
    for _, data in runs:
        all_positive.extend(np.abs(col(data, "vExB_norm")).tolist())
        all_positive.extend(np.abs(col(data, "vparB_norm")).tolist())
    positive = np.asarray(all_positive)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    if positive.size:
        low, high = positive.min(), positive.max()
        axes[1, 0].plot([low, high], [low, high], "k--", linewidth=0.8,
                        label="equal magnitude")

    for axis in axes.flat:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(which="both", alpha=0.25)
        axis.legend(loc="best")
    fig.suptitle("Grazing-incidence Mach diagnostics")
    fig.tight_layout()


def plot_rz_maps(label: str, data: dict[str, np.ndarray]) -> None:
    if "r" not in data or "z" not in data:
        return

    r = col(data, "R")
    z = col(data, "Z")
    bn = np.abs(col(data, "Bnorm") / col(data, "B_abs"))
    quantities = (
        (np.log10(np.maximum(bn, 1.0e-16)), r"$\log_{10}(|B_n|/B)$"),
        (np.log10(np.maximum(np.abs(col(data, "vpar")), 1.0)),
         r"$\log_{10}(|V_\parallel|/[\mathrm{m/s}])$"),
        (col(data, "vExB_norm") / 1.0e3, r"$v_{E,n}$ [km/s]"),
        (col(data, "vflow_norm") / 1.0e3,
         r"$v_{\mathrm{flow},n}$ [km/s]"),
    )

    fig, axes = plt.subplots(2, 2, figsize=(11, 10), sharex=True, sharey=True)
    for axis, (values, title) in zip(axes.flat, quantities):
        points = axis.scatter(r, z, c=values, s=8, cmap="coolwarm")
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(title)
        axis.set_xlabel("R [m]")
        axis.set_ylabel("Z [m]")
        fig.colorbar(points, ax=axis)
    fig.suptitle(f"Boundary location: {label}")
    fig.tight_layout()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path,
                        help="boundary_quantities text files")
    args = parser.parse_args()

    runs = []
    for path in args.files:
        data, time_label = read_boundary_file(path)
        label = path.stem
        if time_label:
            match = re.search(r"#?(\d+)", time_label)
            if match:
                label += f" (step {int(match.group(1))})"
        print_summary(label, data)
        runs.append((label, data))

    plot_profiles(runs)
    plot_grazing_correlations(runs)
    for label, data in runs:
        plot_rz_maps(label, data)
    plt.show()


if __name__ == "__main__":
    main()
