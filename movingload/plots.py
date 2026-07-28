"""
Figures from the extracted CSVs.

    python -m movingload.plots pulse   runs/semi_rigid_vertical
    python -m movingload.plots heatmap runs/semi_rigid_vertical
    python -m movingload.plots compare runs/flexible_vertical runs/semi_rigid_vertical \\
                                       runs/rigid_vertical

Each run directory is expected to hold what `extract_strain.py` produced —
`<case_id>_e33_history.csv` and `<case_id>_e33_centre.csv` — next to the
`<case_id>_mesh.json` sidecar.  Strains are plotted in microstrain
(tension positive, Abaqus convention).

Needs numpy and matplotlib; the rest of the package does not.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_case(run_dir: str | Path) -> dict:
    """Read one run directory: metadata + both CSVs."""
    run_dir = Path(run_dir)
    metas = sorted(run_dir.glob("*_mesh.json"))
    if not metas:
        raise FileNotFoundError(f"no *_mesh.json in {run_dir}")
    meta = json.loads(metas[0].read_text())
    case_id = meta.get("case_id", metas[0].name.replace("_mesh.json", ""))

    centre_csv = run_dir / f"{case_id}_e33_centre.csv"
    if not centre_csv.is_file():
        raise FileNotFoundError(
            f"{centre_csv} missing — run extract_strain.py on the ODB first")
    centre = np.genfromtxt(centre_csv, delimiter=",", names=True)

    hist_csv = run_dir / f"{case_id}_e33_history.csv"
    history = (np.genfromtxt(hist_csv, delimiter=",", names=True)
               if hist_csv.is_file() else None)
    return {"case_id": case_id, "meta": meta, "centre": centre,
            "history": history, "dir": run_dir}


def _annotate_case(ax, meta: dict) -> None:
    t = meta["tyre_geometry"]
    m = meta.get("motion", {})
    bits = [f"R = {t['R']*1e3:.1f} mm", f"p = {t['P0_MPa']:.2f} MPa",
            f"F = {t['force_per_tyre_kN']:.1f} kN"]
    if m.get("mu", 0):
        bits.append(f"mu = {m['mu']:+.2f}" if m.get("mu_sign", 1) > 0
                    else f"mu = -{m['mu']:.2f}")
    ax.text(0.02, 0.02, "   ".join(bits), transform=ax.transAxes,
            fontsize=8, color="0.35", va="bottom")


def plot_pulse(run_dir: str | Path, out: str | Path | None = None) -> Path:
    """Strain-time pulse at the mid-span observation element."""
    c = load_case(run_dir)
    t = c["centre"]["time"]
    e = c["centre"]["E33"] * 1e6
    meta = c["meta"]
    obs = meta["observation"]
    y_mid = 0.5 * (obs["observation_y_top"] + obs["observation_y_bottom"])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, e, "-", lw=2, color="C0")
    ax.axhline(0, color="0.6", lw=0.6)

    i_pk = int(np.argmax(np.abs(e)))
    ax.plot(t[i_pk], e[i_pk], "o", color="C3", ms=5)
    ax.annotate(f"{e[i_pk]:+.1f} " + r"$\mu\varepsilon$" + f"  @ t = {t[i_pk]:.3f}",
                xy=(t[i_pk], e[i_pk]), xytext=(8, 8),
                textcoords="offset points", fontsize=9, color="C3")

    ax.set_xlabel("step time  (tyre centre sweeps the model)")
    ax.set_ylabel(r"$E_{33}$  ($\mu\varepsilon$, tension +)")
    ax.set_title(f"{c['case_id']} — longitudinal strain at the surface-layer "
                 f"bottom\n(x = {obs['observation_x']:.3f} m, "
                 f"y = {y_mid:.3f} m, mid-span)")
    _annotate_case(ax, meta)
    fig.tight_layout()

    out = Path(out) if out else c["dir"] / f"{c['case_id']}_pulse.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def plot_heatmap(run_dir: str | Path, out: str | Path | None = None) -> Path:
    """E33 over (time, longitudinal position) on the whole observation column."""
    c = load_case(run_dir)
    if c["history"] is None:
        raise FileNotFoundError(
            f"{c['case_id']}_e33_history.csv missing — needed for the heatmap")
    h = c["history"]
    times = np.unique(h["time"])
    zs = np.unique(h["z_centre"])
    grid = np.full((len(zs), len(times)), np.nan)
    ti = {v: i for i, v in enumerate(times)}
    zi = {v: i for i, v in enumerate(zs)}
    for row in h:
        grid[zi[row["z_centre"]], ti[row["time"]]] = row["E33"] * 1e6

    fig, ax = plt.subplots(figsize=(8, 4.5))
    lim = np.nanmax(np.abs(grid))
    im = ax.pcolormesh(times, zs, grid, cmap="RdBu_r",
                       vmin=-lim, vmax=lim, shading="nearest")
    fig.colorbar(im, ax=ax, label=r"$E_{33}$  ($\mu\varepsilon$)")

    m = c["meta"].get("motion")
    if m:
        tt = np.linspace(0, m["t_final"], 50)
        zc = m["z_start"] + (m["z_end"] - m["z_start"]) * tt / m["t_final"]
        ax.plot(tt, zc, "k--", lw=1, label="tyre centre")
        ax.legend(loc="upper left", fontsize=8)
        ax.set_ylim(zs.min(), zs.max())

    ax.set_xlabel("step time")
    ax.set_ylabel("longitudinal position  z  (m)")
    ax.set_title(f"{c['case_id']} — strain wave travelling with the tyre")
    fig.tight_layout()

    out = Path(out) if out else c["dir"] / f"{c['case_id']}_heatmap.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def plot_compare(run_dirs: List[str | Path], out: str | Path | None = None,
                 title: str = "") -> Path:
    """Overlay the mid-span pulses of several cases and tabulate the peaks."""
    cases = [load_case(d) for d in run_dirs]

    fig, ax = plt.subplots(figsize=(9, 5))
    rows = []
    for i, c in enumerate(cases):
        t, e = c["centre"]["time"], c["centre"]["E33"] * 1e6
        ax.plot(t, e, "-", lw=2, color=f"C{i}", label=c["case_id"])
        rows.append((c["case_id"], float(e.max()), float(t[np.argmax(e)]),
                     float(e.min()), float(t[np.argmin(e)])))
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("step time")
    ax.set_ylabel(r"$E_{33}$  ($\mu\varepsilon$, tension +)")
    ax.set_title(title or "Longitudinal strain at the surface-layer bottom")
    ax.legend(fontsize=9)
    fig.tight_layout()

    out = Path(out) if out else Path(cases[0]["dir"]).parent / "compare_e33.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)

    width = max(len(r[0]) for r in rows)
    print(f"{'case':<{width}}  {'peak +':>9} {'@t':>6}   {'peak -':>9} {'@t':>6}")
    for name, emax, tmax, emin, tmin in rows:
        print(f"{name:<{width}}  {emax:>9.2f} {tmax:>6.3f}   "
              f"{emin:>9.2f} {tmin:>6.3f}")
    return out


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="python -m movingload.plots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pulse");   p.add_argument("run_dir"); p.add_argument("--out")
    p = sub.add_parser("heatmap"); p.add_argument("run_dir"); p.add_argument("--out")
    p = sub.add_parser("compare")
    p.add_argument("run_dirs", nargs="+")
    p.add_argument("--out")
    p.add_argument("--title", default="")
    args = ap.parse_args(argv)

    if args.cmd == "pulse":
        print("wrote", plot_pulse(args.run_dir, args.out))
    elif args.cmd == "heatmap":
        print("wrote", plot_heatmap(args.run_dir, args.out))
    else:
        print("wrote", plot_compare(args.run_dirs, args.out, args.title))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
