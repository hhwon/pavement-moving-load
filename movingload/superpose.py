"""
Moving dual-tyre load from the axisymmetric step response, by superposition.

    python -m movingload.superpose runs/axi_semi_rigid_viscoelastic --speed 7.6 --speed 16.7

Input: `<case>_axi_response.csv` from `extract_axi.py` — the strain R(rho, s)
at the bottom of the surface layer, at radius rho from a stationary print
that was switched on at s = 0 and held.

Why this is exact (not an approximation of the moving load)
------------------------------------------------------------
The model is linear (elastic or linear viscoelastic), quasi-static, and its
layers are homogeneous in plan, so

  * translating a load translates its response (spatial invariance), and
  * the response to a sum of loads is the sum of the responses, each with
    its own time history (Boltzmann superposition).

Let the tyre centre advance in steps of ds = v*dt: during [t_k, t_k+1) it
sits at z_k.  That is a sum of prints switched on at t_k and off at t_k+1,
so at the observation point P, at t_n,

    eps(t_n) = sum_{k<=n} R(P - z_k, t_n - t_k)  -  sum_{k<n} R(P - z_k, t_n - t_k+1)

which, writing G[m, k] = R(P - z_k, m*dt), is a difference of
anti-diagonal sums:  eps_n = D_n - D_{n-1},  D_n = sum_{m+k=n} G[m, k].
The 3-D model's DLOAD makes the same discretisation (the print jumps once
per increment).  For an elastic layer R does not depend on s, the sum
telescopes, and eps_n = R(P - z_n): the familiar "stationary field read as
a moving pulse".

Tensor rotation and the second tyre
-----------------------------------
R gives radial (E11), hoop (E33) and vertical (E22) strain.  For a print
whose centre is offset (dx, dz) from P, with rho = |(dx, dz)| and
c = dz/rho, s = dx/rho:

    eps_zz (longitudinal) = c^2 E11 + s^2 E33
    eps_xx (transverse)   = s^2 E11 + c^2 E33
    eps_yy (vertical)     = E22

Each tyre of the dual pair is rotated this way and the results are added.
The output uses the 3-D model's global names — E11 transverse, E22 vertical,
E33 longitudinal — so it lines up with `<case>_e33_centre.csv`.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np


@dataclass(frozen=True)
class Response:
    """R(rho, s): one array per component, shape (n_s, n_r)."""
    s: np.ndarray
    r: np.ndarray
    E11: np.ndarray      # radial
    E22: np.ndarray      # vertical
    E33: np.ndarray      # hoop

    def __post_init__(self):
        n = (len(self.s), len(self.r))
        for name in ("E11", "E22", "E33"):
            if getattr(self, name).shape != n:
                raise ValueError(f"{name} has shape {getattr(self, name).shape}, "
                                 f"expected {n}")
        if len(self.s) == 0 or self.s[0] != 0.0:
            raise ValueError("the response must start at s = 0 (the INSTANT step)")
        if np.any(np.diff(self.s) <= 0) or np.any(np.diff(self.r) <= 0):
            raise ValueError("s and r must be strictly increasing")

    @property
    def elastic(self) -> bool:
        return len(self.s) == 1


def read_response(run_dir: str | Path) -> tuple[Response, dict]:
    """Load `<case>_axi_response.csv` and its metadata."""
    run_dir = Path(run_dir)
    metas = sorted(run_dir.glob("*_mesh.json"))
    if not metas:
        raise FileNotFoundError(f"no *_mesh.json in {run_dir}")
    meta = json.loads(metas[0].read_text())
    if meta.get("model") != "axisymmetric":
        raise ValueError(f"{metas[0]} is not an axisymmetric case")
    path = run_dir / f"{meta['case_id']}_axi_response.csv"
    if not path.is_file():
        raise FileNotFoundError(f"{path} missing — run extract_axi.py on the ODB first")

    data = np.genfromtxt(path, delimiter=",", names=True)
    s = np.unique(data["s"])
    r = np.unique(data["r"])
    if len(data) != len(s) * len(r):
        raise ValueError(f"{path}: {len(data)} rows is not {len(s)} times x "
                         f"{len(r)} radii — some frames are incomplete")
    expected = meta.get("time", {}).get("output_times")
    if expected is not None and (len(s) != len(expected)
                                 or not np.allclose(s, expected, rtol=1e-6, atol=1e-12)):
        raise ValueError(
            f"{path}: {len(s)} distinct times, the deck defines {len(expected)}; "
            f"re-run extract_axi.py (frames missing, duplicated or mistimed)")
    order = np.lexsort((data["r"], data["s"]))
    comp = {c: data[c][order].reshape(len(s), len(r)) for c in ("E11", "E22", "E33")}
    return Response(s=s, r=r, **comp), meta


def _at_lags(resp: Response, lags: np.ndarray) -> Dict[str, np.ndarray]:
    """Each component linearly interpolated in s at the given lags:
    arrays of shape (len(lags), n_r)."""
    if resp.elastic:
        return {c: np.repeat(getattr(resp, c), len(lags), axis=0)
                for c in ("E11", "E22", "E33")}
    if lags[-1] > resp.s[-1] * (1 + 1e-9):
        raise ValueError(
            f"the passage lasts {lags[-1]:.4g} s but the creep analysis only "
            f"reached s = {resp.s[-1]:.4g} s; raise CreepSteps.t_max or the speed")
    i = np.clip(np.searchsorted(resp.s, lags, side="right") - 1, 0, len(resp.s) - 2)
    w = ((lags - resp.s[i]) / (resp.s[i + 1] - resp.s[i]))[:, None]
    return {c: (1 - w) * getattr(resp, c)[i] + w * getattr(resp, c)[i + 1]
            for c in ("E11", "E22", "E33")}


def _radial_weights(r: np.ndarray, rho: np.ndarray):
    """Index/weight pairs for linear interpolation in rho.  Inside the first
    centroid the normal strains are flat (they are even in r), so clamp."""
    if rho.max() > r[-1]:
        raise ValueError(
            f"the tyre reaches rho = {rho.max():.3f} m from the observation point "
            f"but the response stops at r = {r[-1]:.3f} m; shorten the travel or "
            f"raise AxiMesh.R_out")
    rr = np.maximum(rho, r[0])
    i = np.clip(np.searchsorted(r, rr, side="right") - 1, 0, len(r) - 2)
    w = (rr - r[i]) / (r[i + 1] - r[i])
    return i, w


def moving_history(resp: Response, speed: float,
                   tyre_x: Sequence[float] = (0.16, -0.16),
                   x_obs: float = 0.16,
                   half_length: float = 3.5,
                   ds: float = 0.005) -> Dict[str, np.ndarray]:
    """Strain history at P = (x_obs, z = 0) while the tyres travel along +z
    from z = -half_length to +half_length at constant `speed` [m/s].

    `tyre_x` lists the tyre centres across the lane: the default is the
    BZZ-100 dual pair (±0.16 m) with P under one of them, the same layout
    the 3-D model reaches through its symmetry plane.  Pass one value for a
    single print.

    Returns t, z_rel (tyre centre minus P, along z) and E11 / E22 / E33 in
    the 3-D global convention (transverse / vertical / longitudinal).
    """
    if speed <= 0 or ds <= 0 or half_length <= 0:
        raise ValueError("speed, ds and half_length must be > 0")
    n = int(round(2 * half_length / ds)) + 1
    z = -half_length + ds * np.arange(n)
    dt = ds / speed
    lags = dt * np.arange(n)
    at = _at_lags(resp, lags)                       # (n_lag, n_r) each

    G = {c: np.zeros((n, n)) for c in ("E11", "E22", "E33")}   # [lag, k]
    for xt in tyre_x:
        dx = x_obs - xt
        dz = -z                                     # P minus tyre centre
        rho = np.hypot(dx, dz)
        safe = np.where(rho > 0, rho, 1.0)
        c2 = np.where(rho > 0, (dz / safe) ** 2, 1.0)
        s2 = 1.0 - c2
        i, w = _radial_weights(resp.r, rho)
        v = {c: (1 - w) * at[c][:, i] + w * at[c][:, i + 1] for c in at}
        G["E33"] += c2 * v["E11"] + s2 * v["E33"]   # longitudinal
        G["E11"] += s2 * v["E11"] + c2 * v["E33"]   # transverse
        G["E22"] += v["E22"]

    diag = np.add.outer(np.arange(n), np.arange(n)).ravel()
    out = {"t": lags, "z_rel": z}
    for c, g in G.items():
        D = np.bincount(diag, weights=g.ravel(), minlength=2 * n)[:n]
        out[c] = D - np.concatenate(([0.0], D[:-1]))
    return out


def peaks(hist: Dict[str, np.ndarray]) -> dict:
    e = hist["E33"]
    k0 = int(np.argmin(np.abs(hist["z_rel"])))
    lo, hi = int(np.argmin(e)), int(np.argmax(e))
    return {"overhead": (float(hist["t"][k0]), float(e[k0])),
            "min": (float(hist["t"][lo]), float(e[lo]), float(hist["z_rel"][lo])),
            "max": (float(hist["t"][hi]), float(e[hi]), float(hist["z_rel"][hi]))}


def write_history(path: str | Path, hist: Dict[str, np.ndarray]) -> Path:
    path = Path(path)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time", "z_rel", "E11", "E22", "E33"])
        for row in zip(hist["t"], hist["z_rel"], hist["E11"], hist["E22"], hist["E33"]):
            w.writerow([f"{x:.9g}" for x in row])
    return path


def history_name(case_id: str, speed: float, single: bool = False) -> str:
    return f"{case_id}_moving_v{speed:g}{'_single' if single else ''}.csv"


def run(run_dir: str | Path, speeds: Sequence[float], ds: float = 0.005,
        half_length: float = 3.5, single: bool = False,
        x_obs: float | None = None) -> List[Path]:
    """Synthesise and write one history per speed; print the peak table."""
    resp, meta = read_response(run_dir)
    xc = meta["tyre_geometry"]["X_TYRE"]
    tyres = (xc,) if single else (xc, -xc)
    x_obs = xc if x_obs is None else x_obs
    print(f"{meta['case_id']}: {meta['material']}, response at {len(resp.s)} "
          f"times x {len(resp.r)} radii, depth {meta['observation']['depth']:.4f} m")
    print(f"  tyres at x = {', '.join(f'{x:+.3f}' for x in tyres)} m, "
          f"observation at x = {x_obs:+.3f} m")
    print(f"  {'v [m/s]':>8s} {'overhead':>10s} {'min E33':>18s} {'max E33':>18s}   (µε)")
    written = []
    for v in speeds:
        h = moving_history(resp, v, tyre_x=tyres, x_obs=x_obs,
                           half_length=half_length, ds=ds)
        p = peaks(h)
        print(f"  {v:8.3g} {p['overhead'][1]*1e6:10.2f} "
              f"{p['min'][1]*1e6:8.2f} @z={p['min'][2]:+.3f} "
              f"{p['max'][1]*1e6:8.2f} @z={p['max'][2]:+.3f}")
        written.append(write_history(Path(run_dir) / history_name(meta["case_id"], v, single), h))
    return written


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m movingload.superpose")
    ap.add_argument("run_dir")
    ap.add_argument("--speed", type=float, action="append",
                    help="tyre speed [m/s]; repeat for several (default 7.6, "
                         "the 3-D reference: 7.6 m in 1 s)")
    ap.add_argument("--ds", type=float, default=0.005,
                    help="distance the tyre advances per load step [m]")
    ap.add_argument("--half-length", type=float, default=3.5,
                    help="travel from -L to +L about the observation point [m]")
    ap.add_argument("--single", action="store_true",
                    help="one print instead of the dual pair")
    ap.add_argument("--x-obs", type=float, default=None,
                    help="transverse position of the observation point "
                         "(default: under the tyre centre)")
    args = ap.parse_args(argv)
    for p in run(args.run_dir, args.speed or [7.6], ds=args.ds,
                 half_length=args.half_length, single=args.single, x_obs=args.x_obs):
        print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
