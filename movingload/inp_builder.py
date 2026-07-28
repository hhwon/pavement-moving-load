"""
Write the Abaqus `.inp` deck and the metadata sidecar the extractor reads.

    from movingload.config import make_case
    from movingload.inp_builder import build

    info = build(make_case("semi_rigid", "vertical"), out_dir="runs/demo")

produces `runs/demo/<case_id>.inp` and `runs/demo/<case_id>_mesh.json`.

The deck is deliberately plain: C3D8R elements, `*ELASTIC` isotropic
materials, one `*STATIC, DIRECT` step, and a user-subroutine surface load.
Nothing depends on a particular Abaqus version beyond `DLOAD` / `UTRACLOAD`.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List, Sequence

from .config import CaseConfig
from .mesh import Grid, build_grid


def _chunks(seq: Sequence[int], n: int = 8) -> Iterable[Sequence[int]]:
    for k in range(0, len(seq), n):
        yield seq[k:k + n]


def _id_block(keyword: str, name: str, ids: Sequence[int]) -> List[str]:
    out = [f"**", f"*{keyword}, {keyword}={name}"]
    out.extend(", ".join(str(i) for i in c) for c in _chunks(ids))
    return out


def _nodes(grid: Grid) -> List[str]:
    out = ["**", "*NODE"]
    for k, z in enumerate(grid.z):
        for j, y in enumerate(grid.y):
            for i, x in enumerate(grid.x):
                out.append(f"{grid.node_id(i, j, k)}, "
                           f"{x:.6f}, {y:.6f}, {z:.6f}")
    return out


def _elements(grid: Grid) -> List[str]:
    out = ["**", "*ELEMENT, TYPE=C3D8R, ELSET=ALL_ELEMS"]
    for k in range(grid.NZ):
        for j in range(grid.NY):
            for i in range(grid.NX):
                n = [grid.node_id(i,     j,     k),
                     grid.node_id(i + 1, j,     k),
                     grid.node_id(i + 1, j + 1, k),
                     grid.node_id(i,     j + 1, k),
                     grid.node_id(i,     j,     k + 1),
                     grid.node_id(i + 1, j,     k + 1),
                     grid.node_id(i + 1, j + 1, k + 1),
                     grid.node_id(i,     j + 1, k + 1)]
                out.append(f"{grid.elem_id(i, j, k)}, "
                           + ", ".join(str(v) for v in n))
    return out


def _layer_elsets(cfg: CaseConfig, grid: Grid) -> List[str]:
    out: List[str] = []
    for ilayer, name in enumerate(cfg.structure.names):
        ids = [grid.elem_id(i, j, k)
               for k in range(grid.NZ)
               for j in range(grid.NY)
               for i in range(grid.NX)
               if grid.layer_of_row[j] == ilayer]
        out.extend(_id_block("ELSET", name, ids))
    return out


def _loaded_and_observed(grid: Grid) -> List[str]:
    top = [grid.elem_id(i, 0, k) for k in range(grid.NZ) for i in range(grid.NX)]
    out = _id_block("ELSET", "TOP_ELEMS", top)
    out += _id_block("ELSET", "OBS_ELEMS", grid.obs_elem_ids)
    return out


def _boundary_nsets(grid: Grid) -> List[str]:
    nx, ny, nz = grid.NX, grid.NY, grid.NZ
    sets = {
        "SYM_X":  [grid.node_id(0,  j, k) for k in range(nz + 1) for j in range(ny + 1)],
        "FAR_X":  [grid.node_id(nx, j, k) for k in range(nz + 1) for j in range(ny + 1)],
        "FAR_Y":  [grid.node_id(i, ny, k) for k in range(nz + 1) for i in range(nx + 1)],
        "FAR_Z0": [grid.node_id(i, j,  0) for j in range(ny + 1) for i in range(nx + 1)],
        "FAR_ZL": [grid.node_id(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)],
    }
    out: List[str] = []
    for name, ids in sets.items():
        out.extend(_id_block("NSET", name, ids))
    return out


def _materials(cfg: CaseConfig) -> List[str]:
    out: List[str] = []
    for name, E, nu in zip(cfg.structure.names, cfg.structure.E, cfg.structure.nu):
        out += ["**", f"*SOLID SECTION, ELSET={name}, MATERIAL={name}_MAT",
                "**", f"*MATERIAL, NAME={name}_MAT",
                "*ELASTIC", f"{E:.4e}, {nu}"]
    return out


def _boundary(cfg: CaseConfig) -> List[str]:
    return [
        "**", "*BOUNDARY",
        "SYM_X,  1, 1",     # symmetry plane: midline between the dual tyres
        "FAR_X,  1, 1",     # lateral roller
        "FAR_Y,  1, 3",     # fixed base
        "FAR_Z0, 3, 3",     # longitudinal rollers at both ends
        "FAR_ZL, 3, 3",
    ]


def _step(cfg: CaseConfig) -> List[str]:
    m, s = cfg.motion, cfg.step
    out = [
        "**", f"*STEP, NLGEOM=NO, INC={s.max_inc}",
        f"Moving tyre traversing z = {m.z_start} -> {m.z_end} m "
        f"in step time 0 -> {m.t_final}  ({cfg.case_id})",
        "*STATIC, DIRECT",
        f"{s.dt}, {m.t_final}",
        "**", "*DLOAD, OP=NEW",
        "TOP_ELEMS, P3NU, 1.0",
    ]
    if m.horizontal:
        # Tangential traction on the same top faces; magnitude and direction
        # both come from UTRACLOAD, so the .inp only names the load type.
        out += ["**", "*DSLOAD, OP=NEW",
                "TOP_ELEMS, TRVEC3NU, 1.0"]
    out += [
        "**", f"*OUTPUT, FIELD, FREQUENCY={s.field_frequency}",
        "*NODE OUTPUT", "U,",
        "*ELEMENT OUTPUT, POSITION=INTEGRATION POINTS", "S, E",
        "**", f"*OUTPUT, HISTORY, FREQUENCY={s.history_frequency}",
        "*ELEMENT OUTPUT, ELSET=OBS_ELEMS", "E33,",
        "*END STEP",
    ]
    return out


def deck_lines(cfg: CaseConfig, grid: Grid) -> List[str]:
    """The complete `.inp` as a list of lines (no trailing newlines)."""
    out = [
        "*HEADING",
        f"Moving single-tyre load on a 3-D half pavement — {cfg.case_id}",
        "** generated by movingload.inp_builder",
        "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO",
    ]
    out += _nodes(grid)
    out += _elements(grid)
    out += _layer_elsets(cfg, grid)
    out += _loaded_and_observed(grid)
    out += _boundary_nsets(grid)
    out += _materials(cfg)
    out += _boundary(cfg)
    out += _step(cfg)
    return out


def metadata(cfg: CaseConfig, grid: Grid) -> dict:
    """Sidecar consumed by `extract_strain.py` and the plotting helpers.

    Everything the postprocessor needs to interpret the ODB without
    re-deriving the mesh: node coordinates, the observation element ids and
    where they sit, and the tyre / motion parameters used for the load.
    """
    t, m = cfg.tyre, cfg.motion
    return {
        "case_id": cfg.case_id,
        "structure": cfg.structure.name,
        "X_NODES": grid.x.tolist(),
        "Y_NODES": grid.y.tolist(),
        "Z_NODES": grid.z.tolist(),
        "NX": grid.NX, "NY": grid.NY, "NZ": grid.NZ,
        "H_LAYERS": cfg.structure.thicknesses(cfg.mesh.LY),
        "E_LAYERS": list(cfg.structure.E),
        "NU_LAYERS": list(cfg.structure.nu),
        "NEL_PER_LAYER": list(cfg.structure.n_elems),
        "NAMES": list(cfg.structure.names),
        "tyre_geometry": {
            "X_TYRE": t.x_centre,
            "R": t.radius,
            "P0_MPa": t.pressure / 1e6,
            "force_per_tyre_kN": t.force / 1e3,
            "dual_tyre_spacing_m": t.dual_spacing,
        },
        "motion": {
            "z_start": m.z_start, "z_end": m.z_end, "t_final": m.t_final,
            "mu": m.mu, "mu_sign": m.mu_sign,
        },
        "observation": {
            "i_obs": grid.i_obs,
            "j_obs": grid.j_obs,
            "OBS_ELEM_IDS": grid.obs_elem_ids,
            "observation_x": grid.obs_x,
            "observation_y_top": float(grid.y[grid.j_obs]),
            "observation_y_bottom": float(grid.y[grid.j_obs + 1]),
        },
    }


def build(cfg: CaseConfig, out_dir: str | Path) -> dict:
    """Write `<case_id>.inp` and `<case_id>_mesh.json` into `out_dir`.

    Returns a small info dict (paths and mesh size).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grid = build_grid(cfg)

    inp = out_dir / f"{cfg.case_id}.inp"
    inp.write_text("\n".join(deck_lines(cfg, grid)) + "\n", encoding="utf-8")

    meta = out_dir / f"{cfg.case_id}_mesh.json"
    meta.write_text(json.dumps(metadata(cfg, grid), indent=2), encoding="utf-8")

    return {
        "case_id": cfg.case_id,
        "inp": str(inp),
        "mesh_json": str(meta),
        "NX": grid.NX, "NY": grid.NY, "NZ": grid.NZ,
        "n_elements": grid.n_elements,
        "n_nodes": grid.n_nodes,
        "observation_x": grid.obs_x,
        "horizontal_traction": cfg.motion.horizontal,
    }
