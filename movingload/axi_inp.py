"""
Write the axisymmetric step-response deck and its metadata sidecar.

    from movingload.axi_config import make_axi_case
    from movingload.axi_inp import build
    build(make_axi_case("semi_rigid", "viscoelastic"), "runs/axi_demo")

The deck applies the tyre pressure on a *stationary* print centred on the
axis and records the strain along the bottom row of the surface layer:

    step INSTANT    *STATIC          load switched on; the instantaneous
                                     (glassy, for viscoelastic layers) response
    steps CREEP_nn  *VISCO           load held; creep of the asphalt layers
                                     (viscoelastic cases only)

`*VISCOELASTIC` is inert in a *STATIC step — Abaqus uses the instantaneous
moduli there — which is exactly what makes INSTANT the s = 0+ response.
The moving load is synthesised from this afterwards by `movingload.superpose`.

No user subroutine: the print does not move, so a plain *DSLOAD does.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List

from .axi_config import AxiCase
from .axi_mesh import AxiGrid, axi_structure, build_axi_grid
from .inp_builder import _id_block

#: Period of the INSTANT step.  Its frames are mapped to s = 0 by the
#: extractor, so the value only has to be known, not meaningful.
STATIC_PERIOD = 1.0


def _nodes(grid: AxiGrid) -> List[str]:
    out = ["**", "*NODE"]
    for j, d in enumerate(grid.depth):
        for i, r in enumerate(grid.r):
            out.append(f"{grid.node_id(i, j)}, {r:.6f}, {-d:.6f}")
    return out


def _elements(case: AxiCase, grid: AxiGrid) -> List[str]:
    out = ["**", f"*ELEMENT, TYPE={case.mesh.element}, ELSET=ALL_ELEMS"]
    for j in range(grid.NY):
        for i in range(grid.NR):
            out.append(f"{grid.elem_id(i, j)}, "
                       + ", ".join(str(n) for n in grid.connectivity(i, j)))
    return out


def _sets(case: AxiCase, grid: AxiGrid) -> List[str]:
    out: List[str] = []
    for ilayer, name in enumerate(case.structure.names):
        ids = [grid.elem_id(i, j) for j in range(grid.NY) for i in range(grid.NR)
               if grid.layer_of_row[j] == ilayer]
        out += _id_block("ELSET", name, ids)
    out += _id_block("ELSET", "PRINT_ELEMS",
                     [grid.elem_id(i, 0) for i in range(grid.n_print)])
    out += _id_block("ELSET", "OBS_ROW", grid.obs_elem_ids)
    out += ["**", "*SURFACE, TYPE=ELEMENT, NAME=TYRE_PRINT", "PRINT_ELEMS, S3"]

    nr, ny = grid.NR, grid.NY
    out += _id_block("NSET", "AXIS", [grid.node_id(0, j) for j in range(ny + 1)])
    out += _id_block("NSET", "FAR_R", [grid.node_id(nr, j) for j in range(ny + 1)])
    out += _id_block("NSET", "BASE", [grid.node_id(i, ny) for i in range(nr + 1)])
    return out


def _materials(case: AxiCase) -> List[str]:
    out: List[str] = []
    s = case.structure
    for name, E, nu in zip(s.names, s.E, s.nu):
        # No data line: for CAX elements the section is the full revolution.
        out += ["**", f"*SOLID SECTION, ELSET={name}, MATERIAL={name}_MAT",
                "**", f"*MATERIAL, NAME={name}_MAT"]
        ve = case.visco.get(name)
        if ve is None:
            out += ["*ELASTIC", f"{E:.4e}, {nu}"]
        else:
            # MODULI must be explicit: with Prony data Abaqus defaults to
            # LONG TERM and would read E0 as E_inf.
            out += ["*ELASTIC, MODULI=INSTANTANEOUS", f"{ve.E0:.6e}, {ve.nu}",
                    "*VISCOELASTIC, TIME=PRONY"]
            out += [f"{g:.6e}, {g:.6e}, {tau:.6e}" for g, tau in zip(ve.g, ve.tau)]
    return out


def _boundary() -> List[str]:
    return ["**", "*BOUNDARY",
            "AXIS,  1, 1",      # u_r = 0 on the axis of symmetry
            "FAR_R, 1, 1",      # radial roller at R_out
            "BASE,  1, 2"]      # fixed base


def _steps(case: AxiCase) -> List[str]:
    p = case.tyre.pressure
    out = ["**", "*STEP, NAME=INSTANT, NLGEOM=NO, INC=10",
           "Tyre print switched on: instantaneous response",
           "*STATIC", f"{STATIC_PERIOD}, {STATIC_PERIOD}",
           "**", "*DSLOAD, OP=NEW", f"TYRE_PRINT, P, {p:.6e}",
           # Observation row every increment, at the element centroid (the
           # superposition samples one value per radius).  The whole-model
           # request writes only the last increment of each step and carries
           # no E, so the E field in the ODB is the observation row alone.
           "**", "*OUTPUT, FIELD, FREQUENCY=1",
           "*ELEMENT OUTPUT, ELSET=OBS_ROW, POSITION=CENTROIDAL", "E,",
           "**", "*OUTPUT, FIELD, FREQUENCY=100000",
           "*NODE OUTPUT", "U,",
           "*ELEMENT OUTPUT", "S,",
           "*END STEP"]
    if not case.viscoelastic:
        return out
    for name, dt, period in case.creep.steps():
        n = int(round(period / dt))
        out += ["**", f"*STEP, NAME={name}, NLGEOM=NO, INC={n + 10}",
                f"Print held: creep over {period:.6g} s",
                # No CETOL: fixed increments, so every output time is known
                # in advance (CreepSteps.output_times).
                "*VISCO", f"{dt:.9e}, {period:.9e}",
                "*END STEP"]
    return out


def deck_lines(case: AxiCase, grid: AxiGrid) -> List[str]:
    out = ["*HEADING",
           f"Axisymmetric step response of a stationary tyre print - {case.case_id}",
           "** generated by movingload.axi_inp",
           "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO"]
    out += _nodes(grid)
    out += _elements(case, grid)
    out += _sets(case, grid)
    out += _materials(case)
    out += _boundary()
    out += _steps(case)
    return out


def metadata(case: AxiCase, grid: AxiGrid) -> dict:
    """Sidecar for `extract_axi.py` and `movingload.superpose`."""
    s = axi_structure(case.structure, case.mesh)
    t = case.tyre
    return {
        "case_id": case.case_id,
        "model": "axisymmetric",
        "structure": case.structure.name,
        "material": case.material,
        "element": case.mesh.element,
        "R_NODES": grid.r.tolist(),
        "DEPTH_NODES": grid.depth.tolist(),
        "NR": grid.NR, "NY": grid.NY,
        "H_LAYERS": case.structure.thicknesses(case.mesh.LY),
        "E_LAYERS": list(case.structure.E),
        "NU_LAYERS": list(case.structure.nu),
        "NEL_PER_LAYER": list(s.n_elems),
        "NAMES": list(case.structure.names),
        "viscoelastic_layers": {
            name: {"E0": ve.E0, "nu": ve.nu, "g": list(ve.g), "tau": list(ve.tau),
                   "E_inf": ve.E_inf}
            for name, ve in case.visco.items()},
        "tyre_geometry": {
            "X_TYRE": t.x_centre, "R": t.radius, "P0_MPa": t.pressure / 1e6,
            "force_per_tyre_kN": t.force / 1e3,
            "dual_tyre_spacing_m": t.dual_spacing,
            "load_face": "S3", "n_print_faces": grid.n_print,
        },
        "time": {
            "static_step": "INSTANT",
            "static_period": STATIC_PERIOD,
            "creep_steps": [list(x) for x in case.creep.steps()]
                           if case.viscoelastic else [],
            "s_max": case.creep.t_max if case.viscoelastic else 0.0,
            # Every s the extractor must find, in double precision.  The ODB
            # stores step and frame times in single precision, so the
            # extractor snaps to these instead of trusting its own sums.
            "output_times": case.creep.output_times() if case.viscoelastic else [0.0],
        },
        "observation": {
            "j_obs": grid.j_obs,
            "OBS_ELEM_IDS": grid.obs_elem_ids,
            "r_centres": grid.r_centres,
            "depth": grid.obs_depth,
            "components": {"E11": "radial", "E22": "axial (vertical)",
                           "E33": "hoop", "E12": "shear r-y (engineering)"},
        },
    }


def build(case: AxiCase, out_dir: str | Path) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grid = build_axi_grid(case)
    inp = out_dir / f"{case.case_id}.inp"
    inp.write_text("\n".join(deck_lines(case, grid)) + "\n", encoding="utf-8")
    meta = out_dir / f"{case.case_id}_mesh.json"
    meta.write_text(json.dumps(metadata(case, grid), indent=2), encoding="utf-8")
    return {"case_id": case.case_id, "inp": str(inp), "mesh_json": str(meta),
            "NR": grid.NR, "NY": grid.NY, "n_elements": grid.n_elements,
            "n_nodes": grid.n_nodes, "obs_depth": grid.obs_depth,
            "n_print": grid.n_print}
