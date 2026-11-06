"""
Driver for the axisymmetric model: build the step-response deck, solve it,
then synthesise the moving dual-tyre load by superposition.

    # build only
    python -m movingload.axi --structure semi_rigid --material viscoelastic --build-only

    # solve here, then superpose at two speeds
    python -m movingload.axi --structure semi_rigid --material elastic --local --speed 7.6

    # solve on the cluster (same environment variables as movingload.submit)
    python -m movingload.axi --structure semi_rigid --material viscoelastic --speed 7.6 --speed 16.7

One Abaqus run serves every speed and tyre layout: the speed enters only in
`movingload.superpose`, which is plain numpy and runs in seconds.  Re-run it
by hand for other speeds:

    python -m movingload.superpose runs/axi_semi_rigid_viscoelastic --speed 20

Each case directory holds

    <case_id>.inp, <case_id>_mesh.json, extract_axi.py, submit.sh | run.sh

and, after the solve, <case_id>_axi_response.csv and one
<case_id>_moving_v<speed>.csv per speed.  There is no dload.f: the print
does not move in the axisymmetric model.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from . import axi_inp, slurm, superpose
from .axi_config import MATERIALS, make_axi_case
from .config import STRUCTURES

_HERE = Path(__file__).resolve().parent
EXTRACTOR = "extract_axi.py"


def prepare_axi(structure: str, material: str, out_root: str | Path = "runs",
                case_id: str | None = None,
                resources: slurm.Resources | None = None,
                local: bool = False,
                abaqus_cmd: str | None = None,
                module_load: str | None = None) -> dict:
    """Materialise one self-contained run directory.  No network, no Abaqus."""
    case = make_axi_case(structure, material, case_id=case_id)
    run_dir = Path(out_root) / case.case_id
    run_dir.mkdir(parents=True, exist_ok=True)
    info = axi_inp.build(case, run_dir)
    shutil.copy(_HERE / EXTRACTOR, run_dir / EXTRACTOR)

    inp_name = Path(info["inp"]).name
    if local:
        script = slurm.write(run_dir / "run.sh", slurm.local_script(
            case.case_id, inp_name, user_sub=None,
            cpus=(resources.cpus if resources else 4),
            abaqus_cmd=abaqus_cmd, module_load=module_load, extractor=EXTRACTOR))
    else:
        script = slurm.write(run_dir / "submit.sh", slurm.slurm_script(
            case.case_id, inp_name, user_sub=None, resources=resources,
            abaqus_cmd=abaqus_cmd, module_load=module_load, extractor=EXTRACTOR))
    return {"config": case, "run_dir": run_dir, "script": script, "inp": info}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m movingload.axi")
    ap.add_argument("--structure", default="semi_rigid", choices=sorted(STRUCTURES))
    ap.add_argument("--material", default="elastic", choices=MATERIALS,
                    help="asphalt layers: elastic, or linear viscoelastic (Prony)")
    ap.add_argument("--case-id", default=None,
                    help="default: axi_<structure>_<material>")
    ap.add_argument("--out", default="runs")
    ap.add_argument("--build-only", action="store_true")
    ap.add_argument("--local", action="store_true")
    ap.add_argument("--speed", type=float, action="append",
                    help="speed(s) to superpose after the solve [m/s] (default 7.6)")
    ap.add_argument("--single", action="store_true",
                    help="superpose one print instead of the dual pair")

    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--mem", default="8g")
    ap.add_argument("--walltime", default="01:00:00")
    ap.add_argument("--partition", default="")
    ap.add_argument("--abaqus-cmd", default=None)
    ap.add_argument("--module-load", default=None)
    ap.add_argument("--no-wait", action="store_true")
    ap.add_argument("--poll", type=int, default=20)
    args = ap.parse_args(argv)

    res = slurm.Resources(cpus=args.cpus, mem=args.mem,
                          walltime=args.walltime, partition=args.partition)
    out = prepare_axi(args.structure, args.material, args.out,
                      case_id=args.case_id, resources=res, local=args.local,
                      abaqus_cmd=args.abaqus_cmd, module_load=args.module_load)
    case, info, run_dir = out["config"], out["inp"], out["run_dir"]
    print(case.summary())
    print(f"\n[build] {info['n_elements']:,} {case.mesh.element} elements, "
          f"{info['n_nodes']:,} nodes; print = {info['n_print']} faces")
    print(f"        observation row at depth {info['obs_depth']:.4f} m")
    print(f"[ready] {run_dir}")
    if args.build_only:
        print("\n--build-only: nothing was run.")
        return 0

    if args.local:
        import subprocess
        rc = subprocess.run(["bash", str(out["script"])], cwd=str(run_dir)).returncode
    else:
        from .submit import _run_cluster
        rc = _run_cluster(case.case_id, run_dir, args.poll, args.no_wait)
        if args.no_wait:
            return rc
    if rc != 0:
        return rc
    try:
        superpose.run(run_dir, args.speed or [7.6], single=args.single)
    except (FileNotFoundError, ValueError) as exc:
        sys.stderr.write(f"superposition failed: {exc}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
