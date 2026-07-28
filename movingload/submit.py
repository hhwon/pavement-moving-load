"""
Top-level driver: build a case, then either run it here or ship it to a cluster.

    # just build the deck and look at it
    python -m movingload.submit --structure semi_rigid --motion vertical --build-only

    # build and run with a local Abaqus
    python -m movingload.submit --structure rigid --motion brake --local

    # build, ship to a SLURM cluster, wait, bring the CSVs back
    export MOVINGLOAD_SSH_HOST=user@cluster.example.org
    export MOVINGLOAD_REMOTE_ROOT=/scratch/$USER/movingload
    python -m movingload.submit --structure semi_rigid --motion accel

Each case gets its own directory under `--out` (default `runs/`):

    runs/<case_id>/
        <case_id>.inp           deck
        <case_id>_mesh.json     mesh + observation metadata
        dload.f                 user subroutine for this motion
        extract_strain.py       postprocessor (copied so the run dir is self-contained)
        submit.sh | run.sh      scheduler or shell script
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from . import dload, inp_builder, slurm
from .config import MOTIONS, STRUCTURES, make_case

_HERE = Path(__file__).resolve().parent


def prepare(structure: str, motion: str, out_root: str | Path = "runs",
            case_id: str | None = None,
            resources: slurm.Resources | None = None,
            local: bool = False,
            abaqus_cmd: str | None = None,
            module_load: str | None = None) -> dict:
    """Materialise one self-contained run directory.  No network, no Abaqus."""
    cfg = make_case(structure, motion, case_id=case_id)
    run_dir = Path(out_root) / cfg.case_id
    run_dir.mkdir(parents=True, exist_ok=True)

    info = inp_builder.build(cfg, run_dir)
    sub = dload.build(cfg, run_dir, filename="dload.f")

    # Ship the extractor with the case so the run directory stands alone.
    shutil.copy(_HERE / "extract_strain.py", run_dir / "extract_strain.py")

    job = cfg.case_id
    inp_name = Path(info["inp"]).name
    if local:
        script = slurm.write(run_dir / "run.sh",
                             slurm.local_script(job, inp_name, "dload.f",
                                                cpus=(resources.cpus if resources else 4),
                                                abaqus_cmd=abaqus_cmd,
                                                module_load=module_load))
    else:
        script = slurm.write(run_dir / "submit.sh",
                             slurm.slurm_script(job, inp_name, "dload.f",
                                                resources=resources,
                                                abaqus_cmd=abaqus_cmd,
                                                module_load=module_load))

    return {"config": cfg, "run_dir": run_dir, "script": script,
            "inp": info, "dload": sub}


def _run_local(run_dir: Path, script: Path) -> int:
    import subprocess
    print(f"[run] bash {script}")
    return subprocess.run(["bash", str(script)], cwd=str(run_dir)).returncode


def _run_cluster(case_id: str, run_dir: Path, poll: int, no_wait: bool) -> int:
    from .cluster import Cluster, ClusterNotConfigured
    try:
        cl = Cluster()
        remote = cl.remote_dir_for(case_id)
    except ClusterNotConfigured as exc:
        sys.stderr.write(f"\n{exc}\n")
        return 2

    probe = cl.check()
    if probe["rc"] != 0:
        sys.stderr.write(f"cannot reach {cl.host}: {probe['stderr']}\n")
        return 2
    print(f"[cluster] {probe['stdout'].splitlines()[0] if probe['stdout'] else cl.host}")

    print(f"[push] -> {remote}")
    cl.push(run_dir, remote)
    job = cl.sbatch(remote)
    print(f"[sbatch] job {job.job_id}")

    if no_wait:
        print(f"    later:  python -m movingload.cluster pull {remote} {run_dir}")
        return 0

    print(f"[wait] polling every {poll} s …")
    state = cl.wait(job.job_id, poll_seconds=poll)
    print(f"    final state: {state}")
    print(f"[pull] -> {run_dir}")
    cl.pull(remote, run_dir)
    return 0 if state in ("COMPLETED", "UNKNOWN") else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m movingload.submit")
    ap.add_argument("--structure", default="semi_rigid", choices=sorted(STRUCTURES))
    ap.add_argument("--motion", default="vertical", choices=sorted(MOTIONS))
    ap.add_argument("--case-id", default=None, help="default: <structure>_<motion>")
    ap.add_argument("--out", default="runs", help="parent directory for run dirs")

    ap.add_argument("--build-only", action="store_true",
                    help="write the run directory and stop (no Abaqus, no network)")
    ap.add_argument("--local", action="store_true",
                    help="run Abaqus on this machine instead of a cluster")

    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--mem", default="32g")
    ap.add_argument("--walltime", default="02:00:00")
    ap.add_argument("--partition", default="",
                    help="scheduler partition (default: site default)")
    ap.add_argument("--abaqus-cmd", default=None,
                    help="Abaqus executable name (default: $MOVINGLOAD_ABAQUS_CMD or 'abaqus')")
    ap.add_argument("--module-load", default=None,
                    help="module command to run first, e.g. 'module load abaqus/2023'")

    ap.add_argument("--no-wait", action="store_true",
                    help="submit and return; pull the results later")
    ap.add_argument("--poll", type=int, default=30, help="squeue poll interval (s)")
    args = ap.parse_args(argv)

    res = slurm.Resources(cpus=args.cpus, mem=args.mem,
                          walltime=args.walltime, partition=args.partition)
    out = prepare(args.structure, args.motion, args.out, case_id=args.case_id,
                  resources=res, local=args.local,
                  abaqus_cmd=args.abaqus_cmd, module_load=args.module_load)

    cfg = out["config"]
    print(cfg.summary())
    print(f"\n[build] {out['inp']['n_elements']:,} C3D8R elements, "
          f"{out['inp']['n_nodes']:,} nodes")
    print(f"        observation column at x = {out['inp']['observation_x']:.4f} m")
    print(f"        tyre passes mid-span (z = {cfg.mesh.LZ / 2:g} m) at "
          f"t = {cfg.motion.time_over(cfg.mesh.LZ / 2):.3f}")
    print(f"        user subroutine: {'DLOAD + UTRACLOAD' if out['dload']['has_utracload'] else 'DLOAD only'}")
    print(f"[ready] {out['run_dir']}")

    if args.build_only:
        print("\n--build-only: nothing was run.")
        return 0
    if args.local:
        return _run_local(out["run_dir"], out["script"])
    return _run_cluster(cfg.case_id, out["run_dir"], args.poll, args.no_wait)


if __name__ == "__main__":
    raise SystemExit(main())
