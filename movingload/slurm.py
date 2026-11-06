"""
Emit the shell script that runs one case: solve, then reduce the ODB to CSV.

Two flavours, same body:

    submit.sh   SLURM batch script (`sbatch submit.sh`)
    run.sh      plain shell script for a workstation (`bash run.sh`)

Site specifics — the module-load line, the partition, the Abaqus command —
are all parameters with neutral defaults, because they differ at every
installation.  Set them once via the environment:

    MOVINGLOAD_ABAQUS_MODULE   e.g. "module load abaqus/2023"   (default: none)
    MOVINGLOAD_ABAQUS_CMD      e.g. "abq2023"                   (default: abaqus)
    MOVINGLOAD_PARTITION       e.g. "compute"                   (default: none)
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Resources:
    """What one job asks the scheduler for."""
    cpus: int = 8
    mem: str = "32g"
    walltime: str = "02:00:00"
    partition: str = ""     # empty -> omit the directive, use the site default

    def __post_init__(self):
        if self.cpus < 1:
            raise ValueError("cpus must be >= 1")


_BODY = """set -euo pipefail
cd "{workdir}"

{module_load}
echo "[$(date)] abaqus: $({abaqus} whereami 2>/dev/null || command -v {abaqus} || echo NOT_FOUND)"

echo "[$(date)] solving {job}"
{abaqus} job={job} input={inp}{user_opt} cpus={cpus} interactive
echo "[$(date)] solver finished"

echo "[$(date)] reducing the ODB to CSV"
{abaqus} python {extractor} {job}.odb
echo "[$(date)] done"
"""


def _user_opt(user_sub: str | None) -> str:
    # The axisymmetric deck loads a fixed footprint with *DSLOAD and needs
    # no subroutine; `user=` with nothing behind it would make Abaqus try
    # (and fail) to compile one.
    return f" user={user_sub}" if user_sub else ""


def _resolve(abaqus_cmd: str | None, module_load: str | None) -> tuple[str, str]:
    cmd = abaqus_cmd or os.environ.get("MOVINGLOAD_ABAQUS_CMD", "abaqus")
    mod = module_load if module_load is not None else os.environ.get(
        "MOVINGLOAD_ABAQUS_MODULE", "")
    mod_block = f"{mod}\n" if mod else "# (no module system: abaqus assumed on PATH)\n"
    return cmd, mod_block


def slurm_script(job: str, inp: str, user_sub: str | None = "dload.f",
                 resources: Resources | None = None,
                 abaqus_cmd: str | None = None,
                 module_load: str | None = None,
                 extractor: str = "extract_strain.py") -> str:
    """SLURM batch script for one case.

    `interactive` on the abaqus line is mandatory: without it the command
    returns as soon as the solver is spawned, the batch script ends, and the
    scheduler tears the allocation down mid-solve.
    """
    r = resources or Resources()
    partition = r.partition or os.environ.get("MOVINGLOAD_PARTITION", "")
    cmd, mod_block = _resolve(abaqus_cmd, module_load)

    directives = [
        f"#SBATCH --job-name={job}",
        f"#SBATCH --output={job}_%j.log",
        f"#SBATCH --error={job}_%j.err",
        f"#SBATCH --time={r.walltime}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={r.cpus}",
        f"#SBATCH --mem={r.mem}",
    ]
    if partition:
        directives.insert(3, f"#SBATCH --partition={partition}")

    return ("#!/bin/bash\n" + "\n".join(directives) + "\n\n"
            + _BODY.format(workdir="$SLURM_SUBMIT_DIR", module_load=mod_block,
                           abaqus=cmd, job=job, inp=inp,
                           user_opt=_user_opt(user_sub), cpus=r.cpus,
                           extractor=extractor))


def local_script(job: str, inp: str, user_sub: str | None = "dload.f",
                 cpus: int = 4, abaqus_cmd: str | None = None,
                 module_load: str | None = None,
                 extractor: str = "extract_strain.py") -> str:
    """Plain shell script for running the same case on a workstation."""
    cmd, mod_block = _resolve(abaqus_cmd, module_load)
    return ("#!/bin/bash\n"
            + _BODY.format(workdir="$(dirname \"$0\")", module_load=mod_block,
                           abaqus=cmd, job=job, inp=inp,
                           user_opt=_user_opt(user_sub), cpus=cpus,
                           extractor=extractor))


def write(path: str | Path, text: str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    p.chmod(0o755)
    return p
