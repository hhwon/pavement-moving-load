"""
Optional SLURM cluster transport: rsync the run directory over, sbatch it,
wait, rsync the results back.

**Nothing is hard-coded.** There is no built-in host, user, port or remote
path — the target must be configured, or every call raises. Configure it
either through the environment:

    export MOVINGLOAD_SSH_HOST=user@cluster.example.org
    export MOVINGLOAD_SSH_PORT=22                 # optional, default 22
    export MOVINGLOAD_REMOTE_ROOT=/scratch/$USER/movingload

or in code:

    from movingload.cluster import Cluster
    cl = Cluster(host="user@cluster.example.org", remote_root="/scratch/me/ml")

Only `ssh` and `rsync` are required locally; Abaqus runs on the cluster.
If you run Abaqus on your own machine instead, skip this module entirely —
`movingload.submit --local` never imports it.
"""
from __future__ import annotations

import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

#: What `pull` brings back by default.  The ODB is deliberately excluded:
#: it is typically hundreds of MB, and `extract_strain.py` has already
#: reduced it to CSV on the cluster.
RESULT_PATTERNS = ("*.csv", "*.json", "*.log", "*.err", "*.sta", "*.msg", "*.dat")

#: Output a previous run may have left behind, cleared before a fresh push.
#: Inputs — the deck, `dload.f`, the scripts, the mesh sidecar — are not here.
#: `fort.*` catches whatever a user subroutine wrote to a bare Fortran unit.
STALE_GLOBS = ("*.csv", "*.odb", "*.dat", "*.msg", "*.sta", "*.fil", "*.log",
               "*.err", "*.lck", "*.prt", "*.sim", "*.com", "*.res", "*.stt",
               "*.mdl", "*.abq", "*.pac", "fort.*")

_JOB_ID_RE = re.compile(r"Submitted batch job (\d+)")


class ClusterNotConfigured(RuntimeError):
    """Raised when no SSH target has been supplied."""


@dataclass
class SubmitInfo:
    job_id: str
    remote_dir: str
    stdout: str


@dataclass
class Cluster:
    """An SSH/SLURM target.

    Parameters default to the corresponding environment variables; `host`
    has no fallback and must be supplied one way or the other.
    """
    host: Optional[str] = None
    port: str = ""
    remote_root: str = ""

    def __post_init__(self):
        self.host = self.host or os.environ.get("MOVINGLOAD_SSH_HOST") or None
        self.port = self.port or os.environ.get("MOVINGLOAD_SSH_PORT", "22")
        self.remote_root = (self.remote_root
                            or os.environ.get("MOVINGLOAD_REMOTE_ROOT", ""))

    # -- plumbing -------------------------------------------------------
    def _require_host(self) -> str:
        if not self.host:
            raise ClusterNotConfigured(
                "No SSH target. Set MOVINGLOAD_SSH_HOST=user@host (and "
                "optionally MOVINGLOAD_SSH_PORT, MOVINGLOAD_REMOTE_ROOT), or "
                "pass Cluster(host=...). To run Abaqus locally instead, use "
                "`python -m movingload.submit --local`.")
        return self.host

    def _ssh(self, remote_cmd: str) -> list[str]:
        return ["ssh", "-p", self.port,
                "-o", "BatchMode=yes",
                "-o", "StrictHostKeyChecking=accept-new",
                self._require_host(), remote_cmd]

    def _rsync(self, src: str, dst: str,
               includes: Optional[Sequence[str]] = None,
               excludes: Optional[Sequence[str]] = None) -> list[str]:
        cmd = ["rsync", "-av", "-e", f"ssh -p {self.port}"]
        for inc in includes or ():
            cmd += ["--include", inc]
        for exc in excludes or ():
            cmd += ["--exclude", exc]
        return cmd + [src, dst]

    @staticmethod
    def _expand(path: str) -> str:
        """Let `~` expand on the remote shell (shlex.quote would suppress it)."""
        if path == "~":
            return "$HOME"
        if path.startswith("~/"):
            return "$HOME/" + path[2:]
        return path

    # -- operations -----------------------------------------------------
    def run(self, remote_cmd: str, check: bool = True) -> subprocess.CompletedProcess:
        """Run a command on the cluster over ssh."""
        return subprocess.run(self._ssh(remote_cmd), check=check,
                              capture_output=True, text=True)

    def check(self) -> dict:
        """Connectivity probe: hostname, sbatch, abaqus."""
        cp = self.run("hostname; which sbatch || echo NO_SBATCH; "
                      "which abaqus || echo NO_ABAQUS_IN_LOGIN_SHELL",
                      check=False)
        return {"rc": cp.returncode,
                "stdout": cp.stdout.strip(),
                "stderr": cp.stderr.strip()}

    def remote_dir_for(self, case_id: str) -> str:
        if not self.remote_root:
            raise ClusterNotConfigured(
                "No remote root. Set MOVINGLOAD_REMOTE_ROOT=/path/on/cluster "
                "or pass Cluster(remote_root=...).")
        return f"{self.remote_root.rstrip('/')}/{case_id}"

    def mkdir(self, remote_dir: str) -> None:
        self.run(f"mkdir -p {self._expand(remote_dir)}", check=True)

    def clear_results(self, remote_dir: str) -> bool:
        """Delete a previous run's output from `remote_dir`, keeping the inputs.

        rsync is run without `--delete`, so re-using a run directory leaves the
        last run's CSVs in place.  A job that then dies *before* writing its own
        would let `pull` bring the stale ones home, where they are
        indistinguishable from fresh results.  Clearing first turns that silent
        wrong answer into an obvious missing file.

        This issues `rm -f` on the cluster, so it refuses any path that is
        shallow or contains shell metacharacters.  Returns whether it ran.
        """
        d = remote_dir.rstrip("/")
        depth = len([p for p in d.replace("~", "").split("/") if p])
        if depth < 2 or any(c in d for c in " ;&|$`\n*?()<>'\"\\"):
            return False
        pattern = " ".join(f"{self._expand(d)}/{g}" for g in STALE_GLOBS)
        self.run(f"rm -f {pattern}", check=False)
        return True

    def push(self, local_dir: str | Path, remote_dir: str,
             clear_stale: bool = True) -> None:
        """Mirror `local_dir/` to `host:remote_dir/`.

        By default the previous run's output is removed first — see
        `clear_results`.  Pass `clear_stale=False` to add to a directory whose
        existing results you want to keep.
        """
        local = str(Path(local_dir).resolve())
        self.mkdir(remote_dir)
        if clear_stale:
            self.clear_results(remote_dir)
        subprocess.run(
            self._rsync(local + "/", f"{self._require_host()}:{remote_dir}/",
                        excludes=("*.odb", "__pycache__", "*.pyc")),
            check=True)

    def pull(self, remote_dir: str, local_dir: str | Path,
             patterns: Sequence[str] = RESULT_PATTERNS) -> None:
        """Bring selected artefacts back.  Excludes the ODB by default."""
        local = Path(local_dir).resolve()
        local.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            self._rsync(f"{self._require_host()}:{remote_dir}/", str(local) + "/",
                        includes=list(patterns), excludes=("*",)),
            check=True)

    def sbatch(self, remote_dir: str, script: str = "submit.sh") -> SubmitInfo:
        cp = self.run(f"cd {self._expand(remote_dir)} && sbatch {script}", check=True)
        m = _JOB_ID_RE.search(cp.stdout) or _JOB_ID_RE.search(cp.stderr)
        if not m:
            raise RuntimeError(
                f"could not parse sbatch output:\nSTDOUT: {cp.stdout}\n"
                f"STDERR: {cp.stderr}")
        return SubmitInfo(job_id=m.group(1), remote_dir=remote_dir,
                          stdout=cp.stdout)

    def is_queued(self, job_id: str) -> bool:
        cp = self.run(f"squeue -j {shlex.quote(job_id)} -h -o '%T'", check=False)
        # squeue exits non-zero once the job has fully left the queue.
        return cp.returncode == 0 and bool(cp.stdout.strip())

    def wait(self, job_id: str, poll_seconds: int = 30,
             timeout_seconds: int = 6 * 3600, on_tick=None) -> str:
        """Poll until the job leaves the queue; return its final sacct state."""
        start = time.time()
        while self.is_queued(job_id):
            if time.time() - start > timeout_seconds:
                raise TimeoutError(
                    f"job {job_id} still queued after {timeout_seconds} s")
            if on_tick:
                on_tick(job_id, time.time() - start)
            time.sleep(poll_seconds)
        cp = self.run(f"sacct -j {shlex.quote(job_id)} -X -n -o State 2>/dev/null "
                      f"| head -1", check=False)
        out = cp.stdout.strip()
        return out.split()[0] if out else "UNKNOWN"


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(
        prog="python -m movingload.cluster",
        description="Low-level SSH/SLURM helpers. Configure the target with "
                    "MOVINGLOAD_SSH_HOST / _PORT / _REMOTE_ROOT.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="probe ssh, sbatch and abaqus")
    p = sub.add_parser("push"); p.add_argument("local"); p.add_argument("remote")
    p = sub.add_parser("pull"); p.add_argument("remote"); p.add_argument("local")
    p = sub.add_parser("sbatch"); p.add_argument("remote")
    p = sub.add_parser("wait"); p.add_argument("job_id")
    args = ap.parse_args(argv)

    cl = Cluster()
    if args.cmd == "check":
        info = cl.check()
        print(info["stdout"] or info["stderr"])
        return 0 if info["rc"] == 0 else 1
    if args.cmd == "push":
        cl.push(args.local, args.remote)
    elif args.cmd == "pull":
        cl.pull(args.remote, args.local)
    elif args.cmd == "sbatch":
        print(cl.sbatch(args.remote))
    elif args.cmd == "wait":
        print(cl.wait(args.job_id))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
