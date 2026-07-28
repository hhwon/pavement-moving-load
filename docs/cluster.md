# Running on a cluster

Optional. If Abaqus is installed on your own machine, use
`python -m movingload.submit --local` and skip this page entirely —
`movingload.cluster` is never imported in that mode.

---

## Configuration

**Nothing is hard-coded.** There is no built-in host, port, user or remote path;
an unconfigured `Cluster` raises `ClusterNotConfigured` rather than guessing.

| variable | meaning | default |
|---|---|---|
| `MOVINGLOAD_SSH_HOST` | `user@host` | *(none — required)* |
| `MOVINGLOAD_SSH_PORT` | SSH port | `22` |
| `MOVINGLOAD_REMOTE_ROOT` | parent directory for run directories | *(none — required)* |
| `MOVINGLOAD_ABAQUS_MODULE` | command run before Abaqus, e.g. `module load abaqus/2023` | *(none)* |
| `MOVINGLOAD_ABAQUS_CMD` | Abaqus executable name | `abaqus` |
| `MOVINGLOAD_PARTITION` | scheduler partition | *(site default)* |

```bash
export MOVINGLOAD_SSH_HOST=user@cluster.example.org
export MOVINGLOAD_REMOTE_ROOT=/scratch/$USER/movingload
export MOVINGLOAD_ABAQUS_MODULE="module load abaqus/2023"
```

Or in code:

```python
from movingload.cluster import Cluster
cl = Cluster(host="user@cluster.example.org", port="22",
             remote_root="/scratch/me/movingload")
```

Locally you need only `ssh` and `rsync`. Key-based authentication is assumed —
`BatchMode=yes` is set, so a password prompt fails rather than hangs.

---

## Requirements on the cluster

- **Abaqus/Standard**, plus **a Fortran compiler on the path**. Abaqus ships no
  compiler of its own, and this model will not run without one: the moving load
  *is* a user subroutine. At most sites the compiler is a separate module, so
  load both:

  ```bash
  --module-load "module load compiler/intel/2021.3.0 apps/abaqus/2020"
  ```

  Confirm once with `abaqus verify -user_std` on a compute node.
- **SLURM** (`sbatch`, `squeue`, `sacct`). Some sites reject jobs that do not
  name a partition — see the troubleshooting table.
- Nothing else. The run directory ships its own `extract_strain.py`, and that
  runs under Abaqus' bundled Python — no site Python packages required.

Probe the target:

```bash
python -m movingload.cluster check
```

reports the remote hostname and whether `sbatch` and `abaqus` are on the path.
`abaqus` is often unavailable in the login shell but present after a module
load on compute nodes; that is fine, and it is why `MOVINGLOAD_ABAQUS_MODULE`
exists.

---

## One command

```bash
python -m movingload.submit --structure semi_rigid --motion accel \
    --cpus 8 --mem 32g --walltime 02:00:00
```

1. build the run directory locally;
2. `rsync` it to `$MOVINGLOAD_REMOTE_ROOT/<case_id>`;
3. `sbatch submit.sh`;
4. poll `squeue` until the job leaves the queue;
5. `rsync` the CSVs, JSON and logs back.

`--no-wait` stops after step 3; collect later with

```bash
python -m movingload.cluster pull /scratch/$USER/movingload/<case_id> runs/<case_id>
```

**The ODB is never transferred.** `extract_strain.py` already reduced it to a
few kilobytes of CSV on the solver host, and the ODB for this model runs to
hundreds of megabytes. Fetch it by hand if you really need it in CAE:

```bash
rsync -av --progress -e "ssh -p $MOVINGLOAD_SSH_PORT" \
      "$MOVINGLOAD_SSH_HOST:$MOVINGLOAD_REMOTE_ROOT/<case_id>/<case_id>.odb" .
```

---

## The batch script

```bash
#!/bin/bash
#SBATCH --job-name=<case_id>
#SBATCH --output=<case_id>_%j.log
#SBATCH --error=<case_id>_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=<if configured>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

<module load line, if configured>

abaqus job=<case_id> input=<case_id>.inp user=dload.f cpus=8 interactive
abaqus python extract_strain.py <case_id>.odb
```

Two details worth keeping:

**`interactive` is mandatory.** Without it the `abaqus` command returns as soon
as the solver is spawned, the batch script reaches its end, and SLURM tears the
allocation down mid-solve. The failure looks like a job that "completed" in
seconds with no ODB.

**The ODB is reduced in place.** Running `extract_strain.py` as the last step of
the job means only CSVs need to cross the network.

---

## Step-by-step

```bash
python -m movingload.cluster check
python -m movingload.cluster push   runs/<case>  /scratch/$USER/movingload/<case>
python -m movingload.cluster sbatch /scratch/$USER/movingload/<case>
python -m movingload.cluster wait   <job_id>
python -m movingload.cluster pull   /scratch/$USER/movingload/<case>  runs/<case>
```

Python API:

```python
from movingload.cluster import Cluster
cl = Cluster()
remote = cl.remote_dir_for("semi_rigid_accel")
cl.push("runs/semi_rigid_accel", remote)
job = cl.sbatch(remote)                      # -> SubmitInfo(job_id, ...)
state = cl.wait(job.job_id, poll_seconds=30) # -> 'COMPLETED' / 'FAILED' / ...
cl.pull(remote, "runs/semi_rigid_accel")
cl.run("squeue -u $USER")                    # any remote command
```

`push` excludes `*.odb`, `__pycache__` and `*.pyc`. `pull` takes
`RESULT_PATTERNS` (`*.csv *.json *.log *.err *.sta *.msg *.dat`) — override
`patterns=` to widen or narrow it.

### Re-using a run directory

rsync runs **without `--delete`**, so pushing into a directory that already
holds a previous run leaves that run's output in place. If the new job then
dies *before* writing its own CSVs, `pull` brings the **old** ones home, and
stale numbers are indistinguishable from fresh ones.

`push` therefore calls `clear_results()` first: it deletes the output globs
(`*.csv *.odb *.dat *.msg *.sta *.fil *.log …`, plus `fort.*`) and leaves the
inputs — deck, `dload.f`, scripts, mesh sidecar — alone. A failed run then
shows up as a missing file instead of a wrong answer.

Because this is an `rm -f` on the cluster it refuses any path that is shallow
or contains shell metacharacters, returning `False` without acting. To keep
existing results, push with `clear_stale=False`.

When in doubt, check the timestamps yourself:

```bash
ssh -p "$MOVINGLOAD_SSH_PORT" "$MOVINGLOAD_SSH_HOST" \
    "ls -la --time-style=+%m-%d_%H:%M $MOVINGLOAD_REMOTE_ROOT/<case_id>"
```

---

## Troubleshooting

Start with `<case_id>_<jobid>.log` (script stdout) and `.err`, then Abaqus'
own `.msg` (increment-by-increment history) and `.sta` (convergence summary).

| symptom | likely cause |
|---|---|
| `ClusterNotConfigured` | `MOVINGLOAD_SSH_HOST` or `MOVINGLOAD_REMOTE_ROOT` unset |
| ssh fails immediately | key-based auth not set up (`BatchMode=yes` will not prompt) |
| `Invalid account or account/partition combination specified` | the site has no usable default partition for your account — pass `--partition <name>` or set `MOVINGLOAD_PARTITION`. `sinfo -s` lists what is available |
| job "completes" in seconds, no ODB | `interactive` missing from the abaqus line |
| `abaqus: command not found` in the log | set `MOVINGLOAD_ABAQUS_MODULE`, or `MOVINGLOAD_ABAQUS_CMD` if the executable has a versioned name |
| link/compile error mentioning `dload.f` | no Fortran compiler on the path. Abaqus ships no compiler of its own — load one **in addition to** Abaqus, e.g. `--module-load "module load compiler/intel/2021.3.0 apps/abaqus/2020"`. Confirm with `abaqus verify -user_std` |
| `Requested node configuration is not available` | `--mem` / `--cpus` exceed what the partition offers |
| killed at the wall clock, no `done` line | raise `--walltime` |
| no CSVs came back | the solve failed before extraction — read the log |
| `no mesh metadata found` | `<case>_mesh.json` was not shipped next to the ODB |
| warning about elements with no `E33` | the deck lost its history request on `OBS_ELEMS` |

Useful remote commands:

```bash
python -m movingload.cluster --help
ssh -p "$MOVINGLOAD_SSH_PORT" "$MOVINGLOAD_SSH_HOST" 'squeue -u $USER'
ssh -p "$MOVINGLOAD_SSH_PORT" "$MOVINGLOAD_SSH_HOST" \
    'sacct -j <jobid> -X -o JobID,State,Elapsed,MaxRSS,ReqMem'
```

`MaxRSS` from `sacct` is the honest basis for the next `--mem` request: take the
observed peak and add roughly 30 %.
