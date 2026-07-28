# AGENTS.md — notes for AI coding agents

For agents working in this repository. After this page you should be able to
build a case, run it, read the results, and avoid the handful of things that
silently produce wrong answers. Deeper detail lives in
[docs/model.md](docs/model.md), [docs/pipeline.md](docs/pipeline.md) and
[docs/cluster.md](docs/cluster.md).

---

## 0. Thirty-second version

`movingload` drives a **moving wheel load** across a 3-D layered pavement in
Abaqus and extracts the longitudinal strain (`E33`) history at the bottom of the
surface layer. Materials are **isotropic linear elastic** — there is no
tension/compression distinction anywhere in this repository.

```bash
export PYTHONPATH=$PWD:$PYTHONPATH
python -m movingload.submit --structure semi_rigid --motion vertical --build-only
python tests/test_build.py        # 32 offline checks, seconds
```

**Deck generation needs neither Abaqus nor a network.** Always verify with
`--build-only` before spending solver time.

---

## 1. This repository is public

Treat every file as world-readable. Do not introduce:

- personal or site paths (`/home/<someone>/…`, site-specific mounts);
- hostnames, IP addresses, usernames, ports, e-mail addresses;
- references to private internal projects, repositories or documents.

`tests/test_build.py::test_no_private_paths_or_hosts_anywhere` scans every
`.py` / `.md` / `.toml` / `.f` in the repo — **including itself** — and fails on
any of the above. `example.org` / `example.com` are the accepted placeholders.

Anything site-specific belongs in an environment variable with a neutral
default, or with no default at all. `movingload.cluster` has **no built-in SSH
target** and raises `ClusterNotConfigured` rather than guessing; keep it that
way.

---

## 2. Build a case

```python
from movingload.config import make_case
from movingload.submit import prepare

out = prepare("semi_rigid", "accel", out_root="runs")
print(out["config"].summary())
```

or from the shell:

```bash
python -m movingload.submit --structure rigid --motion brake --build-only
```

Structures: `flexible`, `semi_rigid`, `rigid`. Motions: `vertical`, `accel`,
`brake`. Any combination is valid.

For something outside the presets, assemble a `CaseConfig` yourself — the
dataclasses in [movingload/config.py](movingload/config.py) validate on
construction, so a bad case fails at build time, not three hours into a solve:

```python
from movingload.config import CaseConfig, Motion, STRUCTURES, Tyre
cfg = CaseConfig(case_id="slow_heavy",
                 structure=STRUCTURES["semi_rigid"],
                 tyre=Tyre(radius=0.12, pressure=0.9e6),
                 motion=Motion(t_final=1.0, mu=0.2, mu_sign=-1.0))
```

### What lands in the run directory

`<case>.inp`, `<case>_mesh.json`, `dload.f`, `extract_strain.py`, and
`submit.sh` (or `run.sh` with `--local`). It is self-contained on purpose: the
whole folder rsyncs to a cluster and solves there with nothing else installed.

### Offline self-check

```bash
grep -c "^\*MATERIAL, NAME=" runs/<case>/<case>.inp   # == number of layers
grep    "interactive"        runs/<case>/submit.sh    # must be present
grep -c "SUBROUTINE"         runs/<case>/dload.f      # 1 vertical, 2 with traction
```

---

## 3. Run it

```bash
python -m movingload.submit --structure rigid --motion brake --local     # here
python -m movingload.submit --structure rigid --motion brake             # cluster
```

Cluster mode needs `MOVINGLOAD_SSH_HOST` and `MOVINGLOAD_REMOTE_ROOT`; it
pushes, `sbatch`es, polls `squeue`, and pulls the CSVs back. `--no-wait`
submits and returns.

The ODB is **not** pulled back — `extract_strain.py` already reduced it to CSV
on the solver host, and an ODB of this model is hundreds of MB.

**Re-running into the same directory.** rsync has no `--delete`, so `push`
first calls `clear_results()` to remove the previous run's output. Without it,
a job that dies before writing its CSVs lets `pull` bring the *old* ones back,
where they look exactly like fresh results. If you drive `cluster` by hand
rather than through `submit`, clear it yourself or check the remote
timestamps.

---

## 4. Read the results

| file | columns |
|---|---|
| `<case>_e33_history.csv` | `time, element_id, z_centre, E33` — whole observation column |
| `<case>_e33_centre.csv` | `time, E33` — mid-span pulse |

`E33` is in strain (not microstrain) and follows Abaqus' sign convention:
**tension positive**. The plotting helpers multiply by 1e6 for display.

```bash
python -m movingload.plots pulse   runs/<case>
python -m movingload.plots heatmap runs/<case>
python -m movingload.plots compare runs/a runs/b runs/c
```

**Sanity check.** The overhead crossing shows up as a **tension spike** at
`motion.time_over(z_obs)`, where `z_obs` is the observation element's *centroid*
(3.55 m on the default grid, so t = 0.507 — not 0.500). If the spike is
somewhere else, the motion parameters baked into `dload.f` and the mesh have
drifted apart.

**Do not "take the peak tension".** The signal is compression-dominated: a
hog–sag–hog bending wave travels with the tyre, so a fixed point sees
compression, then tension as the tyre passes overhead, then compression again.
On the reference case (semi-rigid, vertical) that is −8.8 µε at t = 0.48,
**+4.8 µε at t = 0.51**, and **−9.2 µε at t = 0.54** — the extreme is the
*second compression trough*, not the tensile spike. See
[docs/model.md](docs/model.md#8-what-the-signal-looks-like).

**The peak magnitude is under-resolved at default settings.** The load advances
76 mm per increment, the print is 213 mm across, and the observation element is
100 mm long, so the passage is three or four samples wide. The shape is
converged; the magnitude is not. Refine `Step.dt` *and* `Mesh.dz` together if
you need the number rather than the shape.

---

## 5. Pitfalls

Roughly in order of how easy they are to hit.

**The mesh must resolve the tyre print.** `Mesh.x_fine_to` (default 0.30 m)
has to reach `tyre.x_centre + tyre.radius`. `CaseConfig` checks this and raises,
but if you widen the tyre or move it outward, widen the fine band too.

**`interactive` on the Abaqus command is mandatory.** Without it the command
returns as soon as the solver is spawned, the batch script ends, and the
scheduler tears the allocation down mid-solve. `movingload.slurm` writes it;
don't remove it.

**`extract_strain.py` runs under Abaqus' own Python**, which is 2.7 in most
releases. No f-strings, no `pathlib`, no type annotations, no `ast`-visible
Python-3-only syntax. A test enforces this by parsing the file. Everything
*else* in the package targets Python 3.10+ and may use whatever it likes.

**Element and node numbering is load-bearing.** The extractor finds the
observation elements through the ids recorded in `<case>_mesh.json`. Changing
`mesh.node_id` / `mesh.elem_id`, the layer subdivision, or the grid invalidates
every existing result file. If you must, regenerate the metadata and re-extract.

**`y` is depth, not height.** Positive `y` points *down*. This is why the
longitudinal strain is `E33` and not `E22`; several output requests depend on it.

**The step must divide the motion.** `step.dt` has to fit a whole number of
times into `motion.t_final`, or the load never reaches `z_end`. `CaseConfig`
raises on this.

**Traction cases need both halves.** `mu > 0` must produce *both* a
`UTRACLOAD` subroutine and a `*DSLOAD, TRVEC3NU` in the deck. `movingload`
keeps them in step; if you hand-edit one, fix the other or the traction is
silently dropped.

**Structures are sized to be comparable.** The three presets deliberately share
the same domain and the same 52 080-element count, so differences between them
are structural, not mesh artefacts. A test enforces it. If you add a structure,
keep `sum(n_elems)` equal to the others or say clearly that it is not
comparable.

---

## 6. Working in this repository

- Run `python tests/test_build.py` after any change — it is seconds, offline,
  and covers config validation, mesh, deck, Fortran generation, scripts, run
  directory layout, and the sanitisation scan.
- New generated output? Add an assertion. The tests are the only thing standing
  between a subtly wrong deck and a wasted cluster run.
- Keep deck generation dependency-light: `config`, `inp_builder`, `dload` and
  `slurm` use only the standard library plus numpy (`mesh` needs it for the
  grid). `matplotlib` is confined to `plots.py`.
- Comments explain *why*, not *what*. The existing files tend to note the
  engineering reason behind a choice (why `interactive`, why `y` is depth, why
  the load starts outside the domain); match that.
