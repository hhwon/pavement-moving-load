# Pipeline

Configuration → deck → solve → CSV → figures, with the API for each stage.
For the physics see [model.md](model.md); for the transport see
[cluster.md](cluster.md).

This page describes the 3-D route. The axisymmetric route has a different
pipeline: a stationary step-response deck with no `dload.f`, `extract_axi.py`
producing R(r, s), then `movingload.superpose` building the moving dual tyre.
It is described in [axisymmetric.md](axisymmetric.md).

---

## Stage 0 — configuration

Five frozen dataclasses in [`movingload/config.py`](../movingload/config.py).
All of them validate on construction, so a malformed case fails immediately
rather than three hours into a solve.

### `Structure`

An isotropic layered pavement, top layer first. `h` lists the **finite** layer
thicknesses; the subgrade fills whatever depth remains up to `Mesh.LY`, so
`E`, `nu`, `names` and `n_elems` each have one entry more than `h`.

```python
Structure(name="semi_rigid", description="AC + binder over a cement-treated base",
          h=[0.15, 0.30, 0.20],                       # subgrade excluded
          E=[10.0e9, 7.0e9, 0.30e9, 0.060e9],         # subgrade last
          nu=[0.25, 0.25, 0.35, 0.40],
          names=["ASPHALT", "BINDER", "GRANULAR", "SUBGRADE"],
          n_elems=[6, 6, 4, 15])
```

Presets: `STRUCTURES["flexible" | "semi_rigid" | "rigid"]`.

### `Tyre`

```python
Tyre(x_centre=0.16, radius=0.1065, pressure=0.70e6)   # metres, pascals
```

`.force` = `π R² p`; `.dual_spacing` = `2 · x_centre`. Rejects
`x_centre < radius` — that print would cross the symmetry plane.

### `Motion`

```python
Motion(z_start=-0.30, z_end=7.30, t_final=1.0, mu=0.0, mu_sign=1.0)
```

`mu` is a magnitude and `mu_sign` the direction (`+1` driving, `−1` braking);
`mu = 0` means vertical load only. `.horizontal` reports whether a `UTRACLOAD`
is needed. `.z_at(t)` and `.time_over(z)` are inverses — use the latter to
predict where the pulse should peak.

Presets: `MOTIONS["vertical" | "accel" | "brake"]`.

### `Mesh` and `Step`

`Mesh` holds the domain (`LX/LY/LZ`) and the grading knobs (`x_fine_to`,
`dx_fine`, `x_med_to`, `n_med`, `n_coarse`, `dz`, `subgrade_bias`).
`Step` holds `dt`, `max_inc` and the two output frequencies.

### `CaseConfig`

Ties them together and enforces two cross-cutting rules:

- the fine `x` band must reach `tyre.x_centre + tyre.radius`;
- `step.dt` must divide `motion.t_final` exactly.

```python
from movingload.config import make_case
cfg = make_case("semi_rigid", "accel")        # case_id "semi_rigid_accel"
print(cfg.summary())
```

`make_case` also takes `case_id=` and passes any other keyword through, so
`make_case("rigid", "brake", mesh=Mesh(dz=0.05))` refines the longitudinal grid.

---

## Stage 1 — deck

### Grid

```python
from movingload.mesh import build_grid, describe
g = build_grid(cfg)
print(describe(g))
```

`Grid` exposes `NX/NY/NZ`, `n_elements`, `n_nodes`, `node_id(i,j,k)`,
`elem_id(i,j,k)`, and the observation bookkeeping: `i_obs`, `j_obs`,
`obs_elem_ids`, `obs_x`, `z_centres`.

### `.inp` and metadata

```python
from movingload.inp_builder import build, deck_lines, metadata
info = build(cfg, "runs/semi_rigid_accel")
```

writes `<case_id>.inp` and `<case_id>_mesh.json`. `deck_lines(cfg, grid)`
returns the deck as a list of strings without touching the filesystem — handy
for tests and for diffing two configurations.

The deck contains, in order: heading, nodes, C3D8R elements, one `*ELSET` per
layer, `TOP_ELEMS` (loaded faces) and `OBS_ELEMS` (observation column), the five
boundary node sets, `*SOLID SECTION` + `*MATERIAL` + `*ELASTIC` per layer, the
`*BOUNDARY` block, and the step with `*DLOAD` (plus `*DSLOAD` for traction
cases) and the field/history output requests.

`<case_id>_mesh.json` is the contract with the postprocessor: node coordinates,
grid sizes, layer properties, tyre and motion parameters, and the observation
element ids with their location. `extract_strain.py` reads nothing else.

### Fortran subroutine

```python
from movingload.dload import build as build_dload, source
build_dload(cfg, "runs/semi_rigid_accel")     # writes dload.f
print(source(cfg.tyre, cfg.motion))           # or just look at it
```

Emits `DLOAD` always, and `UTRACLOAD` only when `cfg.motion.horizontal`. Fixed-form
Fortran, all parameters in a `PARAMETER` block.

### Run script

```python
from movingload.slurm import Resources, slurm_script, local_script, write
write("runs/x/submit.sh", slurm_script("x", "x.inp",
      resources=Resources(cpus=8, mem="32g", walltime="02:00:00")))
```

Both flavours share one body: load the module (if configured), solve with
`interactive`, then `abaqus python extract_strain.py <job>.odb`. Site specifics
come from `Resources` or the environment (`MOVINGLOAD_ABAQUS_MODULE`,
`MOVINGLOAD_ABAQUS_CMD`, `MOVINGLOAD_PARTITION`), all with neutral defaults.

### All at once

```python
from movingload.submit import prepare
out = prepare("semi_rigid", "accel", out_root="runs")
```

or `python -m movingload.submit --structure semi_rigid --motion accel --build-only`.

Produces a self-contained directory:

```
runs/semi_rigid_accel/
├── semi_rigid_accel.inp
├── semi_rigid_accel_mesh.json
├── dload.f
├── extract_strain.py          copied in, so the folder needs nothing else
└── submit.sh                  (or run.sh with --local)
```

---

## Stage 2 — solve

```bash
python -m movingload.submit --structure rigid --motion brake --local
python -m movingload.submit --structure rigid --motion brake      # cluster
```

| flag | effect |
|---|---|
| `--build-only` | write the directory and stop |
| `--local` | run `run.sh` here with a local Abaqus |
| *(neither)* | push to the configured cluster, `sbatch`, wait, pull |
| `--no-wait` | submit and return; pull later |
| `--cpus/--mem/--walltime/--partition` | scheduler resources |
| `--abaqus-cmd`, `--module-load` | site overrides |

---

## Stage 3 — extract

Run automatically by both scripts at the end of the solve:

```bash
abaqus python extract_strain.py <case>.odb
```

It finds `<case>_mesh.json` beside the ODB, reads `OBS_ELEM_IDS`, walks the
history regions of the last step, and writes:

| file | columns |
|---|---|
| `<case>_e33_history.csv` | `time, element_id, z_centre, E33` |
| `<case>_e33_centre.csv` | `time, E33` |

It warns if any observation element carried no `E33` output (usually a deck that
lost its history request) and exits non-zero if nothing was collected at all.

> This file runs under **Abaqus' own Python**, 2.7 in most releases: no
> f-strings, no `pathlib`, no annotations. A test parses it to enforce this.
> The rest of the package targets Python 3.10+.

---

## Stage 4 — figures

```bash
python -m movingload.plots pulse   runs/semi_rigid_vertical
python -m movingload.plots heatmap runs/semi_rigid_vertical
python -m movingload.plots compare runs/flexible_vertical runs/semi_rigid_vertical \
                                   runs/rigid_vertical --title "structure comparison"
```

| command | figure |
|---|---|
| `pulse` | mid-span strain–time curve, peak annotated |
| `heatmap` | `E33` over (time, `z`) with the tyre trajectory overlaid |
| `compare` | several mid-span pulses overlaid, plus a printed peak table |

Library equivalents: `plot_pulse`, `plot_heatmap`, `plot_compare`, and
`load_case(run_dir)` which returns the metadata and both CSVs as arrays.
Strains are converted to microstrain for display; the CSVs stay in strain.

---

## Checklist

| check | expected |
|---|---|
| `*MATERIAL, NAME=` count | number of layers |
| `interactive` in the run script | present |
| `SUBROUTINE` count in `dload.f` | 1 vertical, 2 with traction |
| `*DSLOAD` in the deck | only for traction cases |
| tension-spike time | `motion.time_over(z_obs)`, `z_obs` = observation element *centroid* — 0.507 by default, **not** 0.5 |
| extreme value | compressive, at the second trough — see [model.md §8](model.md#8-what-the-signal-looks-like) |
| element count across the three presets | identical (52 080) |

Everything above except the two signal checks is covered by
`python tests/test_build.py`, which needs neither Abaqus nor a network.
