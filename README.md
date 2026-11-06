# movingload — a moving wheel load on a layered pavement, in Abaqus

A standard dual-tyre axle travels along a three-dimensional half-pavement model
while the **longitudinal strain at the bottom of the surface layer** is recorded
against time. That strain–time pulse is the input to asphalt fatigue analysis,
and its shape and peak are what this package produces.

The whole path is scripted — no Abaqus/CAE, no hand-edited decks:

```
 configuration          .inp deck  +  DLOAD subroutine  +  run script
 (dataclasses)   ──►    ───────────────────────────────────────────
                                        │
                                        ▼
                          solve  (local Abaqus or SLURM cluster)
                                        │
                                        ▼
                        ODB ──► CSV  (abaqus python, on the solver host)
                                        │
                                        ▼
                              strain pulse / heat map / comparison
```

```bash
# build a deck and look at it — no Abaqus needed
python -m movingload.submit --structure semi_rigid --motion vertical --build-only

# solve it here
python -m movingload.submit --structure rigid --motion brake --local

# or on a cluster
export MOVINGLOAD_SSH_HOST=user@cluster.example.org
export MOVINGLOAD_REMOTE_ROOT=/scratch/$USER/movingload
python -m movingload.submit --structure semi_rigid --motion accel
```

Materials are **isotropic and linear** — one Young's modulus per layer; the
axisymmetric model can also make the asphalt **linear viscoelastic**:

```bash
# axisymmetric step response + superposition of the moving dual tyre
python -m movingload.axi --structure semi_rigid --material viscoelastic \
                         --speed 7.6 --speed 16.7
python -m movingload.plots moving runs/axi_semi_rigid_viscoelastic
```

See [docs/axisymmetric.md](docs/axisymmetric.md).
Units are SI throughout (m, Pa, s); strains follow Abaqus' convention,
tension positive.

---

## The model

A tyre print sweeps the surface at constant speed. Abaqus has no built-in
travelling load, so the print is imposed through a `DLOAD` user subroutine
evaluated at every surface integration point and every increment.

**Half model.** The full standard axle has four tyres in two dual-tyre pairs.
Two symmetry reductions leave a single meshed tyre:

1. the vehicle centreline decouples the left and right pairs;
2. the gap between the two tyres of one pair halves it again.

```
   x = 0  ←  symmetry plane (midline between the dual tyres)
          │
          │      ●  meshed tyre, centre x = 0.16 m
          │         R = 0.1065 m,  p = 0.7 MPa,  ≈ 25 kN
          │
          ├───────  travel direction  +z  ───────►
          │
          │      (mirror tyre at x = −0.16 m, produced by u_x = 0 at x = 0)
```

**Coordinates.** `x` transverse, `y` depth (positive *down*), `z` longitudinal.
Depth on `y` is what makes the longitudinal strain the `E33` component.

**Domain and mesh.** 2.25 × 5.00 × 7.00 m, 52 080 C3D8R elements. The `x` grid
is fine (25 mm) across the tyre print and coarsens outward; `z` is uniform at
100 mm; `y` is subdivided per layer and biased towards the surface in the
subgrade.

**Step.** `*STATIC, DIRECT`, Δt = 0.01 over t = 1.0. The analysis is
quasi-static: there is no inertia, and time is purely the parameter that
advances the load. The tyre starts at z = −0.30 m and ends at z = +7.30 m, so it
enters and leaves the domain cleanly and crosses mid-span at exactly t = 0.5.

Full detail — boundary conditions, mesh grading, output requests — in
[docs/model.md](docs/model.md).

---

## Built-in cases

Three pavement structures, all with the same total depth and the same element
count so results stay comparable:

| `--structure` | surface | base | subbase | subgrade |
|---|---|---|---|---|
| `flexible` | AC 0.10 m, 5 GPa | granular 0.20 m, 0.40 GPa | granular 0.20 m, 0.15 GPa | 0.060 GPa |
| `semi_rigid` | AC 0.15 m, 10 GPa | AC binder 0.30 m, 7 GPa | granular 0.20 m, 0.30 GPa | 0.060 GPa |
| `rigid` | AC 0.05 m, 10 GPa | PCC 0.25 m, 30 GPa | CTB 0.20 m, 5 GPa | 0.10 GPa |

Three load motions:

| `--motion` | surface load |
|---|---|
| `vertical` | pressure only |
| `accel` | pressure + tangential traction μ = 0.30 in **+z** (driving) |
| `brake` | pressure + tangential traction μ = 0.30 in **−z** (braking) |

Traction cases additionally emit a `UTRACLOAD` subroutine and a `*DSLOAD` on the
same faces; vertical cases emit neither.

Any combination works, e.g. `--structure rigid --motion brake`. To go beyond the
presets, build a `CaseConfig` directly — see [docs/pipeline.md](docs/pipeline.md).

---

## Installing

```bash
git clone <this repository>
cd pavement-moving-load
pip install -e .            # or just: export PYTHONPATH=$PWD:$PYTHONPATH
```

Requires **Python 3.10+** and **numpy**. `matplotlib` is needed only for the
figures (`pip install -e ".[plots]"`). Cluster submission needs `ssh` and
`rsync`. Solving needs **Abaqus/Standard** with a Fortran compiler configured for
user subroutines (the axisymmetric model needs no subroutine and no compiler).

Deck generation touches neither Abaqus nor the network, so the whole modeling
stage — and the entire test suite — runs anywhere:

```bash
python tests/test_build.py        # 50 checks, a few seconds
```

---

## A run directory

Every case gets a self-contained folder under `runs/`:

```
runs/semi_rigid_accel/
├── semi_rigid_accel.inp          the deck
├── semi_rigid_accel_mesh.json    mesh + observation metadata
├── dload.f                       DLOAD (+ UTRACLOAD for traction cases)
├── extract_strain.py             ODB → CSV, runs under `abaqus python`
└── submit.sh                     SLURM script  (or run.sh with --local)
```

Nothing outside the folder is needed to solve it, which is why it can be rsynced
to a cluster as-is. After the solve, two CSVs appear alongside:

| file | columns | contents |
|---|---|---|
| `<case>_e33_history.csv` | `time, element_id, z_centre, E33` | the whole observation column — the strain wave travelling with the tyre |
| `<case>_e33_centre.csv` | `time, E33` | the mid-span pulse |

---

## Figures

```bash
python -m movingload.plots pulse   runs/semi_rigid_vertical
python -m movingload.plots heatmap runs/semi_rigid_vertical
python -m movingload.plots compare runs/flexible_vertical \
                                   runs/semi_rigid_vertical \
                                   runs/rigid_vertical
```

The mid-span signature is **not** a single bump. It is the hog–sag–hog bending
wave of the pavement travelling past a fixed point:

```
  E33 (µε)
   ▲
 +5│                    ╱╲            ← sag: tyre directly overhead
   │                   ╱  ╲
   │  ┄┄┄┄┄┄┄┄┄┄┄┄╮   ╱    ╲   ╭┄┄┄┄┄┄┄┄┄┄┄┄   ← distant: mild tension
  0│──────────────╲──╱──────╲─╱──────────────
   │               ╲╱        ╲╱
   │               ▼          ▼   ← hog: tyre ~0.2–0.3 m before / after
-10│
   └────────────────────────────────────────►  t
                 0.48  0.51  0.54
```

Measured on the reference case (semi-rigid, vertical): a +1.7 µε plateau while
the tyre is elsewhere, a compression trough of −8.8 µε as it approaches, a
+4.8 µε tension spike as it passes directly overhead, and a second trough of
−9.2 µε as it leaves. **The extreme value is compressive**, and it is the second
trough, not the overhead instant.

The tension spike marks the overhead crossing: it lands at
`motion.time_over(z_obs)` where `z_obs` is the observation element's centroid
(3.55 m by default, so t = 0.507).

> **Resolution caveat.** With the default `dt = 0.01` the load advances 76 mm per
> increment, while the tyre print is 213 mm across and the observation element is
> 100 mm long — the entire passage is captured in three or four samples. The
> *shape* is converged; the *peak magnitude* is not. Refine `Step.dt`, and
> `Mesh.dz` with it, if the magnitude is what you need.

`compare` also prints a peak table (maximum, minimum, and the times they occur).

---

## Cluster configuration

There is **no built-in host, port or remote path** — the transport is configured
entirely by you:

| variable | meaning | default |
|---|---|---|
| `MOVINGLOAD_SSH_HOST` | `user@host` | *(none — required)* |
| `MOVINGLOAD_SSH_PORT` | SSH port | `22` |
| `MOVINGLOAD_REMOTE_ROOT` | where run directories land | *(none — required)* |
| `MOVINGLOAD_ABAQUS_MODULE` | e.g. `module load abaqus/2023` | *(none)* |
| `MOVINGLOAD_ABAQUS_CMD` | Abaqus executable | `abaqus` |
| `MOVINGLOAD_PARTITION` | scheduler partition | *(site default)* |

An unconfigured `Cluster` refuses to do anything rather than guessing.
`--local` never touches this module at all. See
[docs/cluster.md](docs/cluster.md).

---

## Repository layout

```
movingload/
├── config.py          Structure / Tyre / Motion / Mesh / Step + presets
├── mesh.py            structured C3D8R grid, observation bookkeeping
├── inp_builder.py     writes the .inp and its metadata sidecar
├── dload.py           generates the DLOAD / UTRACLOAD Fortran
├── slurm.py           SLURM batch script or plain local shell script
├── cluster.py         optional ssh/rsync/sbatch transport
├── submit.py          top-level driver
├── extract_strain.py  ODB → CSV, under `abaqus python`
├── axi_config.py      axisymmetric: Viscoelastic (Prony) / AxiMesh / CreepSteps / AxiCase
├── axi_mesh.py        structured CAX4 grid, print edge node, observation row
├── axi_inp.py         step-response deck (*STATIC + *VISCO) and metadata
├── extract_axi.py     ODB → R(r, s) CSV, under `abaqus python`
├── superpose.py       moving dual tyre from R(r, s): Boltzmann + tensor rotation
├── axi.py             axisymmetric driver
└── plots.py           pulse, heat map, comparison, moving (axisymmetric)
docs/                  model.md, pipeline.md, cluster.md, axisymmetric.md
tests/test_build.py    50 offline checks
```

---

## License

MIT — see [LICENSE](LICENSE).
