# The model

Geometry, symmetry, mesh, boundary conditions, load and output — what the deck
actually contains and why.

---

## 1. From a full axle to one tyre

A standard heavy-vehicle axle carries four tyres arranged as two dual-tyre
pairs, one pair per side. Two successive symmetry reductions take that down to
the single tyre this model meshes:

1. **Vehicle centreline.** The left and right pairs are far enough apart that
   their stress fields barely interact; a symmetry plane between them leaves one
   dual-tyre pair (half the axle load) to analyse.
2. **Between the dual tyres.** That pair is itself symmetric about the gap
   between its two tyres. A second plane through that gap halves it again.

What remains is one tyre with a symmetry condition on the plane facing its twin:

```
   x = 0  ←  symmetry plane, midline of the dual-tyre pair
          │
          │      ●  meshed tyre, centre x = X_TYRE = 0.16 m
          │         R = 0.1065 m,  p = 0.70 MPa,  F ≈ 24.9 kN
          │
          ├──────  travel direction  +z  ──────►
          │
          │      (mirror tyre at x = −0.16 m — never meshed; the
          │       u_x = 0 condition at x = 0 produces it exactly)
```

The default numbers are the Chinese **BZZ-100** standard axle: four prints of
0.1065 m radius at 0.70 MPa, `π R² p ≈ 24.9 kN` each, ≈ 100 kN per axle, with
0.32 m between the centres of a dual pair — hence `X_TYRE = 0.16 m`.

Change them through `Tyre(x_centre=…, radius=…, pressure=…)`; the constructor
refuses a print that would cross the symmetry plane (`x_centre < radius`),
because the mirror would then overlap the original and double-count the load.

---

## 2. Coordinates

| axis | range | meaning |
|---|---|---|
| `x` | 0 → 2.25 m | transverse; `x = 0` is the symmetry plane |
| `y` | 0 → 5.00 m | **depth**, positive downward; `y = 0` is the road surface |
| `z` | 0 → 7.00 m | longitudinal, the direction of travel |

Depth on the second axis is deliberate: it makes the longitudinal strain the
**`E33`** component, which is what the history output requests and what every
downstream CSV column is named after. Swapping `y` and `z` would silently change
which component is extracted.

---

## 3. Mesh

52 080 **C3D8R** elements (8-node linear brick, reduced integration),
24 × 31 × 70 in `x` × `y` × `z`.

**Transverse (`x`), three zones.** Uniform 25 mm from 0 to 0.30 m — covering the
print `x ∈ [0.054, 0.267] m` with margin — then 7 medium elements to 1.0 m, then
5 coarse elements out to 2.25 m.

**Longitudinal (`z`), uniform 100 mm.** The load sweeps the whole length, so no
station is privileged; a uniform grid also makes the sampled strain history
evenly spaced in `z`.

**Depth (`y`), per layer.** Each finite layer is divided uniformly by its
`n_elems`. The subgrade is biased towards the surface with an exponent
(`Mesh.subgrade_bias`, default 1.8), packing elements where the gradients are
and stretching them towards the fixed base.

The three preset structures use different subdivisions (`n_elems`) chosen so all
of them total 31 rows — the grids stay the same size, so differences between
structures are structural rather than mesh artefacts. A test enforces this.

**Node and element numbering** is `i` fastest, then `j`, then `k`:

```
node id = 1 + i + j·(NX+1) + k·(NX+1)·(NY+1)
elem id = 1 + i + j·NX     + k·NX·NY
```

The postprocessor locates the observation elements purely through the ids stored
in `<case>_mesh.json`, so changing this numbering invalidates existing results.

---

## 4. Boundary conditions

| face | condition | why |
|---|---|---|
| `x = 0` | `u_x = 0` | symmetry — generates the mirror tyre |
| `x = LX` | `u_x = 0` | lateral roller, far from the load |
| `y = LY` | `u_x = u_y = u_z = 0` | fixed base at 5 m depth |
| `z = 0`, `z = LZ` | `u_z = 0` | longitudinal rollers at both ends |

The domain is a **truncated** box, not an infinite-element model. It is sized so
the boundaries are remote compared with the tyre print (2.25 m ≈ 21 R laterally,
5 m ≈ 47 R deep), but the strain of interest is a *local* quantity right under
the load, far less sensitive to the far field than a surface deflection would
be. If you repurpose this model to report **deflections**, check the truncation
error first — surface deflection converges towards a truncated boundary very
slowly.

---

## 5. Moving load

Abaqus has no built-in load that translates across the mesh, so the print is
imposed through user subroutines evaluated at every surface integration point,
every increment.

`movingload.dload` generates them:

```
DLOAD      vertical pressure     P0 inside the print, 0 outside     (P3NU)
UTRACLOAD  tangential traction   magnitude mu·P0, direction (0,0,±1) (TRVEC3NU)
```

`UTRACLOAD` is emitted **only** when `mu > 0`, and the deck then also carries a
matching `*DSLOAD, TRVEC3NU`. A purely vertical case has neither.

The centre moves linearly,

```
z_c(t) = z_start + (z_end − z_start) · t / t_final
```

and a point is loaded when `(x − X_TYRE)² + (z − z_c)² ≤ R²`.

`z_start = −0.30 m` and `z_end = +7.30 m` sit **outside** the 0–7 m domain, so
the tyre enters and leaves cleanly instead of appearing abruptly at the edge.
A useful consequence: the tyre is over mid-span (`z = 3.5 m`) at exactly
`t = 0.5`, which is the easiest sanity check on any result.

All parameters live in a Fortran `PARAMETER` block, so changing the pressure or
the speed needs only a recompile — no new `.inp`, provided the mesh still
resolves the print.

---

## 6. Step

`*STATIC, DIRECT` with a fixed increment `dt = 0.01` up to `t_final = 1.0`:
100 increments.

The analysis is **quasi-static**. There is no inertia and no damping; "time" is
purely the parameter that advances the load along `z`. So `t_final` does not set
a vehicle speed — it is just the parameterisation. Comparing two speeds would
require a dynamic step, which this model does not do.

`CaseConfig` refuses a `dt` that does not divide `t_final`, because the load
would then stop short of `z_end`.

---

## 7. Output

**Field**, every 10 increments: `U` at nodes, `S` and `E` at integration
points — for viewing in Abaqus/CAE. Frequency 10 keeps the ODB manageable.

**History**, every increment: `E33` on `OBS_ELEMS` — the *observation column*.

The observation column is the bottom element row of the surface layer, at the
`x`-column whose centroid is nearest the tyre centre, sampled along the entire
`z` range. That is the classic fatigue location: maximum longitudinal tensile
strain at the underside of the asphalt, directly beneath the wheel path.

With 70 elements along `z` and 100 increments, the history is a 70 × 100 field
in (position, time). `extract_strain.py` writes it whole
(`<case>_e33_history.csv`) and also isolates the mid-span element
(`<case>_e33_centre.csv`), which is the pulse most analyses actually use.

---

## 8. What the pulse looks like

```
  E33
   ▲
   │            ⌒
   │           ╱ ╲
   │          ╱   ╲
 0 │───────╌─╯     ╰─╌───────
   │
   └──────────────────────────►  t
       approach  peak  leave
```

Asymmetric *approach – peak – leave*, peaking at `t ≈ 0.5`. The magnitude and
sign at that instant are the headline result.

Adding traction tilts the pulse: driving (`+z`) and braking (`−z`) shift the
peak and change its magnitude relative to the vertical-only case, which is the
point of the `accel` / `brake` presets.
