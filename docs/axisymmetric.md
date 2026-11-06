# The axisymmetric model

A second route to the same quantity, the strain history at the bottom of the
surface layer under a moving dual tyre. It replaces the 3-D mesh with an
axisymmetric one and supports **linear viscoelastic asphalt** in addition to
linear elastic asphalt.

```bash
python -m movingload.axi --structure semi_rigid --material elastic      --build-only
python -m movingload.axi --structure semi_rigid --material viscoelastic --build-only
python -m movingload.axi --structure semi_rigid --material viscoelastic \
                         --speed 7.6 --speed 16.7          # solve + superpose
python -m movingload.superpose runs/axi_semi_rigid_viscoelastic --speed 20
python -m movingload.plots moving runs/axi_semi_rigid_elastic \
                                  runs/axi_semi_rigid_viscoelastic
```

---

## 1. An axisymmetric model cannot move a load by itself

Only a load centred on the axis is axisymmetric. So the model does not move
anything. It computes the **step response** of one stationary print:

> R(ρ, s): strain at radius ρ from the print centre, at time s after the
> print was switched on and then held.

The moving dual tyre is then **assembled from R by superposition**, in
`movingload.superpose`. The assumptions the 3-D model already makes are what
allow this:

| assumption | what it buys |
|---|---|
| linear material (elastic or linear viscoelastic) | responses add |
| quasi-static (no inertia) | the same as the 3-D `*STATIC` step |
| layers homogeneous in plan | moving a load moves its response |

Under these assumptions the superposition is **exact up to the load-step
discretisation**. It is not the haversine-pulse approximation, where a
stationary pulse at r = 0 cannot tell longitudinal from transverse strain,
because ε_rr = ε_θθ on the axis.

**Not covered:** horizontal traction (`accel` / `brake`). A longitudinal
shear is not axisymmetric and would need Fourier (CAXA) elements. The
axisymmetric route is for vertical loading only.

---

## 2. Time: Boltzmann superposition

The tyre centre advances in steps of `ds` (default 5 mm), so `dt = ds / v`.
During [t_k, t_k+1) the tyre sits at z_k, the same discretisation as the 3-D
`DLOAD`. That is a sum of prints switched on at t_k and off at t_k+1, so at
the observation point P

```
eps(t_n) = Σ_{k≤n} R(P − z_k, t_n − t_k) − Σ_{k<n} R(P − z_k, t_n − t_{k+1})
```

With G[m, k] = R(P − z_k, m·dt) this is a difference of anti-diagonal sums,
ε_n = D_n − D_{n−1}, which numpy evaluates in well under a second.

* **Elastic:** R does not depend on s and the sum telescopes to
  ε_n = R(P − z_n). The stationary field is read as a moving pulse.
* **Viscoelastic:** each point carries the memory of the whole approach. The
  compression built up while the tyre approaches is still present when it
  arrives overhead, and the response lags behind the load.

**The speed enters only here.** One Abaqus run serves every speed and every
tyre layout. The slowest speed that can be synthesised is
`2 · half_length / CreepSteps.t_max` (7 m / 10 s = 0.7 m/s by default). The
code refuses slower speeds instead of extrapolating.

**Travel length.** The print appears at z = −L at full load. For an elastic
layer that is irrelevant. For a viscoelastic one, the far-field creep it
starts is still present when the tyre arrives. With L = 3.5 m (the default,
matching the 3-D domain) the peaks are within 1.5 % of L = 5.5 m; for
example, the overhead strain at 2 m/s is 2.98 µε against 3.02 µε. Pass
`--half-length` up to ≈ 5.5 m (limited by `R_out`) if that matters.

Convergence in `ds` is first order. On the reference viscoelastic case at
7.6 m/s the overhead strain is 2.11 / 2.18 / 2.22 µε for ds = 10 / 5 / 2.5 mm.

## 3. Space: the dual tyre and tensor rotation

R gives radial (E11), vertical (E22) and hoop (E33) strain. Consider a print
whose centre is offset (dx, dz) from P, with ρ = |(dx, dz)|, c = dz/ρ and
s = dx/ρ. Its contribution in the road axes is

```
longitudinal  eps_zz = c² E11 + s² E33
transverse    eps_xx = s² E11 + c² E33
vertical      eps_yy = E22
```

Each tyre of the BZZ-100 dual pair (x = ±0.16 m) is rotated this way and the
contributions are added. P sits under one tyre by default (`--x-obs` moves
it, and `--single` drops the twin). This reproduces the 3-D layout, where
the symmetry plane provides the twin.

The output CSV uses the **3-D model's names**: `E11` transverse, `E22`
vertical, `E33` longitudinal. Its columns are
`time, z_rel, E11, E22, E33`, where `z_rel` is the tyre centre minus P.

---

## 4. The deck

| item | value | why |
|---|---|---|
| element | `CAX4` (or `CAX4R`, `CAX4I`) | fine mesh, so locking is negligible; no hourglass mode at the print edge |
| r grid | 7.5 mm to 1 m, graded to `R_out` = 6 m | a node lies exactly on r = tyre radius; R_out exceeds the travel half-length |
| y grid | ≈ 12.5 mm in finite layers, 30 rows biased in the subgrade | 82 rows; 13 530 elements |
| axial axis | **y = −depth (positive up)** | Abaqus convention; elements counter-clockwise, pressure on face `S3` |
| BCs | u_r = 0 on the axis and at R_out; base fixed | |
| load | `*DSLOAD TYRE_PRINT, P, p0` | stationary, so no user subroutine and no `dload.f` |
| observation | bottom element row of layer 1, centroid (depth 0.144 m on semi-rigid) | the 3-D observation row, finer |

Steps:

```
INSTANT    *STATIC                 load on; instantaneous response (s = 0)
CREEP_00   *VISCO 0 → 1e-5 s       5 fixed increments      (viscoelastic only)
CREEP_01…  *VISCO one per decade   18 fixed increments each, up to 10 s
```

`*VISCOELASTIC` is inactive in `*STATIC`, which uses the instantaneous
moduli. That is exactly why the INSTANT step gives the s = 0+ response. The
`*VISCO` steps use fixed increments (no `CETOL`), so every output time is
known in advance.

## 5. Viscoelastic asphalt

```
E(t) = E0 · [1 − Σ g_i (1 − exp(−t/τ_i))]       constant ν  ⇒  g_i = k_i
```

* The deck writes `*ELASTIC, MODULI=INSTANTANEOUS` with E0. With Prony data,
  Abaqus' default is **LONG TERM**, which would read E0 as E∞ and silently
  stiffen the layer by 1/(1 − Σg).
* Only the asphalt layers are viscoelastic (`ASPHALT_LAYERS`: both AC layers
  of `semi_rigid`, the surface layer of the others). Base, subbase and
  subgrade stay elastic.
* **The shipped Prony series is a placeholder.** It is a smooth shape with
  one term per decade from 0.1 ms to 10 s and E0/E∞ = 50, scaled so that
  E(0.05 s) equals the elastic preset modulus. It is not fitted to any
  mixture or temperature. Replace it with your own fit, for example:

```python
from movingload.axi_config import AxiCase, Viscoelastic
from movingload.config import STRUCTURES
ac = Viscoelastic(E0=28e9, nu=0.30, g=(...), tau=(...))
case = AxiCase(case_id="my_ac", structure=STRUCTURES["semi_rigid"],
               material="viscoelastic", visco={"ASPHALT": ac, "BINDER": ac})
```

Temperature (`*TRS`, WLF shift) is out of scope: the series is assumed to be
at the analysis temperature already.

## 6. Reference results (semi-rigid, dual tyre, placeholder Prony)

Longitudinal strain E33 at the bottom of the surface layer, in µε:

| case | overhead | leading trough | trailing trough |
|---|---|---|---|
| elastic | +7.3 | −8.4 @ −0.21 m | −8.4 @ +0.21 m |
| visco, 2 m/s | +3.0 | −8.9 | −7.8 |
| visco, 7.6 m/s | +2.2 | −6.6 | −5.8 |
| visco, 16.7 m/s | +1.9 | −5.6 | −5.0 |
| visco, 30 m/s | +1.8 | −5.1 | −4.7 |

* **Elastic vs 3-D.** The 3-D model gives −8.8 / +4.8 / −9.2 µε. The troughs
  agree; the 3-D overhead spike is under-resolved (3–4 samples per passage,
  see `docs/model.md` §8). Some difference also remains from the truncated
  box and the 3-D element size.
* **Viscoelastic.** Strain falls as speed rises (stiffer asphalt). The curve
  is asymmetric, with the trailing trough smaller and the tension peak just
  past overhead. Residual strain remains after the passage at low speed.
  These are the qualitative signatures of viscoelastic asphalt. The numbers
  depend entirely on the placeholder Prony series.

## 7. Checks

* **Axis:** E11 = E33 at the first radial centroid (8.77 / 8.76 µε on the
  elastic reference). If they differ, the deck is not axisymmetric about r = 0.
* **v → ∞** reproduces the s = 0 map, and **v → 0** the long-term one. Both
  are tested offline on synthetic kernels.
* **Truncation:** at the largest ρ the tyre reaches (≈ 3.5 m), elastic R is
  +0.37 µε radial, about 5 % of the peaks. That is the far-field level the
  history starts and ends at. Widen `half_length` and `AxiMesh.R_out`
  together if the far field matters.
