"""
movingload — a moving wheel load on a layered elastic pavement, in Abaqus.

A standard dual-tyre axle travels along a three-dimensional half-pavement
model while the longitudinal strain at the bottom of the surface layer is
recorded as a function of time.  That strain–time pulse is the input to
fatigue analysis, and its shape and peak are what this package produces.

The whole path is scripted: no Abaqus/CAE, no hand-edited decks.

    configuration  ->  .inp deck + DLOAD subroutine + run script
                   ->  solve (locally or on a SLURM cluster)
                   ->  ODB reduced to CSV, then plotted

Modules
-------
    config        Structure / Tyre / Motion / Mesh / Step dataclasses and presets
    mesh          structured C3D8R grid and the observation-column bookkeeping
    inp_builder   writes the .inp deck and its metadata sidecar
    dload         generates the DLOAD (and optional UTRACLOAD) Fortran source
    slurm         writes the SLURM batch script or a plain local shell script
    cluster       optional ssh/rsync/sbatch transport (no hard-coded target)
    submit        top-level driver tying the above together
    extract_strain  runs under `abaqus python`: ODB -> CSV
    plots         strain pulse, (t, z) heat map, multi-case comparison

Quick start
-----------
    python -m movingload.submit --structure semi_rigid --motion vertical --build-only

The material model is isotropic linear elasticity throughout — one Young's
modulus per layer, no tension/compression distinction.

Units are SI everywhere: metres, pascals, seconds.  Strains are reported in
Abaqus' native convention, tension positive.
"""

__version__ = "1.0.0"
__all__ = ["config", "mesh", "inp_builder", "dload", "slurm", "cluster",
           "submit", "plots"]
