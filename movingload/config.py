"""
Model configuration: pavement structure, tyre, motion, mesh, domain.

Everything the deck generator needs, as plain dataclasses.  All values are
SI: metres, pascals, seconds.  Nothing here touches the filesystem or the
network, so a configuration can be built and inspected offline.

Coordinate system (Abaqus global axes)

    x  transverse   0 -> LX ;  x = 0 is the symmetry plane
    y  depth        0 -> LY ;  y = 0 is the road surface, y positive DOWN
    z  longitudinal 0 -> LZ ;  direction of travel

Note the unusual choice of y as depth: it makes the longitudinal strain the
E33 component, which is what the observation output requests.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence


# ----------------------------------------------------------------------
# Pavement structure
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Structure:
    """An isotropic layered pavement, top layer first.

    `h` lists the FINITE layer thicknesses; the subgrade is whatever depth
    remains between their sum and the domain depth `LY`, so `E` / `nu` /
    `names` / `n_elems` each carry one more entry than `h`.
    """
    name: str
    description: str
    h: List[float]              # finite layer thicknesses [m] (subgrade excluded)
    E: List[float]              # Young's modulus per layer, subgrade last [Pa]
    nu: List[float]             # Poisson ratio per layer, subgrade last [-]
    names: List[str]            # Abaqus material / elset names, subgrade last
    n_elems: List[int]          # elements through the thickness, subgrade last

    def __post_init__(self):
        n = len(self.h) + 1
        for label, seq in (("E", self.E), ("nu", self.nu),
                           ("names", self.names), ("n_elems", self.n_elems)):
            if len(seq) != n:
                raise ValueError(
                    f"{self.name}: len({label}) = {len(seq)}, expected {n} "
                    f"(= len(h) + 1, the extra entry being the subgrade)")
        if any(x <= 0 for x in self.h):
            raise ValueError(f"{self.name}: all thicknesses must be > 0")
        if any(x <= 0 for x in self.E):
            raise ValueError(f"{self.name}: all moduli must be > 0")
        if any(not (0.0 <= x < 0.5) for x in self.nu):
            raise ValueError(f"{self.name}: all Poisson ratios must be in [0, 0.5)")
        if any(x < 1 for x in self.n_elems):
            raise ValueError(f"{self.name}: every layer needs >= 1 element")

    def thicknesses(self, LY: float) -> List[float]:
        """Layer thicknesses including the subgrade, which fills out to LY."""
        h_sub = LY - sum(self.h)
        if h_sub <= 0:
            raise ValueError(
                f"{self.name}: finite layers ({sum(self.h)} m) do not fit in "
                f"a domain of depth LY = {LY} m")
        return list(self.h) + [h_sub]

    def summary(self, LY: float = 5.0) -> str:
        out = [f"Structure: {self.name} — {self.description}"]
        for nm, h, E, nu, ne in zip(self.names, self.thicknesses(LY),
                                    self.E, self.nu, self.n_elems):
            out.append(f"  {nm:>10s}: h={h:6.3f} m  E={E/1e9:7.3f} GPa  "
                       f"nu={nu:.2f}  ({ne} elements)")
        return "\n".join(out)


#: Three reference structures.  Each keeps the same total depth, so meshes
#: differ only through the layer subdivision.
STRUCTURES: Dict[str, Structure] = {
    "flexible": Structure(
        name="flexible",
        description="thin AC over unbound granular base and subbase",
        h=[0.10, 0.20, 0.20],
        E=[5.0e9, 0.40e9, 0.15e9, 0.060e9],
        nu=[0.30, 0.35, 0.35, 0.40],
        names=["ASPHALT", "GRBASE", "GRSUBBASE", "SUBGRADE"],
        n_elems=[4, 6, 4, 17],
    ),
    "semi_rigid": Structure(
        name="semi_rigid",
        description="AC surface + AC binder over a cement-treated base",
        h=[0.15, 0.30, 0.20],
        E=[10.0e9, 7.0e9, 0.30e9, 0.060e9],
        nu=[0.25, 0.25, 0.35, 0.40],
        names=["ASPHALT", "BINDER", "GRANULAR", "SUBGRADE"],
        n_elems=[6, 6, 4, 15],
    ),
    "rigid": Structure(
        name="rigid",
        description="thin AC overlay on a portland-cement concrete slab",
        h=[0.05, 0.25, 0.20],
        E=[10.0e9, 30.0e9, 5.0e9, 0.10e9],
        nu=[0.30, 0.18, 0.25, 0.40],
        names=["ASPHALT", "PCC", "CTB", "SUBGRADE"],
        n_elems=[3, 8, 4, 16],
    ),
}


# ----------------------------------------------------------------------
# Tyre and motion
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Tyre:
    """One circular tyre print in the half model.

    Defaults are the Chinese BZZ-100 standard axle: four tyres in two
    dual-tyre pairs, 100 kN total, 0.7 MPa contact pressure.  See
    `docs/model.md` for how the full axle reduces to this single tyre.
    """
    x_centre: float = 0.16       # tyre centre, = half the dual-tyre spacing [m]
    radius: float = 0.1065       # contact radius [m]
    pressure: float = 0.70e6     # contact pressure [Pa]

    def __post_init__(self):
        if self.radius <= 0 or self.pressure <= 0:
            raise ValueError("tyre radius and pressure must be > 0")
        if self.x_centre < self.radius:
            raise ValueError(
                f"tyre print would cross the symmetry plane: x_centre="
                f"{self.x_centre} < radius={self.radius}")

    @property
    def force(self) -> float:
        """Vertical force on one tyre print [N]."""
        return math.pi * self.radius ** 2 * self.pressure

    @property
    def dual_spacing(self) -> float:
        """Centre-to-centre spacing of the dual-tyre pair [m]."""
        return 2.0 * self.x_centre


@dataclass(frozen=True)
class Motion:
    """Straight-line travel of the tyre centre along +z, plus optional
    tangential traction.

    `z_start` / `z_end` are deliberately outside the domain so the tyre
    enters and leaves cleanly:  z_c(t) = z_start + (z_end - z_start)·t/t_final.

    `mu` is the tyre/pavement friction coefficient and `mu_sign` its
    direction: +1 traction in +z (acceleration), -1 in -z (braking),
    and `mu = 0` for a purely vertical load.
    """
    z_start: float = -0.30
    z_end: float = 7.30
    t_final: float = 1.00
    mu: float = 0.0
    mu_sign: float = 1.0

    def __post_init__(self):
        if self.t_final <= 0:
            raise ValueError("t_final must be > 0")
        if self.mu < 0:
            raise ValueError("mu is a magnitude; use mu_sign for the direction")
        if self.mu_sign not in (1.0, -1.0):
            raise ValueError("mu_sign must be +1.0 or -1.0")

    @property
    def horizontal(self) -> bool:
        """True when a tangential traction is present (needs UTRACLOAD)."""
        return self.mu > 0.0

    def z_at(self, t: float) -> float:
        return self.z_start + (self.z_end - self.z_start) * t / self.t_final

    def time_over(self, z: float) -> float:
        """Step time at which the tyre centre passes over longitudinal z."""
        return (z - self.z_start) / (self.z_end - self.z_start) * self.t_final


#: Named motion presets used by the reference study.
MOTIONS: Dict[str, Motion] = {
    "vertical": Motion(mu=0.0),
    "accel":    Motion(mu=0.30, mu_sign=+1.0),
    "brake":    Motion(mu=0.30, mu_sign=-1.0),
}


# ----------------------------------------------------------------------
# Domain, mesh and step
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Mesh:
    """Domain extent and grid.

    The x grid is three zones: a uniform fine band covering the tyre print
    with margin, a medium band, then a coarse band out to LX.  z is uniform
    (the load sweeps through it, so no zone is privileged).  y is subdivided
    per layer, with the subgrade biased towards the surface by `subgrade_bias`
    (an exponent > 1 packs elements near the top).
    """
    LX: float = 2.25
    LY: float = 5.00
    LZ: float = 7.00

    x_fine_to: float = 0.30      # end of the uniform fine band [m]
    dx_fine: float = 0.025
    x_med_to: float = 1.00
    n_med: int = 7
    n_coarse: int = 5

    dz: float = 0.10
    subgrade_bias: float = 1.8

    def __post_init__(self):
        if not (0 < self.x_fine_to < self.x_med_to < self.LX):
            raise ValueError("need 0 < x_fine_to < x_med_to < LX")
        if self.dx_fine <= 0 or self.dz <= 0:
            raise ValueError("element sizes must be > 0")
        if self.subgrade_bias < 1.0:
            raise ValueError("subgrade_bias < 1 would coarsen towards the surface")

    def covers_tyre(self, tyre: Tyre) -> bool:
        """True when the fine x band spans the whole tyre print."""
        return self.x_fine_to >= tyre.x_centre + tyre.radius


@dataclass(frozen=True)
class Step:
    """Quasi-static step driving the load along z.

    `*STATIC, DIRECT` with a fixed increment: there is no inertia, time is
    purely the parameter that advances the load.  With dt = 0.01 over
    t_final = 1.0 the strain history has 100 samples.
    """
    dt: float = 0.01
    max_inc: int = 200
    field_frequency: int = 10    # write U/S/E every N increments
    history_frequency: int = 1   # write E33 every increment

    def __post_init__(self):
        if self.dt <= 0:
            raise ValueError("dt must be > 0")


@dataclass(frozen=True)
class CaseConfig:
    """One complete analysis case."""
    case_id: str
    structure: Structure
    tyre: Tyre = field(default_factory=Tyre)
    motion: Motion = field(default_factory=Motion)
    mesh: Mesh = field(default_factory=Mesh)
    step: Step = field(default_factory=Step)

    def __post_init__(self):
        if not self.mesh.covers_tyre(self.tyre):
            raise ValueError(
                f"the fine x band ends at {self.mesh.x_fine_to} m but the tyre "
                f"print reaches {self.tyre.x_centre + self.tyre.radius} m; "
                f"raise Mesh.x_fine_to")
        # The step must run exactly as long as the motion it drives.
        if abs(self.step.dt * round(self.motion.t_final / self.step.dt)
               - self.motion.t_final) > 1e-12:
            raise ValueError(
                f"t_final = {self.motion.t_final} is not a whole number of "
                f"dt = {self.step.dt} increments")

    def summary(self) -> str:
        t = self.tyre
        m = self.motion
        out = [
            f"Case: {self.case_id}",
            self.structure.summary(self.mesh.LY),
            f"  Domain: {self.mesh.LX} x {self.mesh.LY} x {self.mesh.LZ} m "
            f"(x transverse, y depth, z longitudinal)",
            f"  Tyre:   centre x={t.x_centre} m, R={t.radius*1e3:.1f} mm, "
            f"p={t.pressure/1e6:.3f} MPa, F={t.force/1e3:.2f} kN "
            f"(dual spacing {t.dual_spacing*1e3:.0f} mm)",
            f"  Motion: z {m.z_start} -> {m.z_end} m over t {0} -> {m.t_final}",
        ]
        if m.horizontal:
            kind = "acceleration (+z)" if m.mu_sign > 0 else "braking (-z)"
            out.append(f"          traction mu={m.mu} {kind}")
        else:
            out.append("          vertical load only")
        return "\n".join(out)


def make_case(structure: str = "semi_rigid", motion: str = "vertical",
              case_id: str | None = None, **overrides) -> CaseConfig:
    """Build a `CaseConfig` from the named presets.

    >>> cfg = make_case("rigid", "brake")
    >>> cfg.case_id
    'rigid_brake'
    """
    if structure not in STRUCTURES:
        raise ValueError(f"unknown structure {structure!r}; "
                         f"choose from {sorted(STRUCTURES)}")
    if motion not in MOTIONS:
        raise ValueError(f"unknown motion {motion!r}; "
                         f"choose from {sorted(MOTIONS)}")
    return CaseConfig(
        case_id=case_id or f"{structure}_{motion}",
        structure=STRUCTURES[structure],
        motion=MOTIONS[motion],
        **overrides,
    )
