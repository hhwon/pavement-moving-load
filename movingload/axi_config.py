"""
Configuration of the axisymmetric model: materials, mesh, time grid.

An axisymmetric model can only carry a load centred on its axis, so it
cannot move a tyre by itself.  What it can do is compute the *step
response* of a stationary circular print — the strain field R(rho, s) at
every radius rho and every time s after the load is switched on.  Under the
assumptions this package already makes (linear material, quasi-static, layers
homogeneous in plan) that response determines the moving-load response
exactly, by superposition: see `movingload.superpose` and
`docs/axisymmetric.md`.

Coordinates (Abaqus axisymmetric axes)

    1  r  radial     0 -> R_out ;  r = 0 is the axis through the tyre centre
    2  y  axial      0 -> -LY   ;  y = 0 is the road surface, y positive UP
    3  theta  hoop

Unlike the 3-D model, the axial coordinate points *up*.  Abaqus' own
convention, CAE's default view and the sign of E12 all assume it; depths
are still reported as positive numbers in the metadata.

Two asphalt material models:

    elastic       *ELASTIC with the structure's E — one *STATIC step
    viscoelastic  *ELASTIC, MODULI=INSTANTANEOUS + *VISCOELASTIC, TIME=PRONY
                  on the asphalt layers only, then *VISCO creep steps

Everything here is offline: no Abaqus, no filesystem.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

from .config import STRUCTURES, Structure, Tyre


# ----------------------------------------------------------------------
# Linear viscoelastic asphalt
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Viscoelastic:
    """Prony-series relaxation modulus with a constant Poisson ratio.

        E(t) = E0 * [1 - sum_i g_i (1 - exp(-t / tau_i))]

    `E0` is the instantaneous (glassy) modulus and the deck says so
    explicitly with `MODULI=INSTANTANEOUS` — Abaqus' default for `*ELASTIC`
    combined with Prony data is LONG TERM, and silently reading E0 as E_inf
    would stiffen the layer by a factor 1/(1 - sum g).

    A constant Poisson ratio means shear and bulk relax alike, so the
    Abaqus shear and bulk ratios are equal: g_i^P = k_i^P = g_i.
    """
    E0: float
    nu: float
    g: Tuple[float, ...]
    tau: Tuple[float, ...]

    def __post_init__(self):
        if self.E0 <= 0:
            raise ValueError("E0 must be > 0")
        if not (0.0 <= self.nu < 0.5):
            raise ValueError("nu must be in [0, 0.5)")
        if len(self.g) != len(self.tau) or not self.g:
            raise ValueError("g and tau need the same, non-zero length")
        if any(x <= 0 for x in self.g) or any(x <= 0 for x in self.tau):
            raise ValueError("Prony g_i and tau_i must be > 0")
        if sum(self.g) >= 1.0:
            raise ValueError(
                f"sum(g) = {sum(self.g):.4f} >= 1 gives a non-positive "
                f"long-term modulus; Abaqus would reject it")

    @property
    def E_inf(self) -> float:
        return self.E0 * (1.0 - sum(self.g))

    def E_at(self, t: float) -> float:
        """Relaxation modulus at time t [s]."""
        return self.E0 * (1.0 - sum(gi * (1.0 - math.exp(-t / ti))
                                    for gi, ti in zip(self.g, self.tau)))

    @classmethod
    def matched(cls, E_ref: float, t_ref: float, nu: float,
                g: Sequence[float], tau: Sequence[float]) -> "Viscoelastic":
        """Scale a normalised Prony shape so that E(t_ref) = E_ref.

        Used for the presets: the viscoelastic layer then has the elastic
        preset's stiffness at a loading time typical of the passage, and the
        two material models differ by rate dependence alone.
        """
        shape = cls(E0=1.0, nu=nu, g=tuple(g), tau=tuple(tau))
        return cls(E0=E_ref / shape.E_at(t_ref), nu=nu,
                   g=tuple(g), tau=tuple(tau))


#: PLACEHOLDER relaxation shape — smooth, one term per decade from 0.1 ms to
#: 10 s, E0/E_inf = 50.  It is NOT fitted to any mixture or temperature.
#: Replace it with a Prony series fitted to your own dynamic-modulus or
#: relaxation data (at the analysis temperature) before using the numbers.
PLACEHOLDER_PRONY_G = (0.15, 0.25, 0.25, 0.18, 0.10, 0.05)
PLACEHOLDER_PRONY_TAU = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)

#: Loading time at which the placeholder layers match the elastic presets.
#: Roughly the time a point at the bottom of the asphalt is under load at
#: the reference speed (print + spread, ~0.4 m, at 7.6 m/s).
T_MATCH = 0.05

#: Which layers of each preset are asphalt, i.e. which may be viscoelastic.
ASPHALT_LAYERS: Dict[str, Tuple[str, ...]] = {
    "flexible": ("ASPHALT",),
    "semi_rigid": ("ASPHALT", "BINDER"),
    "rigid": ("ASPHALT",),
}


def placeholder_visco(structure: Structure) -> Dict[str, Viscoelastic]:
    """Placeholder viscoelastic models for the asphalt layers of a preset."""
    if structure.name not in ASPHALT_LAYERS:
        raise ValueError(
            f"no asphalt-layer list for structure {structure.name!r}; "
            f"pass AxiCase(visco={{layer_name: Viscoelastic(...)}}) yourself")
    out = {}
    for name in ASPHALT_LAYERS[structure.name]:
        i = structure.names.index(name)
        out[name] = Viscoelastic.matched(structure.E[i], T_MATCH, structure.nu[i],
                                         PLACEHOLDER_PRONY_G, PLACEHOLDER_PRONY_TAU)
    return out


# ----------------------------------------------------------------------
# Mesh
# ----------------------------------------------------------------------

ELEMENT_TYPES = ("CAX4", "CAX4R", "CAX4I")


@dataclass(frozen=True)
class AxiMesh:
    """Radial / axial extent and grid.

    r: uniform `dr_fine` from the axis to the print edge (a node lands
    exactly on r = tyre radius, so the loaded faces are exactly the print),
    uniform again to `r_fine_to`, then `n_coarse` elements graded out to
    `R_out`.  The superposition reads the response at every radius the tyre
    centre passes, so R_out has to exceed the travel half-length — and there
    is no symmetric twin any more, so it is wider than the 3-D model's LX.

    y: each finite layer divided into elements of about `dy_fine`; the
    subgrade gets `n_sub` elements biased towards the surface (as in 3-D).

    CAX4 (full integration) is the default: the mesh is fine enough that
    shear locking is negligible, and there is no hourglass mode to excite at
    the sharp edge of the print.
    """
    R_out: float = 6.0
    LY: float = 5.0
    dr_fine: float = 0.0075
    r_fine_to: float = 1.0
    n_coarse: int = 30
    dy_fine: float = 0.0125
    n_sub: int = 30
    subgrade_bias: float = 1.8
    element: str = "CAX4"

    def __post_init__(self):
        if not (0 < self.r_fine_to < self.R_out):
            raise ValueError("need 0 < r_fine_to < R_out")
        if self.dr_fine <= 0 or self.dy_fine <= 0:
            raise ValueError("element sizes must be > 0")
        if self.n_coarse < 1 or self.n_sub < 1:
            raise ValueError("n_coarse and n_sub must be >= 1")
        if self.subgrade_bias < 1.0:
            raise ValueError("subgrade_bias < 1 would coarsen towards the surface")
        if self.element not in ELEMENT_TYPES:
            raise ValueError(f"element must be one of {ELEMENT_TYPES}")


# ----------------------------------------------------------------------
# Time grid of the creep (step-response) analysis
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class CreepSteps:
    """The step response is needed at every lag between "the load arrived
    here" and "now", from a fraction of a millisecond up to the whole
    passage.  Viscoelastic creep is fastest at the start, so the steps are
    one per decade, each with `n_per_decade` fixed increments:

        CREEP_00   0     -> t_min       n_first increments
        CREEP_01   t_min -> 10 t_min    n_per_decade increments
        ...
        up to t_max

    Fixed increments (no CETOL) make the output times deterministic, which
    the superposition and the tests rely on.  `t_max` bounds the slowest
    speed the post-processing can synthesise: travel length / v <= t_max.
    """
    t_min: float = 1e-5
    t_max: float = 10.0
    n_first: int = 5
    n_per_decade: int = 18

    def __post_init__(self):
        if self.t_min <= 0 or self.t_max <= self.t_min:
            raise ValueError("need 0 < t_min < t_max")
        d = math.log10(self.t_max / self.t_min)
        if abs(d - round(d)) > 1e-9:
            raise ValueError("t_max / t_min must be a whole power of ten")
        if self.n_first < 1 or self.n_per_decade < 1:
            raise ValueError("increment counts must be >= 1")

    @property
    def n_decades(self) -> int:
        return int(round(math.log10(self.t_max / self.t_min)))

    def steps(self) -> List[Tuple[str, float, float]]:
        """(name, fixed increment, step period) for every *VISCO step."""
        out = [("CREEP_00", self.t_min / self.n_first, self.t_min)]
        for d in range(self.n_decades):
            t0 = self.t_min * 10.0 ** d
            period = 9.0 * t0
            out.append((f"CREEP_{d + 1:02d}", period / self.n_per_decade, period))
        return out

    def output_times(self) -> List[float]:
        """Every lag s at which the deck writes the observation row,
        including s = 0 (the *STATIC, instantaneous response)."""
        times = [0.0]
        t = 0.0
        for _, dt, period in self.steps():
            n = int(round(period / dt))
            times.extend(t + dt * (i + 1) for i in range(n))
            t += period
        return times


# ----------------------------------------------------------------------
# Case
# ----------------------------------------------------------------------

MATERIALS = ("elastic", "viscoelastic")


@dataclass(frozen=True)
class AxiCase:
    """One axisymmetric step-response analysis.

    `visco` maps layer names to their viscoelastic model; it must be empty
    for `material="elastic"` and non-empty for `"viscoelastic"`.  Layers not
    listed stay `*ELASTIC` with the structure's E and nu.
    """
    case_id: str
    structure: Structure
    material: str = "elastic"
    visco: Dict[str, Viscoelastic] = field(default_factory=dict)
    tyre: Tyre = field(default_factory=Tyre)
    mesh: AxiMesh = field(default_factory=AxiMesh)
    creep: CreepSteps = field(default_factory=CreepSteps)

    def __post_init__(self):
        if self.material not in MATERIALS:
            raise ValueError(f"material must be one of {MATERIALS}")
        if self.material == "elastic" and self.visco:
            raise ValueError("an elastic case cannot carry viscoelastic layers")
        if self.material == "viscoelastic" and not self.visco:
            raise ValueError("a viscoelastic case needs at least one "
                             "viscoelastic layer in `visco`")
        unknown = set(self.visco) - set(self.structure.names)
        if unknown:
            raise ValueError(f"viscoelastic layers {sorted(unknown)} are not in "
                             f"structure {self.structure.name!r}")
        if self.structure.names[-1] in self.visco:
            raise ValueError("the subgrade cannot be the viscoelastic layer")
        if self.mesh.r_fine_to <= self.tyre.radius:
            raise ValueError("the fine r band must extend past the tyre print")
        self.structure.thicknesses(self.mesh.LY)     # raises if they do not fit

    @property
    def viscoelastic(self) -> bool:
        return self.material == "viscoelastic"

    def summary(self) -> str:
        t = self.tyre
        out = [f"Axisymmetric case: {self.case_id}  ({self.material})",
               f"Structure: {self.structure.name} — {self.structure.description}"]
        s = self.structure
        n_el = [max(1, math.ceil(h / self.mesh.dy_fine - 1e-9)) for h in s.h]
        for nm, h, E, nu, ne in zip(s.names, s.thicknesses(self.mesh.LY), s.E,
                                    s.nu, n_el + [self.mesh.n_sub]):
            out.append(f"  {nm:>10s}: h={h:6.3f} m  E={E/1e9:7.3f} GPa  "
                       f"nu={nu:.2f}  ({ne} element rows)")
        for name, ve in self.visco.items():
            out.append(f"  {name:>10s}: viscoelastic  E0={ve.E0/1e9:.2f} GPa  "
                       f"E({T_MATCH:g} s)={ve.E_at(T_MATCH)/1e9:.2f} GPa  "
                       f"E_inf={ve.E_inf/1e9:.3f} GPa  ({len(ve.g)} Prony terms)")
        out.append(f"  Domain: r 0 -> {self.mesh.R_out} m, depth {self.mesh.LY} m, "
                   f"{self.mesh.element}")
        out.append(f"  Print:  R={t.radius*1e3:.1f} mm, p={t.pressure/1e6:.3f} MPa "
                   f"at r = 0 (step load)")
        if self.viscoelastic:
            out.append(f"  Creep:  {len(self.creep.steps())} *VISCO steps, "
                       f"s = {self.creep.t_min:g} .. {self.creep.t_max:g} s")
        return "\n".join(out)


def make_axi_case(structure: str = "semi_rigid", material: str = "elastic",
                  case_id: str | None = None, **overrides) -> AxiCase:
    """Build an `AxiCase` from the 3-D structure presets.

    >>> make_axi_case("semi_rigid", "viscoelastic").case_id
    'axi_semi_rigid_viscoelastic'
    """
    if structure not in STRUCTURES:
        raise ValueError(f"unknown structure {structure!r}; "
                         f"choose from {sorted(STRUCTURES)}")
    s = STRUCTURES[structure]
    if "visco" not in overrides and material == "viscoelastic":
        overrides["visco"] = placeholder_visco(s)
    return AxiCase(case_id=case_id or f"axi_{structure}_{material}",
                   structure=s, material=material, **overrides)
