"""
Structured hexahedral grid for the half-pavement model.

Separated from the deck writer so the grid can be inspected, plotted or
unit-tested without producing an `.inp`.  Node and element numbering is
i-fastest, then j, then k:

    node id  = 1 + i + j*(NX+1) + k*(NX+1)*(NY+1)
    elem id  = 1 + i + j*NX     + k*NX*NY

with i along x (transverse), j along y (depth), k along z (longitudinal).
The extractor relies on this ordering through the element ids recorded in
the metadata sidecar, so changing it invalidates existing result files.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np

from .config import CaseConfig, Mesh, Structure, Tyre


@dataclass(frozen=True)
class Grid:
    """The generated grid plus the bookkeeping the postprocessor needs."""
    x: np.ndarray                # node coordinates along x, length NX+1
    y: np.ndarray                # node coordinates along y, length NY+1
    z: np.ndarray                # node coordinates along z, length NZ+1
    layer_of_row: List[int]      # layer index of each of the NY element rows
    i_obs: int                   # x-column of the observation elements
    j_obs: int                   # y-row  of the observation elements

    @property
    def NX(self) -> int:
        return len(self.x) - 1

    @property
    def NY(self) -> int:
        return len(self.y) - 1

    @property
    def NZ(self) -> int:
        return len(self.z) - 1

    @property
    def n_elements(self) -> int:
        return self.NX * self.NY * self.NZ

    @property
    def n_nodes(self) -> int:
        return (self.NX + 1) * (self.NY + 1) * (self.NZ + 1)

    # -- id helpers -----------------------------------------------------
    def node_id(self, i: int, j: int, k: int) -> int:
        return 1 + i + j * (self.NX + 1) + k * (self.NX + 1) * (self.NY + 1)

    def elem_id(self, i: int, j: int, k: int) -> int:
        return 1 + i + j * self.NX + k * self.NX * self.NY

    # -- observation ----------------------------------------------------
    @property
    def obs_elem_ids(self) -> List[int]:
        """Bottom row of the surface layer at the tyre column, all along z."""
        return [self.elem_id(self.i_obs, self.j_obs, k) for k in range(self.NZ)]

    @property
    def obs_x(self) -> float:
        """Centroid x of the observation column."""
        return float(0.5 * (self.x[self.i_obs] + self.x[self.i_obs + 1]))

    @property
    def z_centres(self) -> List[float]:
        return [float(0.5 * (self.z[k] + self.z[k + 1])) for k in range(self.NZ)]


def build_x(mesh: Mesh) -> np.ndarray:
    """Three-zone transverse grid: uniform fine, medium, coarse."""
    fine = np.arange(0.0, mesh.x_fine_to + 1e-9, mesh.dx_fine)
    med = np.linspace(mesh.x_fine_to, mesh.x_med_to, mesh.n_med + 1)[1:]
    coarse = np.linspace(mesh.x_med_to, mesh.LX, mesh.n_coarse + 1)[1:]
    return np.concatenate((fine, med, coarse))


def build_z(mesh: Mesh) -> np.ndarray:
    """Uniform longitudinal grid — the load sweeps through, so no zone is
    privileged and a uniform dz keeps the sampled strain history evenly spaced."""
    nz = int(round(mesh.LZ / mesh.dz))
    return np.linspace(0.0, mesh.LZ, nz + 1)


def build_y(structure: Structure, mesh: Mesh) -> tuple[np.ndarray, List[int]]:
    """Depth grid: uniform within each finite layer, biased in the subgrade.

    Returns the node coordinates and, for each element row, its layer index.
    """
    thick = structure.thicknesses(mesh.LY)
    nodes: List[float] = [0.0]
    for h, n in zip(thick[:-1], structure.n_elems[:-1]):
        nodes.extend((nodes[-1] + np.linspace(0, h, n + 1)[1:]).tolist())

    # Subgrade: pack elements towards the surface, where the gradients are.
    n_sub = structure.n_elems[-1]
    xi = np.linspace(0, 1, n_sub + 1)[1:]
    nodes.extend((nodes[-1] + thick[-1] * (xi ** mesh.subgrade_bias)).tolist())

    layer_of_row: List[int] = []
    for ilayer, n in enumerate(structure.n_elems):
        layer_of_row.extend([ilayer] * n)
    return np.array(nodes), layer_of_row


def build_grid(cfg: CaseConfig) -> Grid:
    """Assemble the full grid for a case and locate the observation column."""
    x = build_x(cfg.mesh)
    z = build_z(cfg.mesh)
    y, layer_of_row = build_y(cfg.structure, cfg.mesh)

    # Observation column: the x-element whose centroid is nearest the tyre
    # centre; observation row: the bottom element row of the surface layer.
    x_centres = 0.5 * (x[:-1] + x[1:])
    i_obs = int(np.argmin(np.abs(x_centres - cfg.tyre.x_centre)))
    j_obs = cfg.structure.n_elems[0] - 1

    if len(layer_of_row) != len(y) - 1:
        raise AssertionError("layer bookkeeping does not match the y grid")
    return Grid(x=x, y=y, z=z, layer_of_row=layer_of_row,
                i_obs=i_obs, j_obs=j_obs)


def describe(grid: Grid) -> str:
    dx = np.diff(grid.x)
    dy = np.diff(grid.y)
    dz = np.diff(grid.z)
    return "\n".join([
        f"mesh: NX={grid.NX}, NY={grid.NY}, NZ={grid.NZ}, "
        f"elements={grid.n_elements:,}, nodes={grid.n_nodes:,}",
        f"      dx {dx.min():.4f} .. {dx.max():.4f} m",
        f"      dy {dy.min():.4f} .. {dy.max():.4f} m",
        f"      dz {dz.min():.4f} .. {dz.max():.4f} m (uniform)",
        f"      observation column i={grid.i_obs} at x={grid.obs_x:.4f} m, "
        f"row j={grid.j_obs}",
    ])
