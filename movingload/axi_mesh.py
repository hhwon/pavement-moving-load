"""
Structured quadrilateral grid for the axisymmetric model.

Node and element numbering is i-fastest (radial), then j (depth):

    node id = 1 + i + j*(NR+1)
    elem id = 1 + i + j*NR

with j = 0 at the road surface.  Nodes are written with the axial
coordinate y = -depth (positive up), and each element is connected
counter-clockwise in the (r, y) plane:

    n4 (i, j)   ---- n3 (i+1, j)        <- face S3: the upper face
     |                 |
    n1 (i, j+1) ---- n2 (i+1, j+1)

so the tyre pressure goes on face S3 of the surface row.  The extractor
finds the observation row through the ids in the metadata sidecar;
changing this numbering invalidates existing result files.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import List

import numpy as np

from .axi_config import AxiCase, AxiMesh
from .config import Structure, Tyre
from .mesh import build_y


@dataclass(frozen=True)
class AxiGrid:
    r: np.ndarray                # radial node coordinates, length NR+1
    depth: np.ndarray            # node depths (positive down), length NY+1
    layer_of_row: List[int]
    n_print: int                 # elements with r <= tyre radius in the top row
    j_obs: int                   # observation row (bottom row of layer 0)

    @property
    def NR(self) -> int:
        return len(self.r) - 1

    @property
    def NY(self) -> int:
        return len(self.depth) - 1

    @property
    def n_elements(self) -> int:
        return self.NR * self.NY

    @property
    def n_nodes(self) -> int:
        return (self.NR + 1) * (self.NY + 1)

    def node_id(self, i: int, j: int) -> int:
        return 1 + i + j * (self.NR + 1)

    def elem_id(self, i: int, j: int) -> int:
        return 1 + i + j * self.NR

    def connectivity(self, i: int, j: int) -> List[int]:
        """Counter-clockwise in (r, y = -depth): see the module docstring."""
        return [self.node_id(i, j + 1), self.node_id(i + 1, j + 1),
                self.node_id(i + 1, j), self.node_id(i, j)]

    @property
    def obs_elem_ids(self) -> List[int]:
        return [self.elem_id(i, self.j_obs) for i in range(self.NR)]

    @property
    def r_centres(self) -> List[float]:
        return [float(0.5 * (self.r[i] + self.r[i + 1])) for i in range(self.NR)]

    @property
    def obs_depth(self) -> float:
        """Centroid depth of the observation row."""
        return float(0.5 * (self.depth[self.j_obs] + self.depth[self.j_obs + 1]))


def build_r(mesh: AxiMesh, tyre: Tyre) -> np.ndarray:
    """Radial nodes: exact node at the print edge, uniform fine band, graded
    coarse band.  The coarse band is a power law whose first element is
    about twice `dr_fine`, so the size jumps smoothly."""
    a = tyre.radius
    n_a = max(1, math.ceil(a / mesh.dr_fine - 1e-9))
    n_b = max(1, math.ceil((mesh.r_fine_to - a) / mesh.dr_fine - 1e-9))
    inner = np.linspace(0.0, a, n_a + 1)
    band = np.linspace(a, mesh.r_fine_to, n_b + 1)[1:]

    span = mesh.R_out - mesh.r_fine_to
    first = min(2.0 * mesh.dr_fine, span / mesh.n_coarse)
    p = max(1.0, math.log(first / span) / math.log(1.0 / mesh.n_coarse))
    xi = np.linspace(0.0, 1.0, mesh.n_coarse + 1)[1:]
    coarse = mesh.r_fine_to + span * xi ** p
    return np.concatenate((inner, band, coarse))


def axi_structure(structure: Structure, mesh: AxiMesh) -> Structure:
    """The structure with its element counts replaced by the axisymmetric
    mesh density — the 2-D model is cheap enough to be much finer than 3-D."""
    n = [max(1, math.ceil(h / mesh.dy_fine - 1e-9)) for h in structure.h]
    return replace(structure, n_elems=n + [mesh.n_sub])


def build_axi_grid(case: AxiCase) -> AxiGrid:
    r = build_r(case.mesh, case.tyre)
    s = axi_structure(case.structure, case.mesh)
    depth, layer_of_row = build_y(s, case.mesh)
    n_print = int(np.sum(r[1:] <= case.tyre.radius + 1e-12))
    if abs(r[n_print] - case.tyre.radius) > 1e-12:
        raise AssertionError("no node on the print edge")
    return AxiGrid(r=r, depth=depth, layer_of_row=layer_of_row,
                   n_print=n_print, j_obs=s.n_elems[0] - 1)


def describe(grid: AxiGrid) -> str:
    dr, dy = np.diff(grid.r), np.diff(grid.depth)
    return "\n".join([
        f"mesh: NR={grid.NR}, NY={grid.NY}, elements={grid.n_elements:,}, "
        f"nodes={grid.n_nodes:,}",
        f"      dr {dr.min():.4f} .. {dr.max():.4f} m, "
        f"dy {dy.min():.4f} .. {dy.max():.4f} m",
        f"      print: {grid.n_print} faces, edge node at r = {grid.r[grid.n_print]:.4f} m",
        f"      observation row j={grid.j_obs} at depth {grid.obs_depth:.4f} m",
    ])
