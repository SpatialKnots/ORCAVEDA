from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class HessData:
    filename: str
    atoms: List[str]
    masses: np.ndarray
    coords_A: np.ndarray
    frequencies_cm1: np.ndarray
    ir_intensities: np.ndarray
    normal_modes: np.ndarray
    temperature_K: Optional[float] = None
    frequency_scale_factor: Optional[float] = None
    cartesian_hessian: Optional[np.ndarray] = None


@dataclass
class InternalCoordinate:
    name: str
    kind: str
    atoms0: Tuple[int, ...]
    priority: int
    fn: Callable[[np.ndarray], float]
    source: str = "primitive"


@dataclass
class FunctionalGroup:
    group: str
    atoms0: Tuple[int, ...]
    description: str
    confidence: str
    evidence: str


@dataclass
class AtomEnvironmentAnnotation:
    atom: int
    element: str
    degree: int
    h_neighbors: int
    heavy_neighbors: int
    neighbor_elements: str
    environment_label: str


@dataclass
class ChemicalSystemAnnotation:
    formula: str
    system_type: str
    bonds: Tuple[Tuple[int, int, float], ...]
    fragments: Tuple[Tuple[int, ...], ...]
    functional_groups: Tuple[FunctionalGroup, ...]
    interfragment_hbonds: Tuple[Dict[str, object], ...]

    @property
    def fragment_sizes(self) -> Tuple[int, ...]:
        return tuple(len(fragment) for fragment in self.fragments)

    @property
    def functional_group_labels(self) -> Tuple[str, ...]:
        return tuple(sorted({group.group for group in self.functional_groups}))
