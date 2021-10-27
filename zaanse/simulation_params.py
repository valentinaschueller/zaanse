from dataclasses import dataclass
import numpy as np


@dataclass
class SimulationParams:
    """Class to store the simulation parameters"""

    # physics-related parameters
    nu: float = 2.5
    omega: float = 2 * np.pi / (8.64e6)
    a: float = 6.4e6
    g: float = 9.8
    H: float = 8.0e3
    delta_h: float = 1 / 3
    delta_v: float = 1 / 8
    C: float = 0.005
    tau: float = 20 * 8.64e4
    theta_0: float = 287

    # physics-related parameters which should probably not be touched
    p0 = 101325
    cp = 1004.68506
    T0 = 288.16
    M = 0.02896968
    R0 = 8.314462618

    # simulation-related parameters
    dt: float = 2.0
    t_end: float = 60.0
    n_z: int = 50
    n_phi: int = 90
    dphi: float = np.pi / n_phi
    dz: float = H / n_z
    pole_island_factor = 0.1