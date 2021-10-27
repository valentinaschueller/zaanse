import numpy as np

from .simulation_params import SimulationParams

def compute_pressure(params: SimulationParams, z: np.ndarray):
    # using https://en.wikipedia.org/wiki/Atmospheric_pressure#Altitude_variation
    p0 = params.p0
    cp = params.cp
    T0 = params.T0
    M = params.M
    R0 = params.R0
    g = params.g
    p = p0 * (1 - g * z / (cp * T0)) ** (cp * M / R0)
    return p


def compute_temperature(
    params: SimulationParams, theta: np.ndarray, pressure: np.ndarray
):
    # using https://en.wikipedia.org/wiki/Potential_temperature
    temperature = theta * ((pressure / params.p0) ** (0.286))
    return temperature

def P_2(x):
    return 1.5 * x ** 2 - 0.5


def P_2_sin(x):
    return P_2(np.sin(x))


def theta_E(z: float, phi: float, params):
    return params.theta_0 * (
        1
        - 2 / 3 * params.delta_h * P_2_sin(phi)
        + params.delta_v * (z / params.H - 0.5)
    )


def geopotential(z: float, g: float):
    return g * z


def explicit_euler_step(var: np.array, dt: float, rhs: np.array):
    var_next = var + dt * rhs
    return var_next
