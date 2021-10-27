import numpy as np

from zaanse.simulation import SimulationInstance


def test_simulation_creation():
    sim = SimulationInstance()

    # Test grid shapes
    assert sim.solver.z_central.shape == (sim.params.n_z + 2, sim.params.n_phi + 2)
    assert sim.solver.phi_central.shape == (sim.params.n_z + 2, sim.params.n_phi + 2)
    assert sim.solver.z_horizontal.shape == (sim.params.n_z + 1, sim.params.n_phi + 2)
    assert sim.solver.phi_vertical.shape == (sim.params.n_z + 2, sim.params.n_phi + 1)

    assert sim.solver.u_central.shape == sim.solver.z_central.shape
    assert sim.solver.v_vertical.shape == sim.solver.phi_vertical.shape
    assert sim.solver.w_horizontal.shape == sim.solver.z_horizontal.shape
    assert sim.solver.theta_central.shape == sim.solver.z_central.shape

    assert sim.solver.u_vertical.shape == tuple(
        dim - 2 for dim in sim.solver.v_vertical.shape
    )
    assert sim.solver.v_central.shape == tuple(
        dim - 2 for dim in sim.solver.u_central.shape
    )
    assert sim.solver.v_horizontal.shape == tuple(
        dim - 2 for dim in sim.solver.w_horizontal.shape
    )
    assert sim.solver.w_central.shape == tuple(
        dim - 2 for dim in sim.solver.u_central.shape
    )
    assert sim.solver.w_vertical.shape == tuple(
        dim - 2 for dim in sim.solver.v_vertical.shape
    )
    assert sim.solver.theta_horizontal.shape == tuple(
        dim - 2 for dim in sim.solver.w_horizontal.shape
    )


def test_does_timestep():
    sim = SimulationInstance()
    sim.simulate_timestep()
    assert sim.time == sim.params.dt


def test_boundary_conditions():
    sim = SimulationInstance()
    dim_z_central, dim_phi_central = sim.solver.u_central.shape
    dim_z_vertical, dim_phi_vertical = sim.solver.v_vertical.shape
    dim_z_horizontal, dim_phi_horizontal = sim.solver.w_horizontal.shape
    sim.solver.u_central = np.random.rand(dim_z_central, dim_phi_central)
    sim.solver.v_vertical = np.random.rand(dim_z_vertical, dim_phi_vertical)
    sim.solver.w_horizontal = np.random.rand(dim_z_horizontal, dim_phi_horizontal)
    sim.solver.theta_central = np.random.rand(dim_z_central, dim_phi_central)
    sim.solver.apply_boundary_conditions(sim.params)
    assert (
        np.sum(
            sim.solver.u_central[1:-1, [0, -1]] + sim.solver.u_central[1:-1, [1, -2]]
        )
        == 0
    )
    assert np.sum(np.abs(sim.solver.v_vertical[:, [0, -1]])) == 0
    assert (
        np.sum(
            sim.solver.w_horizontal[:, [0, -1]] + sim.solver.w_horizontal[:, [1, -2]]
        )
        == 0
    )
    assert np.sum(np.abs(sim.solver.w_horizontal[[0, -1], :])) == 0


def test_simulation():
    sim = SimulationInstance()
    sim.run_simulation()
    assert sim.time >= sim.params.t_end
    assert (sim.time - sim.params.dt) < sim.params.t_end
