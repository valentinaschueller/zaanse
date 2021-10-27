import numpy as np

from .helpers import geopotential, theta_E, explicit_euler_step


class NSESolver:
    def __init__(self, params):
        self.create_grids(params)
        self.theta_eq_central = theta_E(self.z_central, self.phi_central, params)
        self.f_vertical = 2 * params.omega * np.sin(self.phi_vertical)
        self.f_central = 2 * params.omega * np.sin(self.phi_central)
        self.gz = geopotential(self.z_horizontal, params.g)

        # where the actual dynamics happens
        self.u_central = np.random.rand(params.n_z + 2, params.n_phi + 2)
        self.v_vertical = np.zeros((params.n_z + 2, params.n_phi + 1))
        self.w_horizontal = np.zeros((params.n_z + 1, params.n_phi + 2))
        self.theta_central = self.theta_eq_central.copy()

        # useful interpolated grids
        self.u_vertical = np.zeros((params.n_z, params.n_phi - 1))
        self.v_central = np.zeros((params.n_z, params.n_phi))
        self.v_horizontal = np.zeros((params.n_z - 1, params.n_phi))
        self.w_central = np.zeros((params.n_z, params.n_phi))
        self.w_vertical = np.zeros((params.n_z, params.n_phi - 1))
        self.theta_horizontal = np.zeros((params.n_z - 1, params.n_phi))

        self.u_next = self.u_central[1:-1, 1:-1].copy()
        self.v_next = self.v_vertical[1:-1, 1:-1].copy()
        self.w_next = self.w_horizontal[1:-1, 1:-1].copy()
        self.theta_next = self.theta_central[1:-1, 1:-1].copy()

        self.apply_boundary_conditions(params)

    def create_grids(self, params):
        self.phi_central = np.outer(
            np.ones(params.n_z + 2),
            np.linspace(
                -np.pi / 2 - 0.5 * params.dphi,
                np.pi / 2 + 0.5 * params.dphi,
                params.n_phi + 2,
            ),
        )  # (n_z + 2) * (n_phi + 2)
        self.phi_vertical = np.outer(
            np.ones(params.n_z + 2),
            np.linspace(-np.pi / 2, np.pi / 2, params.n_phi + 1),
        )  # (n_z + 2) * (n_phi + 1)
        self.z_central = np.outer(
            np.linspace(-0.5 * params.dz, params.H + 0.5 * params.dz, params.n_z + 2),
            np.ones(params.n_phi + 2),
        )  # (n_z + 2) * (n_phi + 2)
        self.z_horizontal = np.outer(
            np.linspace(0, params.H, params.n_z + 1), np.ones(params.n_phi + 2)
        )  # (n_z + 1) * (n_phi + 2)

    def initialize_fields(self):
        # initializing
        # ...
        # self.apply_boundary_conditions()
        pass

    def interpolate_fields(self):
        self.u_vertical = 0.5 * (
            self.u_central[1:-1, 1:-2] + self.u_central[1:-1, 2:-1]
        )
        self.v_central = 0.5 * (self.v_vertical[1:-1, 0:-1] + self.v_vertical[1:-1, 1:])
        self.v_horizontal = 0.5 * (self.v_central[0:-1, :] + self.v_central[1:, :])
        self.w_central = 0.5 * (
            self.w_horizontal[0:-1, 1:-1] + self.w_horizontal[1:, 1:-1]
        )
        self.w_vertical = 0.5 * (self.w_central[:, 0:-1] + self.w_central[:, 1:])
        self.theta_horizontal = 0.5 * (
            self.theta_central[1:-2, 1:-1] + self.theta_central[2:-1, 1:-1]
        )

    def apply_boundary_conditions(self, params):
        # Dirichlet at poles for u, v, w, theta:
        self.u_central[:, [0, -1]] = -self.u_central[:, [1, -2]].copy()
        self.v_vertical[:, [0, -1]] = 0.0
        self.w_horizontal[:, [0, -1]] = -self.w_horizontal[:, [1, -2]].copy()
        theta_eq_boundary = theta_E(
            self.z_central[:, [0, -1]], self.phi_central[:, [0, -1]], params
        )
        self.theta_central[:, [0, -1]] = (
            2 * theta_eq_boundary - self.theta_central[:, [1, -2]].copy()
        )
        # Dirichlet at top/bottom for w:
        self.w_horizontal[[0, -1], :] = 0.0
        # von Neumann top for u
        self.u_central[-1, 1:-1] = self.u_central[-2, 1:-1].copy()
        # von Neumann top for v
        self.v_vertical[-1, 1:-1] = self.v_vertical[-2, 1:-1].copy()
        # von Neumann bottom for u
        self.u_central[0, 1:-1] = (
            (params.nu / params.dz)
            * (1 / (params.C + params.nu / params.dz))
            * self.u_central[1, 1:-1]
        )
        # von Neumann bottom for v
        self.v_vertical[0, 1:-1] = (
            (params.nu / params.dz)
            * (1 / (params.C + params.nu / params.dz))
            * self.v_vertical[1, 1:-1]
        )
        # von Neumann top/bottom for Theta
        self.theta_central[[0, -1], 1:-1] = self.theta_central[[1, -2], 1:-1]

    def do_timestep(self, params):
        self.interpolate_fields()
        self.compute_u_next(params)
        self.compute_v_next(params)
        self.compute_w_next(params)
        self.compute_theta_next(params)
        self.update_fields()
        # print(f"u max/min/mean: {np.max(self.u)}, {np.min(self.u)}, {np.mean(self.u)}")
        # print(f"v max/min/mean: {np.max(self.v)}, {np.min(self.v)}, {np.mean(self.v)}")
        # print(f"w max/min/mean: {np.max(self.w)}, {np.min(self.w)}, {np.mean(self.w)}")
        # print(f"theta max/min/mean: {np.max(self.theta)}, {np.min(self.theta)}, {np.mean(self.theta)}")
        self.apply_boundary_conditions(params)

    def compute_u_next(self, params):
        rhs_u = self.eval_fu(params)
        # print(f"rhs_u max/min/mean: {np.max(rhs_u)}, {np.min(rhs_u)}, {np.mean(rhs_u)}")
        self.u_next = explicit_euler_step(self.u_central[1:-1, 1:-1], params.dt, rhs_u)

    def compute_v_next(self, params):
        rhs_v = self.eval_fv(params)
        # print(f"rhs_v max/min/mean: {np.max(rhs_v)}, {np.min(rhs_v)}, {np.mean(rhs_v)}")
        self.v_next = explicit_euler_step(self.v_vertical[1:-1, 1:-1], params.dt, rhs_v)

    def compute_w_next(self, params):
        rhs_w = self.eval_fw(params)
        # print(f"rhs_w max/min/mean: {np.max(rhs_w)}, {np.min(rhs_w)}, {np.mean(rhs_w)}")
        self.w_next = explicit_euler_step(
            self.w_horizontal[1:-1, 1:-1], params.dt, rhs_w
        )

    def compute_theta_next(self, params):
        rhs_theta = self.eval_ftheta(params)
        # print(f"rhs_theta max/min/mean: {np.max(rhs_theta)}, {np.min(rhs_theta)}, {np.mean(rhs_theta)}")
        self.theta_next = explicit_euler_step(
            self.theta_central[1:-1, 1:-1], params.dt, rhs_theta
        )

    def update_fields(self):
        """Copy values from self.var_next to self.var"""
        self.u_central[1:-1, 1:-1] = self.u_next[:].copy()
        self.v_vertical[1:-1, 1:-1] = self.v_next[:].copy()
        self.w_horizontal[1:-1, 1:-1] = self.w_next[:].copy()
        self.theta_central[1:-1, 1:-1] = self.theta_next[:].copy()

    def eval_fu(self, params):
        u_inner = self.u_central[1:-1, 1:-1]
        phi_inner = self.phi_central[1:-1, 1:-1]
        v_advection = (
            -1
            / params.a
            * self.v_central
            * self.ddphi(self.u_central, self.v_central)
            / params.dphi
        )
        w_advection = (
            -self.w_central * self.ddz(self.u_central, self.w_central) / params.dz
        )
        coriolis = self.f_central[1:-1, 1:-1] * self.v_central
        angular_momentum = 1 / params.a * self.v_central * u_inner * np.tan(phi_inner)
        z_mixing = params.nu * self.ddzz(self.u_central) / params.dz ** 2
        return v_advection + w_advection + coriolis + angular_momentum + z_mixing

    def eval_fv(self, params):
        v_inner = self.v_vertical[1:-1, 1:-1]
        phi_inner = self.phi_vertical[1:-1, 1:-1]
        v_advection = (
            -1 / params.a * v_inner * self.ddphi(self.v_vertical, v_inner) / params.dphi
        )
        w_advection = (
            -self.w_vertical * self.ddz(self.v_vertical, self.w_vertical) / params.dz
        )
        coriolis = -self.f_vertical[1:-1, 1:-1] * self.u_vertical
        angular_momentum = (
            -1 / params.a * self.u_vertical * self.u_vertical * np.tan(phi_inner)
        )
        z_mixing = params.nu * self.ddzz(self.v_vertical) / params.dz ** 2
        return v_advection + w_advection + coriolis + angular_momentum + z_mixing

    def eval_fw(self, params):
        v_advection = (
            (-1 / params.a)
            * self.v_horizontal
            * self.ddphi(self.w_horizontal, self.v_horizontal)
            / params.dphi
        )
        p_gradient = -self.ddz(self.gz, self.w_horizontal[1:-1, 1:-1]) / params.dz
        # print(f"p_gradient max/min/mean: {np.max(p_gradient)}, {np.min(p_gradient)}, {np.mean(p_gradient)}")
        buoyancy = params.g * self.theta_horizontal / params.theta_0
        # print(f"buoyancy max/min/mean: {np.max(buoyancy)}, {np.min(buoyancy)}, {np.mean(buoyancy)}")
        return v_advection + p_gradient + buoyancy

    def eval_ftheta(self, params):
        theta_inner = self.theta_central[1:-1, 1:-1]
        v_advection = (
            (-1 / params.a)
            * self.v_central
            * self.ddphi(self.theta_central, self.v_central)
            / params.dphi
        )
        w_advection = (
            -self.w_central * self.ddz(self.theta_central, self.w_central) / params.dz
        )
        heating_cooling = (
            -1 / params.tau * (theta_inner - self.theta_eq_central[1:-1, 1:-1])
        )
        z_mixing = params.nu * self.ddzz(self.theta_central) / params.dz ** 2
        return v_advection + w_advection + heating_cooling + z_mixing

    def ddz(self, array: np.ndarray, velocity_field: np.ndarray):
        bit_array = velocity_field > 0
        result = (
            array[1:-1, 1:-1] * (2 * bit_array - 1)
            - array[0:-2, 1:-1] * bit_array
            + array[2:, 1:-1] * bit_array
        )
        return result

    def ddphi(self, array: np.ndarray, velocity_field: np.ndarray):
        bit_array = velocity_field > 0
        result = (
            array[1:-1, 1:-1] * (2 * bit_array - 1)
            - array[1:-1, 0:-2] * bit_array
            + array[1:-1, 2:] * bit_array
        )
        return result

    def ddzz(self, array: np.ndarray):
        """numerator of a three-point stencil to approximate d^2/dz^2"""
        return array[0:-2, 1:-1] - 2 * array[1:-1, 1:-1] + array[2:, 1:-1]
