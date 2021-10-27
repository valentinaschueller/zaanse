from .helpers import compute_pressure, compute_temperature
from .simulation_params import SimulationParams
from .visualization import plot_scalar_field
from .dynamics import NSESolver


class SimulationInstance:
    def __init__(self, args: dict = {}):
        self.params = SimulationParams()
        if args:
            self.visualize = args["visualize"]
            self.params.dt = args["t_step"]
            self.params.t_end = args["t_end"]
            self.params.n_z = args["nz"]
            self.params.n_phi = args["nphi"]
        else:
            self.visualize = False
        self.solver = NSESolver(self.params)
        self.time = 0.0
        self.pressure = compute_pressure(self.params, self.solver.z_central[1:-1, 1:-1])

    def run_simulation(self):
        while self.time < self.params.t_end:
            self.simulate_timestep()
        # self.check_divergence()
        if self.visualize:
            plot_scalar_field(self.temperature, title=f"T-{round(self.time)}")
            plot_scalar_field(self.solver.u_central, title=f"u-{round(self.time)}")
            plot_scalar_field(self.solver.v_vertical, title=f"v-{round(self.time)}")
            plot_scalar_field(self.solver.w_horizontal, title=f"w-{round(self.time)}")

    def simulate_timestep(self):
        self.solver.do_timestep(self.params)
        self.time = self.time + self.params.dt
        self.temperature = compute_temperature(
            self.params, self.solver.theta_central[1:-1, 1:-1], self.pressure
        )

    def check_divergence(self):
        self.solver.v_vertical[::2] - self.solver.v_vertical[1::2] + self.solver.w_horizontal
