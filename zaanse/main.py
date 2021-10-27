import argparse

from .simulation import SimulationInstance

def main():
    args = parse_arguments()
    sim = SimulationInstance(args)
    sim.run_simulation()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--visualize",
        required=False,
        action="store_true",
        help="Execute the project with a visualization for debugging",
    )
    parser.add_argument(
        "-T",
        "--t_end",
        required=False,
        default=60.0,
        type=float,
        help="Time until the simulation should run",
    )
    parser.add_argument(
        "-dt",
        "--t_step",
        type=float,
        required=False,
        default=5.0,
        help="Time step size (constant)",
    )
    parser.add_argument(
        "-nphi",
        type=int,
        required=False,
        default=90,
        help="Number of cells in phi direction",
    )
    parser.add_argument(
        "-nz",
        type=int,
        required=False,
        default=50,
        help="Number of cells in z direction",
    )
    # add other arguments here if needed
    args = vars(parser.parse_args())
    return args
