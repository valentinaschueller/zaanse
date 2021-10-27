import numpy as np

from zaanse.helpers import explicit_euler_step


def test_explicit_euler_step():
    var = np.array([1, 2, 3])
    dt = 0.1
    rhs = np.array([2, -5, 6.0])
    result = explicit_euler_step(var, dt, rhs)
    expected = np.array([1.2, 1.5, 3.6])
    assert (expected == result).all()
