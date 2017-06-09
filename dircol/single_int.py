#! /usr/bin/python

from collections import OrderedDict

import numpy as np
import sympy as sym
from opty.direct_collocation import Problem
from opty.utils import building_docs
import matplotlib.pyplot as plt
import matplotlib.animation as animation

duration = 10.0
num_nodes = 500

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
t = sym.symbols('t')
x, y, ux, uy = sym.symbols('x, y, ux, uy', cls=sym.Function)

state_symbols = (x(t), y(t))
specified_symbols = (ux(t), uy(t))

eom = sym.Matrix([x(t).diff() - ux(t),
                  y(t).diff() - uy(t)])

# Specify the objective function and it's gradient.

def obj(free):
    """Minimize the sum of the squares of the control input."""
    u = free[2 * num_nodes:]
    return interval_value * np.sum(u**2)

def obj_grad(free):
    grad = np.zeros_like(free)
    grad[2 * num_nodes:] = 2.0 * interval_value * free[2 * num_nodes:]
    return grad

# Specify the symbolic instance constraints, i.e. initial and end
# conditions.
instance_constraints = (x(0.0),
                        x(duration) - 1.0,
                        y(0.0),
                        y(duration))

# Create an optimization problem.
prob = Problem(obj, obj_grad, eom, state_symbols, num_nodes, interval_value,
               instance_constraints=instance_constraints,
               bounds={ux(t): (-10, 10), uy(t): (-10, 10), y(t): (0, 0)})

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

# Make some plots
prob.plot_trajectories(solution)

plt.show()
