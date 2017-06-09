#! /usr/bin/python

import numpy as np
import sympy as sym
from opty.direct_collocation import Problem
import matplotlib.pyplot as plt

duration = 10.0
num_nodes = 500

interval_value = duration / (num_nodes - 1)

# Symbolic equations of motion
t = sym.symbols('t')
x, y, theta, vpx, vpy, L1, L2 = sym.symbols('x, y, theta, vpx, vpy, L1, L2', cls=sym.Function)

state_symbols = (x(t), y(t), theta(t))
specified_symbols = (vpx(t), vpy(t))

def get_pushing_eom(rox1_0, roy1_0):
    c = 1
    r1sq = rox1_0**2 + roy1_0**2
    rox1 = (sym.cos(theta(t))*rox1_0 - sym.sin(theta(t))*roy1_0)
    roy1 = (sym.sin(theta(t))*rox1_0 + sym.cos(theta(t))*roy1_0)
    rox1roy1 = (rox1_0*roy1_0*sym.cos(2*theta(t)) + 0.5*(rox1_0*rox1_0 - roy1_0*roy1_0)*sym.sin(2*theta(t)))
    rox1sq = (0.5*(r1sq + (rox1_0*rox1_0 - roy1_0*roy1_0)*sym.cos(2*theta(t)) - 2*rox1_0*roy1_0*sym.sin(2*theta(t))))
    roy1sq = (0.5*(r1sq + (roy1_0*roy1_0 - rox1_0*rox1_0)*sym.cos(2*theta(t)) + 2*rox1_0*roy1_0*sym.sin(2*theta(t))))
    eom1 = sym.Matrix([x(t).diff() - ((c*c + rox1sq)*vpx(t) + rox1roy1*vpy(t))/(c*c + r1sq),
                      y(t).diff() - (rox1roy1*vpx(t) + (c*c + roy1sq)*vpy(t))/(c*c + r1sq),
                      theta(t).diff() - (rox1*vpy(t) - roy1*vpx(t))/(c*c + r1sq)])
    return eom1

eom1 = get_pushing_eom(-0.5, 0.0)
eom2 = get_pushing_eom(0.0, -0.5)

eom = L1(t)*eom1 + L2(t)*eom2

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
                        y(duration) - 1.0,
                        theta(0.0),
                        theta(duration))

bounds = {
    vpx(t): (-10, 10),
    vpy(t): (-10, 10),
}

known_trajectory_map = {
    L1(t): np.concatenate([np.ones(num_nodes/2),np.zeros(num_nodes/2)]),
    L2(t): np.concatenate([np.zeros(num_nodes/2),np.ones(num_nodes/2)])
}

# Create an optimization problem.
prob = Problem(obj, obj_grad, eom, state_symbols, num_nodes, interval_value,
               instance_constraints=instance_constraints,
               bounds=bounds,
               known_trajectory_map=known_trajectory_map)

# Use a random positive initial guess.
initial_guess = np.random.randn(prob.num_free)

# Find the optimal solution.
solution, info = prob.solve(initial_guess)

# Make some plots
prob.plot_trajectories(solution)

plt.show()
