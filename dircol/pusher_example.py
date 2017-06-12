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
x, y, theta, vpx, vpy, L2, L3, L4, L5 = sym.symbols('x, y, theta, vpx, vpy, L2, L3, L4, L5', cls=sym.Function)

state_symbols = (x(t), y(t), theta(t))
specified_symbols = (vpx(t), vpy(t))

def get_pushing_eom(rox1_0, roy1_0):
    c = 0.38259785823210635
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

eom2 = get_pushing_eom(0.0, -0.5)
eom3 = get_pushing_eom(0.5, 0.0)
eom4 = get_pushing_eom(0.0, 0.5)
eom5 = get_pushing_eom(-0.5, 0.0)

eom = L2(t)*eom2 + L3(t)*eom3 + L4(t)*eom4 + L5(t)*eom5

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
                        x(duration) - 3.0,
                        y(0.0),
                        y(duration) - 1.0,
                        theta(0.0),
                        theta(duration))

bounds = {
    vpx(t): (-1, 1),
    vpy(t): (-1, 1),

    ((vpx(t)*1.0 + vpy(t)*0.0)**2 - 0.030153689607045803*(vpx(t)**2 + vpy(t)**2)*1.0)*L2(t): (-float("inf"), 0),
    (-vpx(t)*0.0 + vpy(t)*1.0)*L2(t): (0, float("inf")),

    ((vpx(t)*0.0 + vpy(t)*1.0)**2 - 0.030153689607045803*(vpx(t)**2 + vpy(t)**2)*1.0)*L3(t): (-float("inf"), 0),
    (-vpx(t)*1.0 + vpy(t)*0.0)*L3(t): (0, float("inf")),

    ((vpx(t)*-1.0 + vpy(t)*0.0)**2 - 0.030153689607045803*(vpx(t)**2 + vpy(t)**2)*1.0)*L4(t): (-float("inf"), 0),
    (-vpx(t)*0.0 + vpy(t)*-1.0)*L4(t): (0, float("inf")),

    ((vpx(t)*0.0 + vpy(t)*-1.0)**2 - 0.030153689607045803*(vpx(t)**2 + vpy(t)**2)*1.0)*L5(t): (-float("inf"), 0),
    (-vpx(t)*-1.0 + vpy(t)*0.0)*L5(t): (0, float("inf")),
}

known_trajectory_map = {
    L2(t): np.concatenate([np.ones(num_nodes/2),np.zeros(num_nodes/2)]),
    L3(t): np.concatenate([np.zeros(num_nodes/2),np.zeros(num_nodes/2)]),
    L4(t): np.concatenate([np.zeros(num_nodes/2),np.zeros(num_nodes/2)]),
    L5(t): np.concatenate([np.zeros(num_nodes/2),np.ones(num_nodes/2)]),
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
