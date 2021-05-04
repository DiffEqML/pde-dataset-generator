from math import degrees
import time
import numpy as np
from fenics import *

t0 = 0.0 # Start time
tn = 2.0 # End time
steps = 30 # Number of time step
dt = tn/steps # Time step size
alpha = 3 # Function parameter
beta = 1.2 # Function parameter

# Create mesh
nx = 8
ny = 8
x0 = 0
xn = 1
y0 = 0
yn = 1

mesh = UnitSquareMesh(nx, ny)
function_space = FunctionSpace(mesh, 'P', 1)

# Boundary Condition
ud = Expression('1 + x[0] * x[0] + alpha * x[1] * x[1] + beta * t', degree=2, alpha=alpha, beta=beta, t=0)
def boundary(x, on_boundary):
    return on_boundary

boundary_condition = DirichletBC(function_space, ud, boundary)

# Initial Condition
un = interpolate(ud, function_space)

# Define Problem
u = TrialFunction(function_space)
v = TestFunction(function_space)
f = Constant(beta - 2 - 2 * alpha)
F = u * v * dx + dt * dot(grad(u), grad(v)) * dx - (un + dt * f) * v * dx
a = lhs(F)
l = rhs(F)

# Define Describe Function
def save_state(f, x0, xn, nx, y0, yn, ny):
    state = np.empty((nx, ny))
    for i in range(nx):
        for j in range(ny):
            state[i, j] = f(x0 + i*(xn-x0)/nx, y0 + j*(yn-y0)/ny)
    return state

# Solve
u = Function(function_space)
t = t0
solution = []
for i in range(steps):
    t += dt
    solve(a == l, u, boundary_condition)
    solution.append(save_state(u, x0, xn, 30, y0, yn, 30))
    un.assign(u)

save_solution = np.stack(solution, axis=0)
print(save_solution.shape)
np.save('../../data/heat.npy', save_solution)
