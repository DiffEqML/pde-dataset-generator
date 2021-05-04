import time
import numpy as np
from fenics import *

t0 = 0.0 # Start time
tn = 2.0 # End time
steps = 256 # Number of time step
dt = tn/steps # Time step size

# Create mesh
nx = 256
ny = 256
x0 = -2
xn = 2
y0 = -2
yn = 2

mesh = RectangleMesh(Point(x0, y0), Point(xn, yn), nx, ny)
function_space = FunctionSpace(mesh, 'P', 1)

# Boundary Condition
def boundary(x, on_boundary):
    return on_boundary

boundary_condition = DirichletBC(function_space, Constant(0), boundary)

# Initial Condition
u0 = Expression('exp(-a*pow(x[0],2) - a*pow(x[1],2))', degree=2, a=5)
un = interpolate(u0, function_space)

# Define Problem
u = TrialFunction(function_space)
v = TestFunction(function_space)
f = Constant(0)
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
    solution.append(save_state(u, x0, xn, nx, y0, yn, ny))
    un.assign(u)

save_solution = np.stack(solution, axis=0)
print(save_solution.shape)
np.save('../../data/gaussian.npy', save_solution)
