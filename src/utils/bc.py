import numpy as np
from fenics import *
from dolfin import *
from mshr import *


def rectangle(x0, 
              xn, 
              y0, 
              yn, 
              ud_top, 
              ud_bottom, 
              ud_left, 
              ud_right,
              cell_size,
              tol,
              t=-1):
    domain = Rectangle(Point(x0, xn), Point(y0, yn))
    mesh = generate_mesh(domain, cell_size)
    function_space = FunctionSpace(mesh, 'P', 1)

    top = YBoundary(yn, tol)
    bottom = YBoundary(y0, tol)
    left = XBoundary(x0, tol)
    right = XBoundary(xn, tol)

    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    top.mark(boundaries, 1)
    bottom.mark(boundaries, 2)
    left.mark(boundaries, 3)
    right.mark(boundaries, 4)
    
    if t >= 0:
        ud_top = Expression(ud_top, degree=2, t=t)
        ud_bottom = Expression(ud_bottom, degree=2, t=t)
        ud_left = Expression(ud_left, degree=2, t=t)
        ud_right = Expression(ud_right, degree=2, t=t)
    else: 
        ud_top = Expression(ud_top, degree=2)
        ud_bottom = Expression(ud_bottom, degree=2)
        ud_left = Expression(ud_left, degree=2)
        ud_right = Expression(ud_right, degree=2)
    bc = []
    bc.append(DirichletBC(function_space, ud_top, boundaries, 1))
    bc.append(DirichletBC(function_space, ud_bottom, boundaries, 2))
    bc.append(DirichletBC(function_space, ud_left, boundaries, 3))
    bc.append(DirichletBC(function_space, ud_right, boundaries, 4))

    return function_space, bc


def circle(x,
           y,
           r,
           ud, 
           cell_size, 
           tol,
           t=-1):
    domain = Circle(Point(x, y), r)
    mesh = generate_mesh(domain, cell_size)
    ud = Expression(ud, degree=2)
    function_space = FunctionSpace(mesh, 'P', 1)
    if t >= 0:
        ud = Expression(ud, degree=2, t=t)
    else:
        bc = DirichletBC(function_space, ud, boundary)
    return function_space, bc


def multi_rectangle(num, 
                    type, 
                    x0, 
                    xn,
                    y0, 
                    yn, 
                    ud_top, 
                    ud_bottom, 
                    ud_left, 
                    ud_right,
                    cell_size,
                    tol,
                    t):
    '''
    '''
    domain = Rectangle(Point(0, 0), Point(0, 0))
    for i in range(num):
        if type[i] == 1:
            domain += Rectangle(Point(x0[i], y0[i]), Point(xn[i], yn[i]))
        else:
            domain -= Rectangle(Point(x0[i], y0[i]), Point(xn[i], yn[i]))
    mesh = generate_mesh(domain, cell_size)
    function_space = FunctionSpace(mesh, 'P', 1)
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    bc = []
    for i in range(num):
        top = YBoundary(yn[i], tol)
        bottom = YBoundary(y0[i], tol)
        left = XBoundary(x0[i], tol)
        right = XBoundary(xn[i], tol)

        top.mark(boundaries, 4 * i + 1)
        bottom.mark(boundaries, 4 * i + 2)
        left.mark(boundaries, 4 * i + 3)
        right.mark(boundaries, 4 * i + 4)

        if t >= 0:
            ud_top = Expression(ud_top[i], degree=2, t=t)
            ud_bottom = Expression(ud_bottom[i], degree=2, t=t)
            ud_left = Expression(ud_left[i], degree=2, t=t)
            ud_right = Expression(ud_right[i], degree=2, t=t)
        else: 
            ud_top = Expression(ud_top[i], degree=2)
            ud_bottom = Expression(ud_bottom[i], degree=2)
            ud_left = Expression(ud_left[i], degree=2)
            ud_right = Expression(ud_right[i], degree=2)
        
        bc.append(DirichletBC(function_space, ud_top, boundaries, 4 * i + 1))
        bc.append(DirichletBC(function_space, ud_bottom, boundaries, 4 * i + 2))
        bc.append(DirichletBC(function_space, ud_left, boundaries, 4 * i + 3))
        bc.append(DirichletBC(function_space, ud_right, boundaries, 4 * i + 4))

    return function_space, bc


def multi_circle(num, 
                 type, 
                 x,
                 y,
                 r,
                 ud, 
                 cell_size, 
                 tol,
                 t):
    domain = Circle(Point(0, 0), 0)
    for i in range(num):
        if type[i] == 1:
            domain += Circle(Point(x[i], y[i]), r[i])
        else:
            domain -= Circle(Point(x[i], y[i]), r[i])
    mesh = generate_mesh(domain, cell_size)
    function_space = FunctionSpace(mesh, 'P', 1)
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    bc = []
    for i in range(num):
        round = CircleBoundary(x[i], y[i], r[i], tol)
        round.mark(boundaries, i)
        if t >= 0:
            ud_temp = Expression(ud[i], degree=2, t=t)
        else: 
            ud_temp = Expression(ud[i], degree=2)
        bc.append(DirichletBC(function_space, ud_temp, boundaries, i))
        
    return function_space, bc


def boundary(x, on_boundary):
    return on_boundary


class CircleBoundary(SubDomain):
    def __init__(self, x, y, r, tol):
        SubDomain.__init__()
        self.x = x
        self.y = y
        self.r = r
        self.tol = tol
    def inside(self, x, on_boundary):
        flag = np.linalg.norm(x - [self.x, self.y])
        return near(flag, self.r, self.tol)


class XBoundary(SubDomain):
    def __init__(self, value, tol):
        SubDomain.__init__()
        self.value = value
        self.tol = tol
    def inside(self, x, on_boundary):
        return near(x[0], self.value, self.tol)


class YBoundary(SubDomain):
    def __init__(self, value, tol):
        SubDomain.__init__()
        self.value = value
        self.tol = tol
    def inside(self, x, on_boundary):
        return near(x[1], self.value, self.tol)  