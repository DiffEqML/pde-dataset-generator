from fenics import *
from dolfin import *
from mshr import *


def rectangle_static(x0, 
                     xn, 
                     y0, 
                     yn, 
                     ud_top, 
                     ud_bottom, 
                     ud_left, 
                     ud_right,
                     cell_size,
                     tol):
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


def rectangle_dynamic(x0, 
                      xn, 
                      y0, 
                      yn, 
                      ud_top, 
                      ud_bottom, 
                      ud_left, 
                      ud_right,
                      cell_size,
                      t, 
                      tol):
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
    
    ud_top = Expression(ud_top, degree=2, t=t)
    ud_bottom = Expression(ud_bottom, degree=2, t=t)
    ud_left = Expression(ud_left, degree=2, t=t)
    ud_right = Expression(ud_right, degree=2, t=t)
    bc = []
    bc.append(DirichletBC(function_space, ud_top, boundaries, 1))
    bc.append(DirichletBC(function_space, ud_bottom, boundaries, 2))
    bc.append(DirichletBC(function_space, ud_left, boundaries, 3))
    bc.append(DirichletBC(function_space, ud_right, boundaries, 4))

    return function_space, bc


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