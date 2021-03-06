import os
from dgl.convert import graph
import torch 
import numpy as np 
import dgl
from fenics import *
from dolfin import *
from mshr import *
from dgl.data.utils import save_graphs
from ufl.mathfunctions import Exp
from src.utils.to_dgl import to_dgl
from src.utils.bc import rectangle


def possion(mesh,  
            ud, 
            f:str='0',
            path:str='data/possion.bin'
            ):
    ''' Create Possion Equation Dataset

    2D Possion queation, static process, custom domain, 
    single dirichlet boundary condition.

    Args:
        mesh: <dolfin.cpp.mesh.Mesh> custom mesh generated by fenics mesh tools
        f: <str> right part of laplace function, in cpp argument format
        ud: <str> boundary condition function, in cpp argument format
        path: <str> path for saving generated dgl graph, in .bin format~
    '''
    function_space = FunctionSpace(mesh, 'P', 1)

    u = TrialFunction(function_space)
    v = TestFunction(function_space)
    f = Expression(f, degree=2)
    ud = Expression(ud, degree=2)
    bc = DirichletBC(function_space, ud, boundary)

    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx
    u = Function(function_space)
    solve(a == L, u, bc)

    graph = to_dgl(function=u, mesh=mesh)
    save_graphs(path, graph)


def possion_square(x0:float,
                   xn:float, 
                   y0:float, 
                   yn:float,
                   f:str='0',
                   ud_top:str='0',
                   ud_bottom:str='0',
                   ud_left:str='0',
                   ud_right:str='0',
                   cell_size:float=5., 
                   tol:float=1e-4, 
                   path:str='data/possion_square.bin'
                   ):
    ''' Create Possion Equation Dataset in Square Domain

    2D Possion queation, static process, rectangle domain, 
    can custom each boundary's condition. 

    Args:
        x0: <float> left boundary for x
        xn: <float> right boundary for x
        y0: <float> left boundary for y
        yn: <float> right boundary for y
        f: <str> right part of laplace function, in cpp argument format
        ud_top: <str> boundary condition on the top of rectangle
        ud_bottom: <str> boundary condition on the bottom of rectangle
        ud_left: <str> boundary condition on the left of rectangle
        ud_right: <str> boundary condition on the right of rectangle
        cell_siez: <float> cell size for created mesh
        tol: <float> boundary bias, e.g. (x-tol, x+tol) is a boundary on x
        path: <str> path for saving generated dgl graph, in .bin format
    '''
    mesh, function_space, bc = rectangle(x0,
                                         xn, 
                                         y0, 
                                         yn, 
                                         ud_top, 
                                         ud_bottom, 
                                         ud_left,
                                         ud_right, 
                                         cell_size, 
                                         tol
                                         )

    u = TrialFunction(function_space)
    v = TestFunction(function_space)
    f = Expression(f, degree=2)

    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx
    u = Function(function_space)
    solve(a == L, u, bc)

    graph = to_dgl(function=u, mesh=mesh)
    save_graphs(path, graph)


def boundary(x, on_boundary):
    return on_boundary
