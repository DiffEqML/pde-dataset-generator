

from dolfin.function.functionspace import FunctionSpace
from dolfin.mesh.meshfunction import MeshFunction


def gaussian(a:float, mesh, boundary_list, constrain_list):
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)
    for i in range(len(boundary_list)):
        boundary_list[i].mark(boundaries, i + 1)

    function_space = FunctionSpace(mesh, 'P', 1)
    bc = []
    for i in range(len(constrain_list)):
        bc.append()