import dgl
import torch
import numpy as np
from fenics import *
from dolfin import *
from mshr import *


def to_dgl(u, mesh):
    '''
    TODO: select premeter you want
    '''
    graph = dgl.DGLGraph()
    src, dst = get_edges(mesh)
    graph.add_nodes(mesh.num_vertices())
    graph.add_edges(src, dst)
    graph.ndata['x'] = torch.tensor(mesh.coordinates()[:,0])
    graph.ndata['y'] = torch.tensor(mesh.coordinates()[:,1])
    graph.ndata['value'] = torch.tensor(get_values(u, mesh))
    return graph


def get_edges(mesh):
    ''' Get edges for graph from mesh

    Get two connectivity lists of nodes, edges are directed from
    src to dst, required by DGL in the graph creation.
    
    Example:
        Edge 0->0, 0->1, 1->2 will be returned by [0, 0, 1], [0, 1, 2]

    Args:
        mesh: <dolfin.cpp.mesh.Mesh> dolfin mesh
        
    Returns:
        src: <list> [edge number] source node indexes for edges
        dst: <list> [edge number] distinate node indexes for edges
    '''
    mesh.init(0,1)
    src = []
    dst = []
    for v in vertices(mesh):
        idx = v.index()
        neighbors = [Edge(mesh, i).entities(0) for i in v.entities(1)]
        neighbors = np.array(neighbors).flatten()
        neighbors = neighbors[np.where(neighbors != idx)[0]]
        for n in neighbors:
            src.append(int(idx))
            dst.append(int(n))
    return src, dst


def get_values(u, mesh):
    ''' Get values for each mesh node

    Provid a function based on coordinate for mesh nodes, this can 
    calculate values in each node by index order. Here the function 
    is mostly get by FEniCS solver. 

    Args:
        u: <dolfin.function> <any other function> function based on 
            nodes' coordinates 
        mesh: <dolfin.cpp.mesh.Mesh> dolfin mesh

    Returns:
        values: <list> [node number] calculated result by index order
    '''
    values = []
    for pos in mesh.coordinates():
        values.append(u(Point(pos)))
    return values