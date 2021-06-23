import torch
import numpy as np
from fenics import *
from dolfin import *
from mshr import *

default_keys = ['x', 'y', 'value']
# Graph Element Networks default keys
gen_keys = ['x', 'y', 'coords', 'value', 'feat',
            'init_feat', 'is_bdd', 'type', 'type_onehot',
            'dist', 'feat', 'init_feat'] 
default_dtype=torch.float32


def fenics_to_graph(mesh, function, **kwargs):
    '''
    Generate graph values from Fenics data
    '''
    if mesh is None or function is None:
        raise RuntimeError('Please provide Fenics mesh and function to create the graph.')
    # Kwargs parsing
    keys = kwargs.get('keys') if 'keys' in kwargs else default_keys
    dtype = kwargs.get('dtype') if 'dtype' in kwargs else default_dtype
    use_gen = kwargs.get('use_gen') if 'use_gen' in kwargs else False
    # Graph creation
    nodes = mesh.num_vertices()
    src, dst = get_edges(mesh)
    edges = [src, dst]
    if use_gen is True:
        keys = gen_keys
        values = _get_gen_values(function, mesh, src, dst)
    else:
        values = []
        for k, i in zip(keys, range(len(keys))):
            # Spatial coordinates
            if k == 'x': values.append(torch.tensor(mesh.coordinates()[:,0], dtype=dtype))
            elif k == 'y': values.append(torch.tensor(mesh.coordinates()[:,1], dtype=dtype))
            elif k == 'z': values.append(torch.tensor(mesh.coordinates()[:,2], dtype=dtype))
            else:
                # Field values
                values.append(torch.tensor(get_values(function, mesh), dtype=dtype)) # get data corresponding to the key
    return [nodes, edges, keys, values]



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


def _get_gen_values(u, mesh, src, dst):
    '''Get the values for the GEN model (Graph Element Networks)
    The order is according to the default GEN keys'''
                              
    coords = mesh.coordinates()
    value = np.array(get_values(u, mesh)).reshape(-1, 1)
    feature = np.hstack([coords, value])
    zero_array = np.zeros((mesh.num_vertices()))
    one_array = np.ones((mesh.num_vertices()))
    onehot_array = np.hstack([zero_array.reshape(-1, 1), one_array.reshape(-1, 1)])

    dist = []
    for i in range(len(src)):
        dist.append(np.linalg.norm(coords[dst[i]] - coords[src[i]]))
    dist = np.array(dist)
    values = []
    values.append(torch.tensor(coords[:,0], dtype=torch.float32))
    values.append(torch.tensor(coords[:,1], dtype=torch.float32))
    values.append(torch.tensor(coords, dtype=torch.float32))
    # WARNNING: value shape here is different from original one
    values.append(torch.tensor(value, dtype=torch.float32))
    values.append(torch.tensor(feature, dtype=torch.float32))
    values.append(torch.tensor(feature, dtype=torch.float32))
    values.append(torch.tensor(zero_array))
    values.append(torch.tensor(one_array))
    values.append(torch.tensor(onehot_array))
    values.append(torch.tensor(dist.reshape(-1, 1), dtype=torch.float32))
    values.append(torch.tensor(dist.reshape(-1, 1), dtype=torch.float32))
    values.append(torch.tensor(dist.reshape(-1, 1), dtype=torch.float32))
    return values
