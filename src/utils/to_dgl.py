import dgl
import torch
import numpy as np

'''
The SIMULATOR_to_graph methods return the following values:
    - nodes: int, number of nodes of the graph
    - edges: list, format [source, destination]
    - keys: list of strings
    - values: list of tensors, corresponding to the keysfrom .arcsim import arcsim_to_graph
'''
from .arcsim import arcsim_to_graph
from .fenics import fenics_to_graph
from .su2 import su2_to_graph

def to_dgl(data=None,
           mesh=None,
           function=None,
           method=None,
           **kwargs):
    '''Obtain DGL graphs from different methods
    Supported file types:
        - arcsim (.obj)
        - su2 (.su2, .vtu, .vtk)
    If the method is not provided, we check automatically
    Example kwargs: 
        - keys: list of specific field keys to look for
        - gen_mode: use graph element networks
    '''
    # Check for method override
    if method is not None:
        if method == 'arcsim':
            graph_data = arcsim_to_graph(data, **kwargs)
        elif method == 'fenics':
            graph_data = fenics_to_graph(mesh, function, **kwargs)
        elif method == 'su2':
            graph_data = su2_to_graph(data, mesh, **kwargs)
        else:
            raise ValueError("Method not supported. Currently available methods are arcsim, fenics and su2.")
    elif data is not None:
        # Automatic file detection
        import os
        filename, file_extension = os.path.splitext(data)
        if file_extension == '.obj':
            graph_data = arcsim_to_graph(data, **kwargs)
        elif file_extension == '.vtu' or file_extension == '.vtk':
            graph_data = su2_to_graph(data, mesh, **kwargs)
        else:
            raise ValueError('Detected file type is {}. Currently readable types for graph data are .obj, .bin, .vtu and .vtk').format(file_extension)
    elif mesh is not None and function is not None:
        graph_data = fenics_to_graph(mesh, function, **kwargs)
    else: 
        raise ValueError('No data provided!')
    # Graph is in a common for all the methods
    nodes, edges, keys, values = graph_data
    return create_dgl_graph(nodes, edges, keys, values)


def create_dgl_graph(nodes, edges, keys, values):
    '''Minimalistic approach to graph creation in DGL'''
    graph = dgl.DGLGraph()
    graph.add_nodes(nodes)
    graph.add_edges(edges[0], edges[1]) # edges: [sources, destinations]
    for k, i in zip(keys, range(len(keys))):
        graph.ndata[k] = values[i]
    return graph

