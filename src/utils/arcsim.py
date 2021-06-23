import torch
import numpy as np
import pywavefront # for loading .obj files
import logging
logging.getLogger().setLevel(logging.CRITICAL) # disable long output due to error ('ms' field not compatible)
default_keys=['x', 'y', 'z']
default_dtype=torch.float32

def arcsim_to_graph(data, **kwargs):
    if data is None:
        raise RuntimeError('Please provide .obj file to create the graph.')
    return obj_to_graph(data, **kwargs)

def obj_to_graph(file, **kwargs):
    '''Transforms pywavefront object into DGL graph
    Compatible with ArcSim .obj files
    '''
    keys = kwargs.get('keys') if 'keys' in kwargs else default_keys
    dtype = kwargs.get('dtype') if 'dtype' in kwargs else default_dtype
    obj = pywavefront.Wavefront(file, collect_faces=True) # load ArcSim .obj file
    vert = np.array(obj.vertices)
    mesh = np.array(obj.mesh_list[0].faces)
    # We add the sources and destination through contiguous nodes
    # We should also repeat by switching indexes to have bidirectional graph
    # BEWARE: this will have lots of duplicates, don't think DGL can deal with them natively
    # extend method similar to append
    src = []
    dst = []
    for i, j in zip([0,1,2], [1,2,0]): 
        src.extend(mesh[:,i].tolist()); dst.extend(mesh[:,j].tolist())
        src.extend(mesh[:,j].tolist()); dst.extend(mesh[:,i].tolist())
    nodes = vert.shape[0] - 1 # this is the number of nodes, so we have to subtract 1 for DGL graph creation
    edges = [src, dst]
    values = []
    for k, i in zip(keys, range(len(keys))):
        # Spatial coordinates
        if k == 'x': values.append(torch.tensor(vert[:,i], dtype=dtype))
        elif k == 'y': values.append(torch.tensor(vert[:,i], dtype=dtype))
        elif k == 'z': values.append(torch.tensor(vert[:,i], dtype=dtype))
        else:
            # Field values
            values.append(torch.tensor(vert[:,i], dtype=dtype)) # get data corresponding to the key
    return [nodes, edges, keys, values]
