# References: https://github.com/locuslab/cfd-gcn/blob/master/mesh_utils.py
import torch
import numpy as np
from os import PathLike
from typing import Sequence, Dict, Union, Tuple, List
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

UnionTensor = Union[torch.Tensor, np.ndarray]

SU2_SHAPE_IDS = {
    'line': 3,
    'triangle': 5,
    'quad': 9,
}

default_dtype = torch.float32
keys = ['x', 'y', 'z']

def su2_to_graph(data, mesh, **kwargs):
    if data is None:
        raise RuntimeError('Please provide .vtu or .vtk data file to create the graph.')
    if mesh is None:
        raise RuntimeError('Please provide .su2 mesh file to create the graph.')
    dtype = kwargs.get('dtype') if 'dtype' in kwargs else default_dtype
    nodes, edges, _, _ = get_mesh_graph(mesh)
    nodes = len(nodes)
    # Read the source file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(data)
    reader.Update()
    output = reader.GetOutput()
    out = dsa.WrapDataObject(output)
    if 'keys' in kwargs:
        keys = kwargs.get('keys') 
    else:
        keys = ['x', 'y', 'z'] # default spatial coordinates
        keys.extend(out.PointData.keys()) # read the fields
    values = []
    for k, i in zip(keys, range(len(keys))):
        # Spatial coordinates
        if k == 'x': values.append(torch.tensor(out.Points[:, 0], dtype=dtype))
        elif k == 'y': values.append(torch.tensor(out.Points[:, 1], dtype=dtype))
        elif k == 'z': values.append(torch.tensor(out.Points[:, 2], dtype=dtype))
        else:
        # Field values
            values.append(torch.tensor(out.PointData[k], dtype=dtype)) # get data corresponding to the key
    return [nodes, edges, keys, values]

def get_mesh_graph(mesh_filename: Union[str, PathLike],
                   dtype: np.dtype = np.float32
                   ) -> Tuple[np.ndarray, np.ndarray, List[List[List[int]]], Dict[str, List[List[int]]]]:
    
    def get_rhs(s: str) -> str:
        return s.split('=')[-1]

    marker_dict = {}
    with open(mesh_filename) as f:
        for line in f:
            if line.startswith('NPOIN'):
                num_points = int(get_rhs(line))
                mesh_points = [[float(p) for p in f.readline().split()[:2]]
                               for _ in range(num_points)]
                nodes = np.array(mesh_points, dtype=dtype)

            if line.startswith('NMARK'):
                num_markers = int(get_rhs(line))
                for _ in range(num_markers):
                    line = f.readline()
                    assert line.startswith('MARKER_TAG')
                    marker_tag = get_rhs(line).strip()
                    num_elems = int(get_rhs(f.readline()))
                    marker_elems = [[int(e) for e in f.readline().split()[-2:]]
                                    for _ in range(num_elems)]
                    # marker_dict[marker_tag] = np.array(marker_elems, dtype=np.long).transpose()
                    marker_dict[marker_tag] = marker_elems

            if line.startswith('NELEM'):
                edges = []
                triangles = []
                quads = []
                num_edges = int(get_rhs(line))
                for _ in range(num_edges):
                    elem = [int(p) for p in f.readline().split()]
                    if elem[0] == SU2_SHAPE_IDS['triangle']:
                        n = 3
                        triangles.append(elem[1:1+n])
                    elif elem[0] == SU2_SHAPE_IDS['quad']:
                        n = 4
                        quads.append(elem[1:1+n])
                    else:
                        raise NotImplementedError
                    elem = elem[1:1+n]
                    edges += [[elem[i], elem[(i+1) % n]] for i in range(n)]
                edges = np.array(edges, dtype=np.int32).transpose()
                # triangles = np.array(triangles, dtype=np.long)
                # quads = np.array(quads, dtype=np.long)
                elems = [triangles, quads]

    return nodes, edges, elems, marker_dict