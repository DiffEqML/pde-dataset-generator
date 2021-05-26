import torch
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_graph(graph, 
               mesh_color='black', 
               mesh_alpha=0.5, 
               mesh_linewidth=1, 
               contour_level=30,
               *args):
    ''' Draw DGLGraph with node value

    Generate two figures, one will show the mesh, another
    one will show the value with node coordinates matched.

    Args:
        graph: <dgl.graph> should have node features: 'x',
            'y', and 'value' and edges. 
        mesh_color: <str> pyplot color option
        mesh_alpha: <float> pyplot alpha option
        mesh_linewidth: <int> pyplot linewidht option
        contour_level: <int> number of classified level in
            contour graph
        args: <float> 2 the boundary of color bar, the first
            value is the lowest value, and the second value
            is the hightest value
    '''
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    x = graph.ndata['x']
    y = graph.ndata['y']
    value = graph.ndata['value']
    plt.scatter(x, y, value)
    src, dst = graph.edges()
    for i in range(len(dst)):
        nodes_x = [x[src[i]], x[dst[i]]]
        nodes_y = [y[src[i]], y[dst[i]]]
        plt.plot(nodes_x, 
                 nodes_y, 
                 color=mesh_color, 
                 alpha=mesh_alpha, 
                 linewidth=mesh_linewidth)

    ax = plt.subplot(1, 2, 2)
    if args:
        norm = matplotlib.colors.Normalize(vmin=args[0], vmax=args[1])
        fig = plt.tricontourf(x, y, value, levels=contour_level, norm=norm)
    else:
        fig = plt.tricontourf(x, y, value, levels=contour_level)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(fig, cax=cax)
