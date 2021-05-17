from fenics import *

import matplotlib.pyplot as plt
import numpy as np
from mshr import *

domain = Rectangle(Point(-2, -2), Point(2, 2))
mesh = generate_mesh(domain, 16)
function_space = FunctionSpace(mesh, 'P', 1)