# 📦 PDEs Dataset Generator

A tool for generating datasets from several PDE processes based on Fenics and ARCSim.

## 💡 Introduction

Research on PDEs needs ground truth datasets. Several tools can meet this requirement, FEniCS is one of the most famous tools. We here are going to provide some scripts to help us to more easily create such datasets and save them in DGL graph format. 

Features:

- Compressed scripts, contain full generating process, including create domain, generate mesh, create boundary constrains, solve function, transfer mesh and result to dgl graph, and then save them. 
- Results will be trainsferred to dgl graph, which is more convenient to use in graph model. We also going to provide numpy tools. 
- In detail tutorial. We provide notebooks to show how our scripts work and knowledge about using FEniCS. 
- Provide useful mini tools, including mesh to dgl transferring, dgl grpah plotting, etc. 

## 🔧 Environment 

```Python
Python Version: 3.6 or later
Python Packages: jupyterlab, fenics, dgl, numpy, torch, matplotlib
```

## 📁 Structure

```
.
├── fig/
├── notebook/
│   └── *
├── src/
│   ├── utils/
│   └── *
└── README.md
```

- `fig`: example figures
- `notebook`: tutorials in jupyter notebook format
- `src`: all source code will be here, including scripts, tools
  - `utils`: mini tools will be here, including dgl transferring, graph plot

Here are guides for Fenics and ARCSim:
- [Guide with Fenics](#guide-with-fenics):
- [Guide with ARCSim](#guide-with-arcsim):

## 💾 Guide with Fenics

**Step 1**. [Download](https://github.com/cbhua/tool-pdeset-generator/archive/refs/heads/main.zip) or [Clone](https://github.com/cbhua/tool-pdeset-generator.git) this repository.  

**Step 2**. Based on your requirement refer to notebooks, there would be tutorials and examples. You can find all methods provided in below list. 

**Step 3**. Modify parameters to generate your datasets. 

Provided methods:

- Poisson process
  - Customize domain & Single boundary control
  - Square domain & Separate boundary control
- Gaussian process
  - Customize domain & Single boundary control (support time dynamic control)
  - Square domain & Separate boundary control (support time dynamic control)
  - Squares in square domain & Separate boundary control (support time dynamic control)
  - Circles in circle domain & Single boundary control (support time dynamic control)

Support methods will keep updating. For more detail, you can refer to the [project manager](https://github.com/cbhua/tool-pdeset-generator/projects/1). 

## 📊 Examples

### Possion process, square domain, single boundary control

<img src="fig/possion_square.png" alt="possion process, square domain, single boundary control" style="zoom:60%;" />

### Possion process, L shape domain, single boundary control

<img src="fig/possion_l.png" alt="possion process, L shape domain, single boundary control" style="zoom:60%;" />

### Possion process, circle shape domain, single boundary control

<img src="fig/possion_circle.png" alt="possion process, circle shape domain, single boundary control" style="zoom:60%;" />

### Gaussian process, rectangle shape domain, single boundary control

<img src="fig/gaussian_square.gif" alt="gaussian process, rectangle shape domain, single boundary control" style="zoom:80%;" />

### Gaussian process, rectangle shape domain, multi & dynamic boundary control

<img src="fig/gaussian_square_dynamic.gif" alt="gaussian process, rectangle shape domain, multi & dynamic boundary control" style="zoom:80%;" />

## 💾 Guide with Fenics
**Step 1**. [Download](https://github.com/cbhua/tool-pdeset-generator/archive/refs/heads/main.zip) or [Clone](https://github.com/cbhua/tool-pdeset-generator.git) this repository. 

**Step 2**. ArcSim installation

You may find the repository with fixes [here](https://github.com/kaist-silab/arcsim) with further instructions.
To install it, run the following:

`git clone https://github.com/kaist-silab/arcsim.git && cd arcsim/`

`sudo chmod +x install.sh && sudo ./install.sh`

At this point, you should be ready to go.

**Step 3**. ArcSim simulation and `.obj` file saving
Let's consider the flag example. In the ArcSim folder, make a new directory called data. Then run:

`bin/arcsim simulate conf/flag.json data/`

(you may also run `simulateoffline` if you cannot visualize on your computer

When the simulation ends (we may do that with `Esc` as well) copy the `conf/flag.json` into the folder where we saved the simulation, in our case`data/` and run:

`bin/arcsim generate data/`

This will generate `.obj` files that we can load into Python with `pywavefront` and the `obj_to_dgl` method we provide.


## 📊 Examples
### Flag simulation with adaptive remeshing

<img src="fig/flag.gif" alt="Flag simulation with adaptive remeshing" style="zoom:80%;" />


## 📜 References

1. FEniCS project: https://fenicsproject.org/
2. ARCSim project: http://graphics.berkeley.edu/resources/ARCSim/
3. DGL project: https://www.dgl.ai/
