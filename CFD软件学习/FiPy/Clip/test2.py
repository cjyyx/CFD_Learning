from fipy import *
from fipy.tools import numerix

import numpy as np

geo = '''
// A mesh consisting of a square

// define the corners of the square

Point(1) = {1, 1, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {0, 0, 0, 1};
Point(4) = {1, 0, 0, 1};

// define the square

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// define the boundary

Line Loop(1) = {1, 2, 3, 4};

// define the domain

Plane Surface(1) = {1};
'''

error = []
bkg = None
from builtins import range
for refine in range(4):
    square = Gmsh2D(geo, background=bkg) 
    x, y = square.cellCenters 
    bkg = CellVariable(mesh=square, value=abs(x / 4) + 0.01) 
    error.append(((2 * numerix.sqrt(square.cellVolumes) / bkg - 1)**2).cellVolumeAverage)

phi = CellVariable(name = "solution variable",
                   mesh = square,
                   value = 0.) 
viewer = Viewer(vars=phi, datamin=-1, datamax=1.)
viewer.plotMesh()

input()