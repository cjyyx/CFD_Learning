# %%
from fipy import *
from fipy.tools import numerix


# %%
cellSize = 0.05
radius = 1.

mesh = Gmsh2D('''
              cellSize = %(cellSize)g;
              radius = %(radius)g;
              Point(1) = {0, 0, 0, cellSize};
              Point(2) = {-radius, 0, 0, cellSize};
              Point(3) = {0, radius, 0, cellSize};
              Point(4) = {radius, 0, 0, cellSize};
              Point(5) = {0, -radius, 0, cellSize};
              Circle(6) = {2, 1, 3};
              Circle(7) = {3, 1, 4};
              Circle(8) = {4, 1, 5};
              Circle(9) = {5, 1, 2};
              Line Loop(10) = {6, 7, 8, 9};
              Plane Surface(11) = {10};
              ''' % locals()) 

phi = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = 0.) 

viewer = Viewer(vars=phi)


D = 1.
eq =  DiffusionTerm(coeff=D)

X, Y = mesh.faceCenters 
phi.constrain(10, mesh.exteriorFaces) 

eq.solve(var=phi)

# %%

viewer.plot() 
