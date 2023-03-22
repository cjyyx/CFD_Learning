# %%
import numpy as np
from fipy import *
from fipy.tools import numerix

# 打印完整数组
np.set_printoptions(threshold=np.inf)

# %%

T_hot = 80.
T_cool = 20.

# %%

mesh = Gmsh2D("clip.msh2")

T = CellVariable(name="T",
                 mesh=mesh,
                 value=T_cool)

T.constrain(T_hot, mesh.physicalFaces["hot"])
T.constrain(T_cool, mesh.physicalFaces["cool"])
T.faceGrad.constrain(0., mesh.physicalFaces["edge"])

# %%

viewer = Viewer(
    vars=T,
    datamin=T_cool,
    datamax =T_hot
    )

# %%

eq = DiffusionTerm(var=T) == 0

eq.solve(var=T)


# %%

viewer.plot("result.png")

