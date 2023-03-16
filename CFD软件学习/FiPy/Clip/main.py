# %%
import numpy as np
from fipy import *
from fipy.tools import numerix
from tqdm import tqdm

# 打印完整数组
np.set_printoptions(threshold=np.inf)

# %%

D = 1.

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

eq = TransientTerm() == DiffusionTerm(coeff=D)

# %%

timeStepDuration = 1

for step in tqdm(range(10)):
    eq.solve(var=T,dt=timeStepDuration)


# %%

viewer.plot("result.png")

