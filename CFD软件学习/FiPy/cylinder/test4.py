# %%
import matplotlib.pyplot as plt
import numpy as np
from fipy import *
from fipy.tools import numerix
from fipy.variables.faceGradVariable import _FaceGradVariable
from tqdm import tqdm

# 打印完整数组
np.set_printoptions(threshold=np.inf)

# %%

filename = "cylinder"
mesh = Gmsh2D("{}.msh2".format(filename))

# %%
L = 1.0
N = 50
dL = L / N
viscosity = 1
U = 10.

pBase = 0.

rho = 100.

Rp = 0.4
Rv = 0.2

# %%

volume = CellVariable(mesh=mesh, value=mesh.cellVolumes)
contrvolume = volume.arithmeticFaceValue
Xc, Yc = mesh.cellCenters
Xf, Yf = mesh.faceCenters
Vc = mesh.cellVolumes

Vx = CellVariable(mesh=mesh, name="x velocity", value=U)
Vy = CellVariable(mesh=mesh, name="y velocity", value=0.)

Vf = FaceVariable(mesh=mesh, rank=1)
Vf.setValue((Vx.faceValue, Vy.faceValue))

p = CellVariable(mesh=mesh, name="pressure", value=-Xc/100.)
pc = CellVariable(mesh=mesh, value=0.)

Temp = CellVariable(mesh=mesh, rank=1)

apx = CellVariable(mesh=mesh, value=1.)
apy = CellVariable(mesh=mesh, value=1.)


inletFace = mesh.physicalFaces["inlet"]
outletFace = mesh.physicalFaces["outlet"]
cylinderFace = mesh.physicalFaces["cylinder"]
top_bottomFace = mesh.physicalFaces["top"] | mesh.physicalFaces["bottom"]

Vx.constrain(U, inletFace)
Vy.constrain(0., inletFace)
p.faceGrad.constrain(0., inletFace)
pc.faceGrad.constrain(0., inletFace)

Vx.faceGrad.constrain(0., outletFace)
Vy.faceGrad.constrain(0., outletFace)
p.constrain(pBase, outletFace)
pc.constrain(0., outletFace)

Vx.constrain(0., cylinderFace)
Vy.constrain(0., cylinderFace)
p.faceGrad.constrain(0., cylinderFace)
pc.faceGrad.constrain(0., cylinderFace)

Vx.faceGrad.constrain(0., top_bottomFace)
Vy.faceGrad.constrain(0., top_bottomFace)
p.constrain(pBase, top_bottomFace)
pc.constrain(0., top_bottomFace)

# %%

Vx_M_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vx) \
    * rho + ImplicitSourceTerm(coeff=1.0, var=p.grad[0]) == 0
Vy_M_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vy) \
    * rho + ImplicitSourceTerm(coeff=1.0, var=p.grad[1]) == 0

coeff = (
    1. / (
        apx.faceValue
        * mesh._faceAreas
        * mesh._cellDistances
    ),
    1. / (
        apy.faceValue
        * mesh._faceAreas
        * mesh._cellDistances
    )
)

# coeff = (
#     1. / apx.faceValue,
#     1. / apy.faceValue
# )

p_Eq = \
    DiffusionTerm(coeff=coeff, var=p) \
    == Temp.faceValue.divergence


# %%
def sweep():

    Vx_M_Eq.cacheMatrix()
    xres = Vx_M_Eq.sweep(var=Vx, underRelaxation=Rv)
    xmat = Vx_M_Eq.matrix
    apx[:] = numerix.asarray(xmat.takeDiagonal())

    Vy_M_Eq.cacheMatrix()
    yres = Vy_M_Eq.sweep(var=Vy, underRelaxation=Rv)
    ymat = Vy_M_Eq.matrix
    apy[:] = numerix.asarray(ymat.takeDiagonal())

    Temp.setValue((Vx-xres/apx, Vy-yres/apy))

    pres = p_Eq.sweep(var=p, underRelaxation=Rp)

    Vx.setValue(Temp[0] - (p.grad[0] / apx))
    Vy.setValue(Temp[1] - (p.grad[1] / apy))
    Vf.setValue((Vx.faceValue, Vy.faceValue))

    return xres, yres, pres

# %%


res_limit = 1.

# info = "xres= %(xres)g, yres= %(yres)g, pres= %(pres)g"
# info = """
# xres= %(xres)g
# yres= %(yres)g
# pres= %(pres)g
#  """

step_num = 10

pbar = tqdm(range(step_num))
# while (1):
for i in pbar:
    xres, yres, pres = sweep()

    # print(info % locals())

    pbar.set_postfix({
        "xres": f'{xres:.2e}',
        "yres": f'{yres:.2e}',
        "pres": f'{pres:.2e}'
    })

    if (max([xres, yres, pres]) < res_limit):
        break


# %%

Viewer(p)
Viewer(Vx, datamin=-U/5, datamax=U*2)
Viewer(Vy)
