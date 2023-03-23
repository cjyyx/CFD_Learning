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

U = 10.
rho = 100.

Rp = 0.8
Rv = 0.8

# %%
Xc, Yc = mesh.cellCenters
Xf, Yf = mesh.faceCenters
Vc = mesh.cellVolumes

Vx = CellVariable(mesh=mesh, name="x velocity", value=U)
Vy = CellVariable(mesh=mesh, name="y velocity", value=0.)

Vf = FaceVariable(mesh=mesh, rank=1)
Vf.setValue((Vx.faceValue, Vy.faceValue))

p = CellVariable(mesh=mesh, name="pressure", value=-Xc/100.)
pc = CellVariable(mesh=mesh, value=0.)

apx = CellVariable(mesh=mesh, value=1.)
apy = CellVariable(mesh=mesh, value=1.)

# %%

inletFace = mesh.physicalFaces["inlet"]
outletFace = mesh.physicalFaces["outlet"]
cylinderFace = mesh.physicalFaces["cylinder"]
top_bottomFace = \
    mesh.physicalFaces["top"] | mesh.physicalFaces["bottom"]

Vx.constrain(U, inletFace)
Vy.constrain(0., inletFace)
p.faceGrad.constrain(0., inletFace)
pc.faceGrad.constrain(0., inletFace)

Vx.faceGrad.constrain(0., outletFace)
Vy.faceGrad.constrain(0., outletFace)
p.constrain(0., outletFace)
pc.constrain(0., outletFace)

Vx.constrain(0., cylinderFace)
Vy.constrain(0., cylinderFace)
p.faceGrad.constrain(0., cylinderFace)
pc.faceGrad.constrain(0., cylinderFace)

Vx.faceGrad.constrain(0., top_bottomFace)
Vy.faceGrad.constrain(0., top_bottomFace)
p.constrain(0., top_bottomFace)
pc.constrain(0., top_bottomFace)

# %%

Vx_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vx) * rho \
    + ImplicitSourceTerm(coeff=1.0, var=p.grad[0])
Vy_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vy) * rho \
    + ImplicitSourceTerm(coeff=1.0, var=p.grad[1])

# coeff = (
#     1. / (
#         apx.faceValue
#         * mesh._faceAreas
#         * mesh._cellDistances
#     ),
#     1. / (
#         apy.faceValue
#         * mesh._faceAreas
#         * mesh._cellDistances
#     )
# )
coeff = (
    1. / (
        apx.faceValue
        * mesh._faceAreas
        * mesh._cellDistances
    )
)
pc_Eq = \
    DiffusionTerm(coeff=coeff, var=pc) \
    - Vf.divergence

# %%


def sweep():
    Vx_Eq.cacheMatrix()
    # Vx_Eq.cacheRHSvector()

    xres = Vx_Eq.sweep(var=Vx, underRelaxation=Rv)

    xmat = Vx_Eq.matrix
    # xrhs = Vx_Eq.RHSvector
    apx[:] = numerix.asarray(xmat.takeDiagonal())

    # Vy_Eq.cacheMatrix()
    # Vy_Eq.cacheRHSvector()

    yres = Vy_Eq.sweep(var=Vy, underRelaxation=Rv)

    # ymat = Vy_Eq.matrix
    # yrhs = Vy_Eq.RHSvector
    # apy[:] = numerix.asarray(ymat.takeDiagonal())

    # 注意
    # xrhs == -p.grad[0].value * Vc
    # yrhs == -p.grad[1].value * Vc

    Vf.setValue((Vx.faceValue, Vy.faceValue))

    pcres = pc_Eq.sweep(var=pc)
    p.setValue(p + Rp * pc)

    Vx.setValue(Vx-(Vc*pc.grad[0])/apx)
    Vy.setValue(Vy-(Vc*pc.grad[1])/apx)

    return xres, yres, pcres

# %%

step_num = 200
res_limit = 66

pbar = tqdm(range(step_num))
for i in pbar:
    xres, yres, pcres = sweep()

    pbar.set_postfix({
        "xres": f'{xres:.2e}',
        "yres": f'{yres:.2e}',
        "pcres": f'{pcres:.2e}'
    })

    if (sum([xres, yres, pcres]) < res_limit):
        break

# %%

Viewer(p)
Viewer(Vx, datamin=-U/5, datamax=U*2)
Viewer(Vy)
