# %%
import matplotlib.pyplot as plt
import numpy as np
from fipy import *
from fipy.tools import numerix
from fipy.variables.faceGradVariable import _FaceGradVariable

# 打印完整数组
np.set_printoptions(threshold=np.inf)

# %%

filename = "cylinder"
mesh = Gmsh2D("{}.msh2".format(filename))

# %%

viscosity = 1
U = 10.

pBase=0.1

rho = 100.

p_R = 0.4
v_R = 0.2

# %%

Vx = CellVariable(mesh=mesh, name="x velocity")
Vy = CellVariable(mesh=mesh, name="y velocity")

Vf = FaceVariable(mesh=mesh, rank=1, value=(U, 0))

p = CellVariable(mesh=mesh, name="pressure")
pc = CellVariable(mesh=mesh)


apx = CellVariable(mesh=mesh, value=1.)
apy = CellVariable(mesh=mesh, value=1.)

volume = CellVariable(
    mesh=mesh, value=mesh.cellVolumes, name='Volume'
)
contrvolume = volume.arithmeticFaceValue

X, Y = mesh.faceCenters

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
pc.constrain(pBase, outletFace)

Vx.constrain(0., cylinderFace)
Vy.constrain(0., cylinderFace)
p.faceGrad.constrain(0., cylinderFace)
pc.faceGrad.constrain(0., cylinderFace)

Vx.constrain(U, top_bottomFace)
Vy.faceGrad.constrain(0., top_bottomFace)
p.constrain(pBase, top_bottomFace)
pc.constrain(pBase, top_bottomFace)

# %%

Vx_M_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vx) \
    * rho + p.grad[0] == 0
Vy_M_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vy) \
    * rho + p.grad[1] == 0

coeff = (
    1. / (
        apx.arithmeticFaceValue
        * mesh._faceAreas
        * mesh._cellDistances
    ),
    1. / (
        apy.arithmeticFaceValue
        * mesh._faceAreas
        * mesh._cellDistances
    )
)

# coeff = (
#     1. / apx.arithmeticFaceValue,
#     1. / apy.arithmeticFaceValue
# )
pc_Eq = \
    DiffusionTerm(coeff=coeff, var=pc) == 0


# %%
def sweep():

    Vx_M_Eq.cacheMatrix()
    Vy_M_Eq.cacheMatrix()

    xres = Vx_M_Eq.sweep(
        var=Vx,
        underRelaxation=v_R)
    xmat = Vx_M_Eq.matrix

    yres = Vy_M_Eq.sweep(
        var=Vy,
        underRelaxation=v_R)
    ymat = Vy_M_Eq.matrix

    apx[:] = -numerix.asarray(xmat.takeDiagonal())
    apy[:] = -numerix.asarray(ymat.takeDiagonal())

    pc_Eq.cacheRHSvector()
    pres = pc_Eq.sweep(var=pc)
    rhs = pc_Eq.RHSvector

    p.setValue(p + p_R * pc)

    # Vx.setValue(Vx - pc.grad[0] /
    #                    apx * mesh.cellVolumes)
    # Vy.setValue(Vy - pc.grad[1] /
    #                    apy * mesh.cellVolumes)

    Vx_M_Eq.solve(var=Vx)
    Vy_M_Eq.solve(var=Vy)

    Vf.setValue(
        (Vx.arithmeticFaceValue, Vy.arithmeticFaceValue))

    return xres, yres, pres

# %%


Vx_M_Eq.cacheMatrix()
Vx_M_Eq.cacheRHSvector()

xres = Vx_M_Eq.sweep(
    var=Vx,
    underRelaxation=v_R)
xmat = Vx_M_Eq.matrix
xrhs = Vx_M_Eq.RHSvector

Vy_M_Eq.cacheMatrix()
Vy_M_Eq.cacheRHSvector()

yres = Vy_M_Eq.sweep(
    var=Vy,
    underRelaxation=v_R)
ymat = Vy_M_Eq.matrix
yrhs = Vy_M_Eq.RHSvector

apx[:] = -numerix.asarray(xmat.takeDiagonal())
apy[:] = -numerix.asarray(ymat.takeDiagonal())

pc_Eq.cacheMatrix()
pc_Eq.cacheRHSvector()

pres = pc_Eq.sweep(var=pc)
pmat = pc_Eq.matrix
prhs = pc_Eq.RHSvector

p.setValue(p + p_R * pc)

# Vx.setValue(Vx - pc.grad[0] /
#                    apx * mesh.cellVolumes)
# Vy.setValue(Vy - pc.grad[1] /
#                    apy * mesh.cellVolumes)

Vx_M_Eq.solve(var=Vx)
Vy_M_Eq.solve(var=Vy)

Vf.setValue(
    (Vx.arithmeticFaceValue, Vy.arithmeticFaceValue))


# %%

xrhs/(p.grad[0].value*volume.value)

# %%

res_limit = 1.

# while (1):
for i in range(100):
    xres, yres, pres = sweep()

    info = """ 
    xres= %(xres)g
    yres= %(yres)g
    pres= %(pres)g
     """
    print(info % locals())

    if (max([xres, yres, pres]) < res_limit):
        break


# %%

Viewer(p)
Viewer(Vx)
Viewer(Vy)
