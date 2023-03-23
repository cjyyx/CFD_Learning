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
L = 1.0
N = 50
dL = L / N
viscosity = 1
U = 10.

pressureRelaxation = 0.9
velocityRelaxation = 0.8

# %%

xVelocity = CellVariable(mesh=mesh, name='X velocity')
yVelocity = CellVariable(mesh=mesh, name='Y velocity')
pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)

velocity = FaceVariable(mesh=mesh, rank=1)

ap = CellVariable(mesh=mesh, value=1.)

volume = CellVariable(
    mesh=mesh, value=mesh.cellVolumes, name='Volume'
)
contrvolume = volume.arithmeticFaceValue

X, Y = mesh.faceCenters

xVelocityEq = \
    DiffusionTerm(coeff=viscosity, var=xVelocity) \
    - pressure.grad.dot([1., 0.])
yVelocityEq = \
    DiffusionTerm(coeff=viscosity, var=yVelocity) \
    - pressure.grad.dot([0., 1.])

coeff = \
    1. / (
        ap.arithmeticFaceValue
        * mesh._faceAreas
        * mesh._cellDistances
    )
pressureCorrectionEq = \
    DiffusionTerm(coeff=coeff) \
    - velocity.divergence

inletFace = mesh.physicalFaces["inlet"]
outletFace = mesh.physicalFaces["outlet"]
cylinderFace = mesh.physicalFaces["cylinder"]
top_bottomFace = mesh.physicalFaces["top"] | mesh.physicalFaces["bottom"]

# xVelocity.constrain(U, inletFace)
# yVelocity.constrain(0., inletFace)
# pressure.grad.constrain(0., inletFace)
# pressureCorrection.grad.constrain(0., inletFace)

# xVelocity.grad.constrain(0., outletFace)
# yVelocity.grad.constrain(0., outletFace)
# pressure.constrain(0., outletFace)
# pressureCorrection.constrain(0., outletFace)

# xVelocity.constrain(0., cylinderFace)
# yVelocity.constrain(0., cylinderFace)
# pressure.grad.constrain(0., cylinderFace)
# pressureCorrection.grad.constrain(0., cylinderFace)

# xVelocity.constrain(U, top_bottomFace)
# # xVelocity.faceGrad.constrain(0., top_bottomFace)
# yVelocity.grad.constrain(0., top_bottomFace)
# pressure.constrain(0., top_bottomFace)
# pressureCorrection.constrain(0., top_bottomFace)


xVelocity.constrain(U, inletFace)
yVelocity.constrain(0., inletFace)
pressure.faceGrad.constrain(0., inletFace)
pressureCorrection.faceGrad.constrain(0., inletFace)

xVelocity.faceGrad.constrain(0., outletFace)
yVelocity.faceGrad.constrain(0., outletFace)
pressure.constrain(0., outletFace)
pressureCorrection.constrain(0., outletFace)

xVelocity.constrain(0., cylinderFace)
yVelocity.constrain(0., cylinderFace)
pressure.faceGrad.constrain(0., cylinderFace)
pressureCorrection.faceGrad.constrain(0., cylinderFace)

xVelocity.constrain(U, top_bottomFace)
# xVelocity.faceGrad.constrain(0., top_bottomFace)
yVelocity.faceGrad.constrain(0., top_bottomFace)
pressure.constrain(0., top_bottomFace)
pressureCorrection.constrain(0., top_bottomFace)


# %%


def sweep():

    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(
        var=xVelocity,
        underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix

    yres = yVelocityEq.sweep(
        var=yVelocity,
        underRelaxation=velocityRelaxation)

    ap[:] = -numerix.asarray(xmat.takeDiagonal())

    presgrad = pressure.grad
    facepresgrad = _FaceGradVariable(pressure)

    velocity[0] = xVelocity.arithmeticFaceValue \
        + contrvolume / ap.arithmeticFaceValue * \
        (presgrad[0].arithmeticFaceValue-facepresgrad[0])
    velocity[1] = yVelocity.arithmeticFaceValue \
        + contrvolume / ap.arithmeticFaceValue * \
        (presgrad[1].arithmeticFaceValue-facepresgrad[1])

    pressureCorrectionEq.cacheRHSvector()
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector

    pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] /
                       ap * mesh.cellVolumes)
    yVelocity.setValue(yVelocity - pressureCorrection.grad[1] /
                       ap * mesh.cellVolumes)

    return xres, yres, pres


# %%

res_limit = 0.5

while(1):
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

Viewer(pressure)
Viewer(xVelocity)
Viewer(yVelocity)
Viewer(velocity)
