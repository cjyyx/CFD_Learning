# %%
import matplotlib.pyplot as plt
import numpy as np
from fipy import *
from fipy.tools import numerix
from tqdm import tqdm

# %%

filename = "cylinder"
mesh = Gmsh2D("{}.msh2".format(filename))

# %%

mu = 0.1
rho = 1.

U = 10.

sum_res_list = []


# %%
Xc, Yc = mesh.cellCenters
Vc = mesh.cellVolumes
Vcf = CellVariable(mesh=mesh, value=Vc).faceValue

Vx = CellVariable(mesh=mesh, name="x velocity", value=U)
Vy = CellVariable(mesh=mesh, name="y velocity", value=0.)

Vf = FaceVariable(mesh=mesh, rank=1)
Vf.setValue((Vx.faceValue, Vy.faceValue))

p = CellVariable(mesh=mesh, name="pressure", value=0.)
pc = CellVariable(mesh=mesh, value=0.)

apx = CellVariable(mesh=mesh, value=1.)

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
    UpwindConvectionTerm(coeff=Vf, var=Vx) * rho == \
    DiffusionTerm(coeff=mu, var=Vx) - \
    ImplicitSourceTerm(coeff=1.0, var=p.grad[0])
Vy_Eq = \
    UpwindConvectionTerm(coeff=Vf, var=Vy) * rho == \
    DiffusionTerm(coeff=mu, var=Vy) - \
    ImplicitSourceTerm(coeff=1.0, var=p.grad[1])


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
V_limit = 1e2
p_limit = 2e3


def OverflowPrevention():
    Vx[Vx.value > V_limit] = V_limit
    Vx[Vx.value < -V_limit] = -V_limit

    Vy[Vy.value > V_limit] = V_limit
    Vy[Vy.value < -V_limit] = -V_limit

    Vf[Vf.value > V_limit] = V_limit
    Vf[Vf.value < -V_limit] = -V_limit

    p[p.value > p_limit] = p_limit
    p[p.value < -p_limit] = -p_limit


def sweep(Rp, Rv):
    OverflowPrevention()

    Vx_Eq.cacheMatrix()
    xres = Vx_Eq.sweep(var=Vx, underRelaxation=Rv)
    xmat = Vx_Eq.matrix
    apx[:] = numerix.asarray(xmat.takeDiagonal())

    yres = Vy_Eq.sweep(var=Vy, underRelaxation=Rv)

    # Vf.setValue((Vx.faceValue, Vy.faceValue))

    presgrad = p.grad
    facepresgrad = presgrad.faceValue
    Vf[0] = Vx.faceValue + Vcf / apx.faceValue * \
        (presgrad[0].faceValue-facepresgrad[0])
    Vf[1] = Vy.faceValue + Vcf / apx.faceValue * \
        (presgrad[1].faceValue-facepresgrad[1])

    pcres = pc_Eq.sweep(var=pc)

    p.setValue(p + Rp * pc)

    Vx.setValue(Vx-(Vc*pc.grad[0])/apx)
    Vy.setValue(Vy-(Vc*pc.grad[1])/apx)

    # Vf.setValue((Vx.faceValue, Vy.faceValue))

    presgrad = p.grad
    facepresgrad = presgrad.faceValue
    Vf[0] = Vx.faceValue + Vcf / apx.faceValue * \
        (presgrad[0].faceValue-facepresgrad[0])
    Vf[1] = Vy.faceValue + Vcf / apx.faceValue * \
        (presgrad[1].faceValue-facepresgrad[1])

    return xres, yres, pcres

# %%


def is_increasing(arr):
    for i in range(len(arr)-1):
        if arr[i] >= arr[i+1]:
            return False
    return True


def valueRange(val, a, b):
    """ (a,b] """
    if (val > a and val <= b):
        return True
    else:
        return False

# %%


MaxSweep = 120
res_limit = 1e-5
sum_res = 1e10

Rp = 0.8
Rv = 0.5

pbar = tqdm(range(MaxSweep))
for i in pbar:

    if sum_res > 2e2:
        Rp = 0.8
        Rv = 0.5
    elif valueRange(sum_res, 28., 2e2):
        Rp = 0.92
        Rv = 0.8
    elif valueRange(sum_res, 10., 28.):
        Rp = 0.95
        Rv = 0.8
    elif valueRange(sum_res, 1., 10.):
        Rp = 0.98
        Rv = 0.9
    elif valueRange(sum_res, 0, 1.):
        Rp = 0.99
        Rv = 0.95

    xres, yres, pcres = sweep(Rp, Rv)
    sum_res = sum([xres, yres, pcres])

    sum_res_list.append(sum_res)

    pbar.set_postfix({
        "sum res": f'{sum_res:.2e}'
    })

    if (sum_res < res_limit):
        print("残差收敛")
        break

    if (sum_res > 1e4):
        if (len(sum_res_list) > 20):
            if ((sum_res_list[-1] - sum_res_list[-3]) > 1e3):
                print("残差爆炸")
                break
            if (is_increasing(sum_res_list[-18:])):
                print("残差不收敛")
                break

plt.plot(np.log10(np.array(sum_res_list)))
plt.xlabel('sweep num')
plt.ylabel('log10 sum res')
plt.show()


# %%
OverflowPrevention()

Viewer(p)
plt.show(block=True)

Viewer(Vx)
plt.show(block=True)

Viewer(Vy)
plt.show(block=True)
