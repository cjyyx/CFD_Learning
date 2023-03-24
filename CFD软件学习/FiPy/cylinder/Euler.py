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

mu = 0.01
rho = 1000.
U = 10.

sum_res_list = []


# %%
Xc, Yc = mesh.cellCenters
Vc = mesh.cellVolumes

Vx = CellVariable(mesh=mesh, name="x velocity", value=U)
Vy = CellVariable(mesh=mesh, name="y velocity", value=0.)

Vf = FaceVariable(mesh=mesh, rank=1)
Vf.setValue((Vx.faceValue, Vy.faceValue))

p = CellVariable(mesh=mesh, name="pressure", value= -Xc/100.)
pc = CellVariable(mesh=mesh, value=0.)

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
p.constrain(0., outletFace)
pc.constrain(0., outletFace)

Vx.constrain(0., cylinderFace)
Vy.constrain(0., cylinderFace)
p.faceGrad.constrain(0., cylinderFace)
pc.faceGrad.constrain(0., cylinderFace)

Vx.faceGrad.constrain(0., top_bottomFace)
Vy.faceGrad.constrain(0., top_bottomFace)
# Vx.constrain(U, top_bottomFace)
# Vy.constrain(0., top_bottomFace)
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

# Vx_Eq = \
#     UpwindConvectionTerm(coeff=Vf, var=Vx) * rho == \
#     ImplicitSourceTerm(coeff=1.0, var=p.grad[0])
# Vy_Eq = \
#     UpwindConvectionTerm(coeff=Vf, var=Vy) * rho == \
#     ImplicitSourceTerm(coeff=1.0, var=p.grad[1])

coeff = (
    1. / (
        apx.faceValue
        * mesh._faceAreas
        * mesh._cellDistances
    )
)
pc_Eq = \
    DiffusionTerm(coeff=coeff, var=pc) \
    - Vf.divergence == 0


# %%
def sweep(Rp, Rv):

    Vx_Eq.cacheMatrix()
    Vy_Eq.cacheMatrix()

    xres = Vx_Eq.sweep(var=Vx, underRelaxation=Rv)
    xmat = Vx_Eq.matrix

    yres = Vy_Eq.sweep(var=Vy, underRelaxation=Rv)
    ymat = Vy_Eq.matrix

    apx[:] = numerix.asarray(xmat.takeDiagonal())
    apy[:] = numerix.asarray(ymat.takeDiagonal())

    Vf.setValue((Vx.faceValue, Vy.faceValue))

    pc_Eq.cacheRHSvector()
    pcres = pc_Eq.sweep(var=pc)
    prhs = pc_Eq.RHSvector

    p.setValue(p + Rp * pc)

    Vx.setValue(Vx-(Vc*pc.grad[0])/apx)
    Vy.setValue(Vy-(Vc*pc.grad[1])/apx)

    return xres, yres, pcres

# %%


MaxSweep = 1000
res_limit = 0.5
sum_res = 1e10

Rp = 0.8
Rv = 0.8

pbar = tqdm(range(MaxSweep))
for i in pbar:

    # if (sum_res > 6666):
    #     Rp = Rv = 0.9
    # elif (sum_res > 3666):
    #     Rp = Rv = 0.8
    # elif (sum_res > 666):
    #     Rp = Rv = 0.4
    # elif (sum_res > 66):
    #     Rp = Rv = 0.2
    # elif (sum_res > 0.6):
    #     Rp = Rv = 0.1

    Rp=Rv=0.001

    xres, yres, pcres = sweep(Rp, Rv)
    sum_res = sum([xres, yres, pcres])

    sum_res_list.append(sum_res)

    # pbar.set_postfix({
    #     "xres": f'{xres:.2e}',
    #     "yres": f'{yres:.2e}',
    #     "pcres": f'{pcres:.2e}',
    #     "sum": f'{s:.2e}'
    # })
    pbar.set_postfix({
        "sum res": f'{sum_res:.2e}'
    })

    if (sum_res < res_limit):
        print("残差收敛")
        break
    
    if(len(sum_res_list)>10):
        if ((sum_res_list[-1] - sum_res_list[-3])>1e3):
            print("残差爆炸")
            break

# %%

plt.plot(sum_res_list[:])

# %%

sweep(0.1, 0.1)

# %%

Viewer(p)
Viewer(Vx)
Viewer(Vy)
