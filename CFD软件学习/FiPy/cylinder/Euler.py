# %%
from scipy.linalg import lu
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

# mu = 0.01
# rho = 1000.

mu = 0.1
rho = 1.

U = 10.

sum_res_list = []


# %%
Xc, Yc = mesh.cellCenters
Vc = mesh.cellVolumes
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes)
Vcf = volume.arithmeticFaceValue

Vx = CellVariable(mesh=mesh, name="x velocity", value=U)
Vy = CellVariable(mesh=mesh, name="y velocity", value=0.)

Vf = FaceVariable(mesh=mesh, rank=1)
Vf.setValue((Vx.faceValue, Vy.faceValue))

p = CellVariable(mesh=mesh, name="pressure", value=-Xc/1000.)
# p = CellVariable(mesh=mesh, name="pressure", value= 0.)
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
#     DiffusionTerm(coeff=mu, var=Vx) - \
#     ImplicitSourceTerm(coeff=1.0, var=p.grad[0])
# Vy_Eq = \
#     DiffusionTerm(coeff=mu, var=Vy) - \
#     ImplicitSourceTerm(coeff=1.0, var=p.grad[1])

# Vx_Eq = \
#     UpwindConvectionTerm(coeff=Vf, var=Vx) * rho \
#     + ImplicitSourceTerm(coeff=1.0, var=p.grad[0])
# Vy_Eq = \
#     UpwindConvectionTerm(coeff=Vf, var=Vy) * rho \
#     + ImplicitSourceTerm(coeff=1.0, var=p.grad[1])

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
    facepresgrad = _FaceGradVariable(p)
    Vf[0] = Vx.faceValue + Vcf / apx.faceValue * \
        (presgrad[0].faceValue-facepresgrad[0])
    Vf[1] = Vy.faceValue + Vcf / apx.faceValue * \
        (presgrad[1].faceValue-facepresgrad[1])

    pcres = pc_Eq.sweep(var=pc)

    p.setValue(p + Rp * pc)

    Vx.setValue(Vx-(Vc*pc.grad[0])/apx)
    Vy.setValue(Vy-(Vc*pc.grad[1])/apx)

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


MaxSweep = 100
res_limit = 1e-5
sum_res = 1e10

Rp = 0.8
Rv = 0.5

pbar = tqdm(range(MaxSweep))
for i in pbar:

    if sum_res > 1e2:
        Rp = 0.8
        Rv = 0.5
    elif valueRange(sum_res, 28., 1e2):
        Rp = 0.9
        Rv = 0.6
    elif valueRange(sum_res, 10., 28.):
        Rp = 0.95
        Rv = 0.8
    elif valueRange(sum_res, 0.5, 10.):
        Rp = 0.98
        Rv = 0.92
    elif valueRange(sum_res, 0, 0.5):
        Rp = 0.998
        Rv = 0.95

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

    if (sum_res > 1e4):
        if (len(sum_res_list) > 20):
            # if ((sum_res_list[-1] - sum_res_list[-3]) > 1e3):
            #     print("残差爆炸")
                # break
            if (is_increasing(sum_res_list[-18:])):
                print("残差不收敛")
                break

plt.plot(np.log10(np.array(sum_res_list)))
plt.xlabel('sweep num')
plt.ylabel('log10 sum res')
plt.show()

# %%

n = len(sum_res_list)
y = np.array(sum_res_list)
y = np.log10(y)

a = 12720
b = n

plt.plot(range(a,b), y[a:b])
plt.xlabel('sweep num')
plt.ylabel('log10 sum res')
plt.show()


# %%

sweep(1e-2, 1e-2)

# %%

OverflowPrevention()
Viewer(p)
Viewer(Vx)
Viewer(Vy)

# %%
OverflowPrevention()
Viewer(p).fig.savefig('p.png', dpi=600)
Viewer(Vx).fig.savefig('Vx.png', dpi=600)
Viewer(Vy).fig.savefig('Vy.png', dpi=600)

# %%

Rv = 0.2

Vx_Eq.cacheMatrix()
Vx_Eq.cacheRHSvector()
xres = Vx_Eq.sweep(var=Vx, underRelaxation=Rv)
xmat = Vx_Eq.matrix
xrhs = Vx_Eq.RHSvector
apx[:] = numerix.asarray(xmat.takeDiagonal())

Vy_Eq.cacheMatrix()
Vy_Eq.cacheRHSvector()
yres = Vy_Eq.sweep(var=Vy, underRelaxation=Rv)
ymat = Vy_Eq.matrix
yrhs = Vy_Eq.RHSvector
apy[:] = numerix.asarray(ymat.takeDiagonal())

# %%

""" 
当 Rv = 1. 时，有
xrhs == (-p.grad[0].value * Vc)
yrhs == (-p.grad[1].value * Vc)

恒有
xmat.matrix.data == ymat.matrix.data

降低`RuntimeError: Factor is exactly singular`风险的方法
1. 使用更精细网格
2. 降低 Rv
3. 降低 rho

残差爆炸出现原因
1. Rp 过高
2. 流场中出现奇点，形成正反馈

降低梯度爆炸风险的方法
1. 降低 Rp
2. 设置阈值，避免流场变量溢出

最终残差处于缓慢下降时，使用更高的 Rp,Rv 效果意外得好。
 """

# %%
sparse_matrix = Vx_Eq.matrix.matrix

# convert sparse matrix to dense matrix
dense_matrix = sparse_matrix.toarray()

# compute the determinant of the matrix
determinant = np.linalg.det(dense_matrix)
print("Determinant:", determinant)

# %%


# 构造稀疏矩阵A
A = Vx_Eq.matrix.matrix

# 将稀疏矩阵转化为密集矩阵
A_dense = A.todense()

# 进行LU分解
P, L, U = lu(A_dense)

# 计算行列式值
det_A = np.linalg.det(U) * ((-1) ** P[np.diag_indices_from(U)].sum())

# 输出行列式值
print('determinant of A:', det_A)

