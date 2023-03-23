# %%
from fipy import Grid2D, CellVariable, FaceVariable, DiffusionTerm, ImplicitSourceTerm, Viewer

# 创建网格
nx = 50
ny = 50
dx = 0.1
dy = 0.1
mesh = Grid2D(nx=nx, ny=ny, dx=dx, dy=dy)

# 定义变量
pressure = CellVariable(name='pressure', mesh=mesh, value=0.0)
velocity = FaceVariable(name='velocity', mesh=mesh, value=(0.0, 0.0))
density = CellVariable(name='density', mesh=mesh, value=1.0)
viscosity = 0.01

# 定义边界条件
velocity.constrain((0.0, 0.0), mesh.facesLeft)
pressure.faceGrad.constrain(0.0, mesh.facesRight)
velocity.constrain((0.0, 0.0), mesh.facesTop)
velocity.constrain((1.0, 0.0), mesh.facesBottom)

# 定义方程
eq1 = DiffusionTerm(coeff=density, var=velocity) - ImplicitSourceTerm(coeff=1.0, var=pressure.grad)
eq2 = DiffusionTerm(coeff=viscosity, var=velocity.divergence)

# %%

# SIMPLE算法求解
for i in range(100):
    # 计算压力
    pressure.updateOld()
    eq1.solve(var=velocity, dt=1.0)
    res = eq1.sweep(var=velocity) + eq2.sweep(var=velocity)
    pressureCorrection = (pressure.faceGrad * mesh.faceNormals).divergence
    pressure.setValue(pressure + pressureCorrection)

    # 更新速度和密度
    velocity.setValue(velocity - pressureCorrection.grad / density)
    density.setValue(1.0 / (1.0 - velocity.divergence * dx * dy))

# %%

# 查看结果
viewer = Viewer(vars=(velocity, pressure, density), datamin=-1.0, datamax=1.0)
viewer.plot()