# %%
# 导入必要的库
from fipy import CellVariable, FaceVariable, Grid2D, DiffusionTerm, Viewer
from fipy.tools import numerix
from fipy.variables.faceGradVariable import _FaceGradVariable

# %%
# 设置流体参数与仿真参数
L = 1.0
N = 50
dL = L / N
viscosity = 1
U = 1.
# 0.8 for pressure and 0.5 for velocity are typical relaxation values for SIMPLE
pressureRelaxation = 0.8
velocityRelaxation = 0.5

mesh = Grid2D(nx=N, ny=N, dx=dL, dy=dL)

# %%

""" 
pressure, xVelocity, yVelocity: 分别是代表压力、x方向速度、y方向速度的 CellVariable 对象，这些变量在网格的单元格上定义。

pressureCorrection: 代表压力修正的 CellVariable 对象，也在网格的单元格上定义。

velocity: 代表速度的 FaceVariable 对象，这个变量在网格的面上定义。

xVelocityEq, yVelocityEq: 分别是代表 x 方向速度和 y 方向速度的方程。它们通过 DiffusionTerm 类来描述黏性的贡献，并通过 pressure.grad.dot([1., 0.]) 和 pressure.grad.dot([0., 1.]) 来描述压力的贡献。

ap: 代表系数 a_p 的 CellVariable 对象，它在网格的单元格上定义，并初始化为1。

coeff: 代表系数的 FaceVariable 对象，这些系数在网格的面上定义。它通过 1./ ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances 计算得出。

pressureCorrectionEq: 代表压力修正方程的 DiffusionTerm，它描述了压力修正随时间的变化。它通过 velocity.divergence 来描述速度的贡献。

volume, contrvolume: 分别是代表单元体积和对应面的面积的 CellVariable 对象，这些变量在网格的单元格和面上定义。其中，contrvolume 是 volume 的算术平均数，用于计算通量。

xVelocity.constrain(0., mesh.facesRight | mesh.facesLeft | mesh.facesBottom): 在网格的右侧、左侧和底部边界处，x 方向速度被固定为零。

xVelocity.constrain(U, mesh.facesTop): 在网格的顶部边界处，x 方向速度被固定为一个常数 U。

yVelocity.constrain(0., mesh.exteriorFaces): 在网格的外部面上，y 方向速度被固定为零。

pressureCorrection.constrain(0., mesh.facesLeft & (Y < dL)): 在网格的左侧且位于 Y < dL 的面上，压力修正被固定为零。

 """

# 在网格点上声明变量
pressure = CellVariable(mesh=mesh, name='pressure')
xVelocity = CellVariable(mesh=mesh, name='X velocity')
yVelocity = CellVariable(mesh=mesh, name='Y velocity')

pressureCorrection = CellVariable(mesh=mesh)
velocity = FaceVariable(mesh=mesh, rank=1)

xVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.grad.dot([1., 0.])
yVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.grad.dot([0., 1.])

ap = CellVariable(mesh=mesh, value=1.)
coeff = 1. / ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume = volume.arithmeticFaceValue

xVelocity.constrain(0., mesh.facesRight | mesh.facesLeft | mesh.facesBottom)
xVelocity.constrain(U, mesh.facesTop)
yVelocity.constrain(0., mesh.exteriorFaces)
X, Y = mesh.faceCenters
pressureCorrection.constrain(0., mesh.facesLeft & (Y < dL))

# %%

viewer = Viewer(vars=(pressure, xVelocity, yVelocity, velocity),
                xmin=0., xmax=1., ymin=0., ymax=1., colorbar='vertical', scale=5)

# %%

""" 
针对 $x$ 和 $y$ 方向上的速度分量，利用 FiPy 中的 DiffusionTerm 类定义扩散项，利用CellVariable 定义速度和压力，以及 FaceVariable 定义对应面的速度（rank=1表示变量是向量）；

通过 cacheMatrix 方法缓存系数矩阵（$x$ 方向和 $y$ 方向分别计算）；

利用 sweep 方法求解扩散方程，其中 var 参数表示求解的变量，underRelaxation 参数表示松弛因子；

利用系数矩阵的对角线更新系数 ap，该系数是在压力修正方程中使用的；

利用 Rhie-Chow 校正，更新面速度（arithmeticFaceValue 表示面变量的算术平均值）；

利用 cacheRHSvector 方法缓存右端项向量；

利用 sweep 方法求解压力修正方程；

利用修正的压力值更新压力；

利用修正的压力值更新速度，其中 $\text{grad}$ 表示梯度算子。
 """

""" 
xVelocityEq.cacheMatrix(): 缓存 $x$ 速度方程的系数矩阵
xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation): 执行一次 $x$ 速度方程的迭代计算
xmat = xVelocityEq.matrix: 获取 $x$ 速度方程的系数矩阵
yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation): 执行一次 $y$ 速度方程的迭代计算
ap[:] = -numerix.asarray(xmat.takeDiagonal()): 从矩阵对角线更新 $ap$ 系数
presgrad = pressure.grad: 计算压力梯度
facepresgrad = _FaceGradVariable(pressure): 计算面上的压力梯度
velocity[0]: 更新 $x$ 速度面上的值，通过星号值与 Rhie-Chow 修正
velocity[1]: 更新 $y$ 速度面上的值，通过星号值与 Rhie-Chow 修正
velocity[..., mesh.exteriorFaces.value] = 0.: 边界处速度为零
velocity[0, mesh.facesTop.value] = U: 顶部 $x$ 速度设为常数 $U$
pressureCorrectionEq.cacheRHSvector(): 缓存压力修正方程的右侧向量
pres = pressureCorrectionEq.sweep(var=pressureCorrection): 解压力修正方程
rhs = pressureCorrectionEq.RHSvector: 获取压力修正方程的右侧向量
pressure.setValue(pressure + pressureRelaxation * pressureCorrection): 更新压力值
xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes): 根据压力修正值更新 $x$ 速度值
yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / ap * mesh.cellVolumes): 根据压力修正值更新 $y$ 速度值
 """

# xVelocityEq.cacheMatrix()是将线性方程组系数矩阵缓存到内存中，以便于在多次迭代中重复使用。这样可以节省计算时间，特别是对于大型线性方程组的求解来说，缓存系数矩阵是非常重要的优化手段。
xVelocityEq.cacheMatrix()

# 返回求解后的残差。这是一个表示方程误差的向量，用于判断解是否收敛。残差通过将已更新的解决方案与先前的解决方案之间的差异加权来计算。其中，权重是通过指定的relaxation参数计算得出的。
xres = xVelocityEq.sweep(var=xVelocity,
                         underRelaxation=velocityRelaxation)
xmat = xVelocityEq.matrix

yres = yVelocityEq.sweep(var=yVelocity,
                         underRelaxation=velocityRelaxation)

# update the ap coefficient from the matrix diagonal
ap[:] = -numerix.asarray(xmat.takeDiagonal())

# update the face velocities based on starred values with the
# Rhie-Chow correction.
# cell pressure gradient
presgrad = pressure.grad
# face pressure gradient
facepresgrad = _FaceGradVariable(pressure)

velocity[0] = xVelocity.arithmeticFaceValue \
    + contrvolume / ap.arithmeticFaceValue * \
    (presgrad[0].arithmeticFaceValue-facepresgrad[0])
velocity[1] = yVelocity.arithmeticFaceValue \
    + contrvolume / ap.arithmeticFaceValue * \
    (presgrad[1].arithmeticFaceValue-facepresgrad[1])
velocity[..., mesh.exteriorFaces.value] = 0.
velocity[0, mesh.facesTop.value] = U

# solve the pressure correction equation
pressureCorrectionEq.cacheRHSvector()
# left bottom point must remain at pressure 0, so no correction
pres = pressureCorrectionEq.sweep(var=pressureCorrection)
rhs = pressureCorrectionEq.RHSvector

# update the pressure using the corrected value
pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
# update the velocity using the corrected pressure
xVelocity.setValue(xVelocity - pressureCorrection.grad[0] /
                   ap * mesh.cellVolumes)
yVelocity.setValue(yVelocity - pressureCorrection.grad[1] /
                   ap * mesh.cellVolumes)

# %%

viewer.viewers[0].plot("temp.png")
