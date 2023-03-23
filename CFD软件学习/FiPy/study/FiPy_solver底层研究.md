
## 面上法向量

normals = FaceVariable(mesh=mesh, rank=1, value=mesh._orientedFaceNormals)



## 边界

mesh.physicalFaces["edge"].globalValue == mesh.exteriorFaces


## sweep

```python
# 缓存方程组的系数矩阵和右侧向量
eq.cacheMatrix()
eq.cacheRHSvector()

# 求解方程
# underRelaxation参数，是一种控制求解器迭代过程的技术，用于加速方程求解的过程。在每次迭代中，uR将当前迭代的解与前一次迭代的解之间进行插值，从而平滑并减缓解的变化速度。这样做可以使求解器更快地收敛，并且可以减少数值振荡和不稳定性。uR通常取值在0和1之间，表示插值的程度，一般来说，uR越小，迭代收敛的速度就越慢，但是数值稳定性会更好
# xres，代表初始值与迭代值的残差？？？反正越小越说明收敛
xres = eq.sweep(var=X, underRelaxation=uR)

# 缓存以后才能读取系数矩阵和右侧向量
xmat = eq.matrix
xrhs = eq.RHSvector

# 矩阵的对角值
ap = numerix.asarray(mat.takeDiagonal())
```

事实上，有

```python
from scipy import linalg
X.value == linalg.solve(xmat.matrix.toarray(), xrhs)
```