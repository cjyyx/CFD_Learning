#! https://zhuanlan.zhihu.com/p/614634508

# FiPy常用API

官方文档(https://www.ctcms.nist.gov/fipy/documentation/API.html)

## 网格

### [AbstractMesh](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.abstractMesh.AbstractMesh)

基类：object


#### [\_\_add__(other)](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.abstractMesh.AbstractMesh.__add__)

用于平移网格

#### .aspect2D

2D网格，y与x纵横比

#### .cellCenters

返回`fipy.variables.cellVariable.CellVariable`，形状为`(2, cellNum)`

cellNum为网格单元的数量

网格中心的坐标

可通过`numpy.array`转化成`numpy`数组

#### .cellDistanceVectors

返回`numpy.ndarray`，形状为`(2, faceNum)`

faceNum是网格面数

该属性表示网格中心与相邻网格中心的矢量差

#### .cellFaceIDs

返回`numpy.ma.core.MaskedArray`，形状为`(3, cellNum)`

该属性表示网格单元所对应的面（或边）的全局面（或边）ID

#### .cellToFaceDistanceVectors

返回`numpy.ma.core.MaskedArray`

该属性表示网格单元和其相邻面（或边）之间的距离向量

对于2维网格，形状为`(2, 2, faceNum)`。`cellToFaceDistanceVectors[0][0][i]`表示该单元与其第`i`个边之间的距离向量的第一个分量，`cellToFaceDistanceVectors[0][1][i]`则表示该距离向量的第二个分量。（ChatGPT说的，不知道对不对）

#### .cellVolumes

返回`numpy.ndarray`，形状为`(cellNum,)`


#### .extents

返回`dict`

表示网格坐标的范围

#### .exteriorFaces

返回`fipy.variables.faceVariable.FaceVariable`，形状为`(faceNum,)`


#### .faceCenters

返回`fipy.variables.faceVariable.FaceVariable`，形状为`(2, faceNum)`

网格对象的每个面的中心坐标

#### .faceNormals

返回`numpy.ma.core.MaskedArray`

对于2维网格，形状为'(2, faceNum)'。此时，`(faceNormals[0][i], faceNormals[1][i])`表示全局边ID为`i`的边的法向量。

#### getNearestCell(point)

用于查找距离某个给定位置最近的网格单元

#### [.x.globalValue & .y.globalValue](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.abstractMesh.AbstractMesh.x)

返回`numpy.ndarray`，表示网格中心的坐标



### [Mesh](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.mesh.Mesh)

基类：AbstractMesh

#### [\_\_mul__(factor)](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.mesh.Mesh.__mul__)

拉伸网格


### [Mesh2D](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##module-fipy.meshes.mesh2D)

基类：Mesh


### [Gmsh2D](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html##fipy.meshes.gmshMesh.Gmsh2D)

基类：Mesh2D

#### 导入`.msh`文件

详见(https://zhuanlan.zhihu.com/p/614003913)

####通过标签获取网格面

详见(https://zhuanlan.zhihu.com/p/614003913)

#### 增加网格密度

官方文档内提了几句，没看懂

## 变量

## 边界条件

## 求解项

## Viewer

### [AbstractViewer](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.viewers.html##fipy.viewers.viewer.AbstractViewer)

#### .fig

返回`matplotlib.figure.Figure`

#### _plot()

基于当前的变量值绘制，更新`fig`属性

### [MultiViewer](https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.viewers.html##module-fipy.viewers.multiViewer)

基类：AbstractViewer

<!-- 一般当代码如下时，会返回`MultiViewer`类 -->

#### .viewers

返回普通`viewer`类的列表

```python
viewer = Viewer(
    vars=(pressure, xVelocity, yVelocity)
    )

## 绘制pressure 的图
viewer.viewers[0].plot()

## 绘制xVelocity 的图
viewer.viewers[1].plot()

## 绘制yVelocity 的图
viewer.viewers[2].plot()
```