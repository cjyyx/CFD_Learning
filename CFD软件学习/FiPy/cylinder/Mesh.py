# %%
import gmsh
import numpy as np
from fipy import *

isTest = True

R = 1.
N = 32

L = 36.
H = 32.

meshSize = 0.5
enMeshSize = 0.2

filename = "cylinder"

# %%


def generate_points_on_circle(n, r=1, center=(0, 0)):
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    x = center[0] + r*np.cos(angles)
    y = center[1] + r*np.sin(angles)
    return list(zip(x, y))


pL = generate_points_on_circle(N, r=R)

# %%

gmsh.initialize()

# 设置输出gmsh22版网格
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# 设置在终端打印日志
gmsh.option.setNumber("General.Terminal", 1)


gmsh.model.add("cylinder_model")
geom = gmsh.model.geo

pTL = []
lTL = []

for p in pL:
    pTL.append(geom.addPoint(
        p[0], p[1], 0
    ))

for i in range(len(pTL)-1):
    lTL.append(geom.addLine(
        pTL[i], pTL[i+1]
    ))

lTL.append(geom.addLine(
    pTL[-1], pTL[0]
))

cylinder = geom.addCurveLoop(lTL)

p1 = geom.addPoint(-L/2, -H/2, 0)
p2 = geom.addPoint(L/2, -H/2, 0)
p3 = geom.addPoint(L/2, H/2, 0)
p4 = geom.addPoint(-L/2, H/2, 0)

bottom = geom.addLine(p1, p2)
outlet = geom.addLine(p2, p3)
top = geom.addLine(p3, p4)
inlet = geom.addLine(p4, p1)

edge = geom.addCurveLoop([bottom, outlet, top, inlet])

s1 = geom.add_plane_surface([edge, cylinder])

gmsh.model.addPhysicalGroup(1, [inlet], name="inlet")
gmsh.model.addPhysicalGroup(1, [outlet], name="outlet")
gmsh.model.addPhysicalGroup(1, [top], name="top")
gmsh.model.addPhysicalGroup(1, [bottom], name="bottom")

gmsh.model.addPhysicalGroup(1, lTL, name="cylinder")

gmsh.model.addPhysicalGroup(2, [edge], name="surface")
gmsh.model.addPhysicalGroup(2, [s1], name="fluid")

geom.synchronize()


# 加密网格，参考 https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields

Box = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(Box, "VIn", enMeshSize)
gmsh.model.mesh.field.setNumber(Box, "VOut", meshSize)
gmsh.model.mesh.field.setNumber(Box, "XMin", -(L/2)*0.3)
gmsh.model.mesh.field.setNumber(Box, "XMax", (L/2)*0.88)
# gmsh.model.mesh.field.setNumber(Box, "XMax", (L/2))
gmsh.model.mesh.field.setNumber(Box, "YMin", -(H/2)*0.5)
gmsh.model.mesh.field.setNumber(Box, "YMax", (H/2)*0.5)
gmsh.model.mesh.field.setNumber(Box, "Thickness", L*0.05)

Min = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(Min, "FieldsList", [Box])

gmsh.model.mesh.field.setAsBackgroundMesh(Min)


# 设置生成四边形网格
gmsh.option.setNumber("Mesh.RecombineAll", 1)

gmsh.model.mesh.generate(2)

# 可视化网格
gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
gmsh.option.setNumber('Mesh.Points', 1)
gmsh.fltk.run()

# 输出
gmsh.write("{}.msh2".format(filename))
gmsh.finalize()

# %%

# 下面测试网格是否可用

mesh = Gmsh2D("{}.msh2".format(filename))

T_in = 128.
T_out = 64.
T_mid = 0.

T = CellVariable(
    name="temperature",
    mesh=mesh,
    value=T_out
)

T.constrain(T_in, mesh.physicalFaces["inlet"])
T.constrain(T_out, mesh.physicalFaces["outlet"])
T.constrain(T_mid, mesh.physicalFaces["cylinder"])

T.faceGrad.constrain(
    0.,
    mesh.physicalFaces["top"] |
    mesh.physicalFaces["bottom"]
)

viewer = Viewer(vars=T)

eq = DiffusionTerm()
eq.solve(var=T)

viewer.plot()
