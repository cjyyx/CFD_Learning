# %%

import matplotlib.pyplot as plt
import numpy as np
import gmsh
import os

# %%

x1 = 5.
x2 = 7.
y1 = 5.
y2 = 8.

a = 2.

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("clip_model")
geom = gmsh.model.geo

p1 = geom.addPoint(0, 0, 0)
p2 = geom.addPoint(-x1, 0, 0)
p3 = geom.addPoint(-x1, y1, 0)
p4 = geom.addPoint(x2, y1, 0)
p5 = geom.addPoint(x2, -y2, 0)
p6 = geom.addPoint(0, -y2, 0)

l1 = geom.addLine(p1, p2)
l2 = geom.addLine(p2, p3)
l3 = geom.addLine(p3, p4)
l4 = geom.addLine(p4, p5)
l5 = geom.addLine(p5, p6)
l6 = geom.addLine(p6, p1)

edge = geom.addCurveLoop([l1, l2, l3, l4, l5, l6])

p7 = geom.addPoint(-a/2, y1/2-a/2, 0)
p8 = geom.addPoint(-a/2, y1/2+a/2, 0)
p9 = geom.addPoint(a/2, y1/2+a/2, 0)
p10 = geom.addPoint(a/2, y1/2-a/2, 0)

l7 = geom.addLine(p7, p8)
l8 = geom.addLine(p8, p9)
l9 = geom.addLine(p9, p10)
l10 = geom.addLine(p10, p7)

cool = geom.addCurveLoop([l7, l8, l9, l10])

p11 = geom.addPoint(x2/2-a/2, -y2/2-a/2, 0)
p12 = geom.addPoint(x2/2-a/2, -y2/2+a/2, 0)
p13 = geom.addPoint(x2/2+a/2, -y2/2+a/2, 0)
p14 = geom.addPoint(x2/2+a/2, -y2/2-a/2, 0)

l11 = geom.addLine(p11, p12)
l12 = geom.addLine(p12, p13)
l13 = geom.addLine(p13, p14)
l14 = geom.addLine(p14, p11)

hot = geom.addCurveLoop([l11, l12, l13, l14])

s1 = geom.addPlaneSurface([edge, cool, hot])

geom.synchronize()

gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6], name="edge")
gmsh.model.addPhysicalGroup(1, [l7, l8, l9, l10], name="cool")
gmsh.model.addPhysicalGroup(1, [l11, l12, l13, l14], name="hot")
gmsh.model.addPhysicalGroup(2, [edge], name="clip")


# %%

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)

gmsh.model.mesh.setOrder(2)
gmsh.model.mesh.generate(2)


filename = 'clip'
gmsh.write(filename+".msh")

print("generate ok")

gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
gmsh.option.setNumber('Mesh.Points', 1)

gmsh.fltk.run()

gmsh.finalize()

os.system("gmsh -2 -format msh2 {}.msh -o {}.msh2".format(
    filename, filename
))
