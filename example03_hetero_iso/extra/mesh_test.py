#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:37:05 2020

@author: ranjeet
"""


from dolfin import *

p0 = Point(-1,-1,-1)
p1 = Point(1,1,1)

mesh = BoxMesh.create([p0,p1],[4,4,4], CellType.Type.hexahedron)

p3 = Point(0,0)
p4 = Point(1,1)

mesh = RectangleMesh.create([p3,p4],[4,4], CellType.Type.quadrilateral)
plot(mesh)
ele = FiniteElement("Q-", mesh.ufl_cell(),1,1)
mesh2 = UnitSquareMesh(4,4)
plot(mesh)