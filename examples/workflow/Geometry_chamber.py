# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import os
import sys #
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
pwd = os.getcwd()
print("Current path is {}".format(pwd))
notebook = salome_notebook.NoteBook(theStudy) #
sys.path.insert(0, pwd)
meshname = 'study_mobini2017_mesh'

### geometry parameters###
r_el = 0.0005
r_well = 0.01689
h_well = 0.01733
h_el = 0.001
dist_el = 0.022
level = 0.8
l_el = 0.05
ratio = 0.44
####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("r_well", r_well)
notebook.set("h_well", h_well)
notebook.set("level", level)
notebook.set("r_el", r_el)
notebook.set("h_el", h_el)
notebook.set("dist_el", dist_el)
notebook.set("h_liquid", "h_well*level")
notebook.set("h_air", "h_well*1.5")
notebook.set("l_el", l_el)
notebook.set("ratio", ratio)
notebook.set("l_part1", "l_el*ratio-1.5*r_el")
notebook.set("l_part2", "l_el-l_part1-1.5*r_el")
notebook.set("dist_el_x_n", "-l_part1/2.+0.75*r_el")
notebook.set("dist_el_y", "dist_el/2.+r_el")
notebook.set("dist_el_y_n", "-dist_el_y")
notebook.set("h_el_pipe", "h_el+r_el*1.5")
notebook.set("dist_el_x_pipe_n", "dist_el_x_n-r_el*1.5")
notebook.set("h_el_pipe_h", "h_el+1.5*(r_el-r_el*0.52532198881)")
notebook.set("dist_el_x_pipe_h_n", "dist_el_x_n-1.5*r_el*0.85090352453")
####################################################
##        End of NoteBook variables section       ##
####################################################
###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Cylinder_1 = geompy.MakeCylinderRH("r_well", "h_liquid")
Vertex_1 = geompy.MakeVertex(0, 0, "h_liquid")
Cylinder_2 = geompy.MakeCylinder(Vertex_1, OZ, "r_well", "h_air")
Vertex_2 = geompy.MakeVertex("dist_el_x_n", "dist_el_y", "h_el")
Cylinder_3 = geompy.MakeCylinder(Vertex_2, OX, "r_el", "l_part1")
Vertex_3 = geompy.MakeVertex("dist_el_x_pipe_h_n", "dist_el_y", "h_el_pipe_h")
Vertex_4 = geompy.MakeVertex("dist_el_x_pipe_n", "dist_el_y", "h_el_pipe")
Arc_1 = geompy.MakeArc(Vertex_2, Vertex_3, Vertex_4)
Cylinder_4 = geompy.MakeCylinder(Vertex_4, OZ, "r_el", "l_part2")
[Face_1] = geompy.SubShapes(Cylinder_4, [12])
geomObj_1 = geompy.MakePipe(Face_1, Arc_1)
listSubShapeIDs = geompy.SubShapeAllIDs(geomObj_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(geomObj_1, geompy.ShapeType["VERTEX"])
Pipe_1 = geompy.MakePipe(Face_1, Arc_1)
Electrode_1 = geompy.MakePartitionNonSelfIntersectedShape([Cylinder_3, Cylinder_4, Pipe_1], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0, True)
[geomObj_2,geomObj_3,geomObj_4,geomObj_5,geomObj_6,geomObj_7,geomObj_8] = geompy.SubShapeAll(Electrode_1, geompy.ShapeType["FACE"])
Vertex_5 = geompy.MakeVertex("dist_el_x_n", "dist_el_y_n", "h_el")
Cylinder_5 = geompy.MakeCylinder(Vertex_5, OX, "r_el", "l_part1")
Vertex_6 = geompy.MakeVertex("dist_el_x_pipe_h_n", "dist_el_y_n", "h_el_pipe_h")
Vertex_7 = geompy.MakeVertex("dist_el_x_pipe_n", "dist_el_y_n", "h_el_pipe")
Arc_2 = geompy.MakeArc(Vertex_5, Vertex_6, Vertex_7)
Cylinder_6 = geompy.MakeCylinder(Vertex_7, OZ, "r_el", "l_part2")
[Face_2] = geompy.SubShapes(Cylinder_6, [12])
Pipe_2 = geompy.MakePipe(Face_2, Arc_2)
Electrode_2 = geompy.MakePartitionNonSelfIntersectedShape([Cylinder_5, Cylinder_6, Pipe_2], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0, True)
Partition_1 = geompy.MakePartitionNonSelfIntersectedShape([Cylinder_1, Cylinder_2], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0, True)
Cut_1 = geompy.MakeCutList(Partition_1, [Electrode_1, Electrode_2], True)
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["SOLID"])
listSameIDs = geompy.GetSameIDs(Cut_1, geomObj_2)
listSameIDs = geompy.GetSameIDs(Cut_1, geomObj_4)
listSameIDs = geompy.GetSameIDs(Cut_1, geomObj_7)
listSameIDs = geompy.GetSameIDs(Cut_1, geomObj_8)
Contact2 = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Contact2, [31, 41, 51, 74, 21, 64])
Contact1 = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Contact1, [26, 36, 46, 53, 69, 76])
OuterIsolation = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(OuterIsolation, [57, 4, 62, 11])
Medium = geompy.CreateGroup(Cut_1, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Medium, [2])
Air = geompy.CreateGroup(Cut_1, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Air, [55])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Cylinder_3, 'Cylinder_3' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Arc_1, 'Arc_1' )
geompy.addToStudy( Cylinder_4, 'Cylinder_4' )
geompy.addToStudyInFather( Cylinder_4, Face_1, 'Face_1' )
geompy.addToStudy( Pipe_1, 'Pipe_1' )
geompy.addToStudy( Electrode_1, 'Electrode_1' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Cylinder_5, 'Cylinder_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Arc_2, 'Arc_2' )
geompy.addToStudy( Cylinder_6, 'Cylinder_6' )
geompy.addToStudyInFather( Cylinder_6, Face_2, 'Face_2' )
geompy.addToStudy( Pipe_2, 'Pipe_2' )
geompy.addToStudy( Electrode_2, 'Electrode_2' )
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Cut_1, Contact2, 'Contact2' )
geompy.addToStudyInFather( Cut_1, Contact1, 'Contact1' )
geompy.addToStudyInFather( Cut_1, OuterIsolation, 'OuterIsolation' )
geompy.addToStudyInFather( Cut_1, Medium, 'Medium' )
geompy.addToStudyInFather( Cut_1, Air, 'Air' )

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
