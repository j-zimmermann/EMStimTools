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
meshname = 'HeatRobinBC_mesh'

### geometry parameters###
w = 0.6
h = 1.0
####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("h", h)
notebook.set("w", w)
notebook.set("h_half", "h/2.0")
notebook.set("w_half", "w/2.0")
notebook.set("max_size", "h/10.")
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
Face_1 = geompy.MakeFaceHW("w", "h", 1)
geompy.TranslateDXDYDZ(Face_1, "w_half", "h_half", 0)
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Face_1, geompy.ShapeType["FACE"])
bottom = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(bottom, [6])
right = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(right, [8])
upper = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(upper, [10])
left = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(left, [3])
Plate = geompy.CreateGroup(Face_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Plate, [1])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, bottom, 'bottom' )
geompy.addToStudyInFather( Face_1, right, 'right' )
geompy.addToStudyInFather( Face_1, upper, 'upper' )
geompy.addToStudyInFather( Face_1, left, 'left' )
geompy.addToStudyInFather( Face_1, Plate, 'Plate' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Face_1)
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( "max_size" )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetMinSize( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
isDone = Mesh_1.Compute()
bottom_1 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.EDGE)
right_1 = Mesh_1.GroupOnGeom(right,'right',SMESH.EDGE)
upper_1 = Mesh_1.GroupOnGeom(upper,'upper',SMESH.EDGE)
left_1 = Mesh_1.GroupOnGeom(left,'left',SMESH.EDGE)
Plate_1 = Mesh_1.GroupOnGeom(Plate,'Plate',SMESH.FACE)
bottom_2 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.NODE)
right_2 = Mesh_1.GroupOnGeom(right,'right',SMESH.NODE)
upper_2 = Mesh_1.GroupOnGeom(upper,'upper',SMESH.NODE)
left_2 = Mesh_1.GroupOnGeom(left,'left',SMESH.NODE)
Plate_2 = Mesh_1.GroupOnGeom(Plate,'Plate',SMESH.NODE)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(pwd + '/' + str(meshname) + '.med', 0, SMESH.MED_V2_2, 1, None ,1)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Plate_1, 'Plate')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(bottom_1, 'bottom')
smesh.SetName(upper_1, 'upper')
smesh.SetName(right_1, 'right')
smesh.SetName(right_2, 'right')
smesh.SetName(upper_2, 'upper')
smesh.SetName(left_1, 'left')
smesh.SetName(bottom_2, 'bottom')
smesh.SetName(left_2, 'left')
smesh.SetName(Plate_2, 'Plate')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
