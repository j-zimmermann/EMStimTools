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
meshname = 'RFAblationFixedPatchstationary_mesh'

### geometry parameters###
r_tissue = 0.1
r_electrode = 0.00115
l_el = 0.03
h_tissue = 0.1
r_patch = 0.01
####################################################
##       Begin of NoteBook variables section      ##
####################################################
notebook.set("r_tissue", r_tissue)
notebook.set("r_electrode", r_electrode)
notebook.set("l_el", l_el)
notebook.set("h_tissue", h_tissue)
notebook.set("h_tip", "h_tissue-l_el")
notebook.set("h_max", "h_tissue/10.")
notebook.set("r_patch", r_patch)
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
Cylinder_1 = geompy.MakeCylinderRH("r_tissue", "h_tissue")
Vertex_1 = geompy.MakeVertex(0, 0, "h_tip")
Cylinder_2 = geompy.MakeCylinder(Vertex_1, OZ, "r_electrode", "l_el")
Sphere_1 = geompy.MakeSpherePntR(Vertex_1, "r_electrode")
Cut_1 = geompy.MakeCutList(Cylinder_1, [Cylinder_2, Sphere_1], True)
Disk_1 = geompy.MakeDiskR("r_patch", 1)
Partition_1 = geompy.MakePartition([Cut_1, Disk_1], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
listSubShapeIDs = geompy.SubShapeAllIDs(Partition_1, geompy.ShapeType["SOLID"])
OuterBoundary = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(OuterBoundary, [3, 10, 15])
geompy.DifferenceIDs(OuterBoundary, [3, 10, 15])
geompy.UnionIDs(OuterBoundary, [3, 10, 15])
Patch = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Patch, [25])
Electrode = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Electrode, [20, 27])
Tissue = geompy.CreateGroup(Partition_1, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Tissue, [1])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
geompy.addToStudy( Sphere_1, 'Sphere_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudyInFather( Partition_1, OuterBoundary, 'OuterBoundary' )
geompy.addToStudyInFather( Partition_1, Patch, 'Patch' )
geompy.addToStudyInFather( Partition_1, Electrode, 'Electrode' )
geompy.addToStudyInFather( Partition_1, Tissue, 'Tissue' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Partition_1)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( "h_max" )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 3 )
NETGEN_3D_Parameters_1.SetMinSize( 5.13575e-08 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
isDone = Mesh_1.Compute()
OuterBoundary_1 = Mesh_1.GroupOnGeom(OuterBoundary,'OuterBoundary',SMESH.FACE)
Patch_1 = Mesh_1.GroupOnGeom(Patch,'Patch',SMESH.FACE)
Electrode_1 = Mesh_1.GroupOnGeom(Electrode,'Electrode',SMESH.FACE)
Tissue_1 = Mesh_1.GroupOnGeom(Tissue,'Tissue',SMESH.VOLUME)
OuterBoundary_2 = Mesh_1.GroupOnGeom(OuterBoundary,'OuterBoundary',SMESH.NODE)
Patch_2 = Mesh_1.GroupOnGeom(Patch,'Patch',SMESH.NODE)
Electrode_2 = Mesh_1.GroupOnGeom(Electrode,'Electrode',SMESH.NODE)
Tissue_2 = Mesh_1.GroupOnGeom(Tissue,'Tissue',SMESH.NODE)
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(pwd + '/' + str(meshname) + '.med', 0, SMESH.MED_V2_2, 1, None ,1)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'


## Set names of Mesh objects
smesh.SetName(Electrode_2, 'Electrode')
smesh.SetName(Patch_2, 'Patch')
smesh.SetName(OuterBoundary_2, 'OuterBoundary')
smesh.SetName(Tissue_2, 'Tissue')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(OuterBoundary_1, 'OuterBoundary')
smesh.SetName(Tissue_1, 'Tissue')
smesh.SetName(Electrode_1, 'Electrode')
smesh.SetName(Patch_1, 'Patch')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
