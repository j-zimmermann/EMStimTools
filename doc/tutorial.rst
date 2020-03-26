How to use EMStimTools
----------------------

Before giving an overview over the different parts of this package, you will find a brief tutorial on how to use the SALOME-GMSH-FEniCS workflow.


SALOME
======
SALOME is needed to create or at least process a CAD geometry and prepare the mesh.
If a MED-formatted mesh already exists, feel free to skip this part and start with the next section.
Currently, no other mesh format is supported.
Nevertheless, support for (g)msh files could be easily implemented.

EMStimTools
===========

YAML-File
^^^^^^^^^

.. warning::
	As pointed out on `github <https://github.com/yaml/pyyaml/pull/174/commits/1e0453f29583ebc96fc22ded1a78d828a4ef87c7>`_, pyyaml is not capable of understand floats in scientific notation without the dot, i.e. it cannot recognize 1e-3 but 1.0e-3!!!

You have to provide a YAML-file. 
This file wraps a Python dictionary to generate the geometry, the mesh and also the FEniCS study.
Moreover, one can specify output quantitities and formats.


ParaView
========

When the actual EMStimTools run is done, you can visualize the results (usually in xdmf-format) by ParaView.
You can also generate snapshots of the visualization state such that other researchers can understand how you generated your figures.
