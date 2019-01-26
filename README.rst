EMStimTools
============

Work on an open-source tool for EM stimulation chambers and other tools.
The entire package is developed and tested on Ubuntu 18.04. 

Installation
------------

Simply run:

.. code::

	python3 setup.py install

or if you want to have a local installation:

.. code::

	python3 setup.py install --user


Third party code
----------------

The files 'meshconvert.py' and 'xml_writer.py' from the site-package 'dolfin_utils' of dolfin. They can be found in 'emstimtools/utils' and the first file was modified slightly (change indicated in file header). It is licensed under the GNU Lesser General Public License (LGPL) version 3, which can be found in the file LICENSE-3RD-PARTY.txt.

External depedencies
--------------------

The package works with:

- Python 3

- Salome 8.3.0 (to be downloaded here_)

.. _here: http://salome-platform.org/downloads/previous-versions/salome-v8.3.0

- GMSH SDK 4.0.2 (to be download from gmsh_ site) 

.. _gmsh: http://gmsh.info/bin/Linux/gmsh-4.0.2-Linux64-sdk.tgz

- FEniCS 2018.1.0 (from official Ubuntu PPA, see details on FEniCS site_)

.. _site: https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages

Optional:

- ParaView > V5.5 (https://www.paraview.org/) for postprocessing (XDMF support needed).


Overview
--------

The package takes a user-prepared geometry and input file, where all details about the (multi-)physics model are described.
Then a fully-automated workflow yields the desired output data and files, which can be readily visualized in e.g. ParaView.

Documentation
-------------

This code is more or less documented by the help of Sphinx.
For a good user experience, install sphinx of your machine and run 

.. code::

	make latexpdf

in the doc directory. It will give you a good overview of the capabilities of this tool.
Also, certain bugs ('features') are explained in case an error occurs during execution.

Examples
--------

In the examples directory, a few published as well as unpublished examples are given. Please feel free to add some on your own!
