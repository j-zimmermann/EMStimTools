EMStimTools
============

Work on an open-source tool for EM stimulation chambers and other tools.
The entire package is developed and tested on Ubuntu 18.04. 

Installation
------------

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

To install the external packages, just run:

.. code::

	sudo install.sh 

This should give you all dependencies and correct setting of the environment.

Install EMStimTools (SALOME and GMSH there)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simply run:

.. code::

	python3 setup.py install

or if you want to have a local installation:

.. code::

	python3 setup.py install --user


Third party code
----------------

The files 'meshconvert.py' and 'xml_writer.py' from the site-package 'dolfin_utils' of dolfin. They can be found in 'emstimtools/utils' and the first file was modified slightly (change indicated in file header). It is licensed under the GNU Lesser General Public License (LGPL) version 3, which can be found in the file LICENSE-3RD-PARTY.txt.

External dependencies
---------------------

The package works with:

- Python 3

- Salome 8.3.0 (to be downloaded here_)

.. _here: http://salome-platform.org/downloads/previous-versions/salome-v8.3.0

- GMSH SDK 4.0.2 (to be download from gmsh_ site) 

.. _gmsh: http://gmsh.info/bin/Linux/gmsh-4.0.2-Linux64-sdk.tgz

- FEniCS 2019.1.0 (from official Ubuntu PPA, see details on FEniCS site_)

.. _site: https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages

Optional:

- ParaView > V5.5 (https://www.paraview.org/) for postprocessing (XDMF support needed).

.. note:: For SALOME and GMSH certain manual installations are needed. Those are specified in the documentation.


Overview
--------

The package takes a user-prepared geometry and input file, where all details about the (multi-)physics model are described.
Then a fully-automated workflow yields the desired output data and files, which can be readily visualized in e.g. ParaView.

Documentation
-------------

This code is documented by the help of Sphinx.
For a good user experience, install sphinx on your machine by (example for Ubuntu 18.04)

.. code::

        sudo apt-get install python3-sphinx
        sudo apt-get install latexmk

and run 

.. code::

	make latexpdf

in the `doc` directory. It will give you a good overview of the capabilities of this tool.
Also, certain bugs are explained in case an error occurs during execution.

If you do not run Ubuntu 16 or 18, you can use Docker to generate the LaTex file locally.
The procedure here is:

1. Set up Docker on your machine
2. Run `docker build . -t=emstimtools:0.1.3.dev0` (or use a different name and tag with `-t=name:tag` format) 
3. Run `docker run -ti -v $(pwd):/home/ -p 8080:8080 emstimtools:0.1.3.dev0 bash` in EMStimTools directory.
4. A bash shell will open. Execute `cd /home/doc` and then execute `make latexpdf`. The documentation will be generated. 

When executing a study, use the same docker image, but proceed as follows:

1. Run `docker run -ti -v $(pwd):/home/ -p 8080:8080 emstimtools:0.1.3.dev0 bash` in the directory, where the Python script for your study is located.
2. A bash shell will open. Execute `cd /home` and then execute the Python file of your study by running `python3 study.py`. 


Examples
--------

In the examples directory, a few published as well as unpublished examples are given. Please feel free to add some on your own!
