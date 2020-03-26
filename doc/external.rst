External depedencies
--------------------

The package works with:

- Salome 8.3.0 (to be downloaded `here <http://salome-platform.org/downloads/previous-versions/salome-v8.3.0>`_)
  Salome requires on Ubuntu 18.04. the manual installation of two libraries

  #. libicu55 (for instance to be downloaded `here <https://packages.ubuntu.com/de/xenial/amd64/libicu55/download>`_)
  #. libpng12 (see this `site <https://packages.ubuntu.com/de/xenial/amd64/libpng12-0/download>`_)

  Furthermore, net-tools and libpcre3 should be installed:
  
  #. run :code:`sudo apt-get install net-tools`
  #. run :code:`sudo apt-get install libpcre3-dev`
    

- GMSH 4.0.0 SDK (to be downloaded `here <http://gmsh.info/bin/Linux/>`_) 

  This requires on Ubuntu 18.04. the manual installation of libgfortran3, e.g. with :code:`sudo apt-get install libgfortran3`
  To use the SDK, set a `GMSH_DIRECTORY`, i.e. the directory where the downloaded SDK archive is located. Then add `$GMSH_DIRECTORY/bin` to your `PATH` variable and `$GMSH_DIRECTORY\lib` to your `PYTHONPATH` variable.
- matplotlib, e.g. with :code:`sudo apt-get install python3-matplotlib`
- Sphinx (:code:`sudo apt-get install python-sphinx`)
