#    EMStimTools is a package that provides a SALOME-GMSH-FEniCS workflow to solve problems
#    related to electromagnetic stimulation of for instance biological tissue or cell cultures.
#
#    Copyright (C) 2018 Julius Zimmermann, julius.zimmermann[AT]uni-rostock.de
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
import gmsh
import os
import errno
import dolfin as d
from . import meshconvert
from .extractsubdomains import extractSubdomainsandFacets


def convertMEDfile(med_file):
    """
    This function opens an MED file by means of gmsh
    and saves it in .msh2 format.

    :param str med_file: has to end with .med

    """
    filename, ext = os.path.splitext(med_file)
    if ext != ".med":
        print("You must provide a med file!")
        sys.exit(0)
    # initialize gmsh
    gmsh.initialize()
    # gmsh writes out messages
    gmsh.option.setNumber("General.Terminal", 1)
    # add a new model
    gmsh.model.add("new_model")
    # open med file
    gmsh.open(med_file)
    # need dirty hack to write in correct format
    gmsh.write(filename + ".msh2")
    gmsh.finalize()
    return


def convertmshfile(msh_file, dimension):
    """
    convert msh file by using meshconvert (distributed in dolfin_utils)

    :param str msh_file: msh2 format
    :param int dimension: dimension of model

    """

    filename, ext = os.path.splitext(msh_file)
    if ext != ".msh2":
        print("You must provide a msh file in format msh2!")
        sys.exit(0)
    xml_file = filename + '.xml'
    iformat = 'gmsh'
    meshconvert.convert2xml(msh_file, xml_file, iformat=iformat)
    subdomains, facets = extractSubdomainsandFacets(msh_file, dimension)
    return subdomains, facets


def createh5file(xml_file, logger=None):
    """
    use dolfin to create hdf5 file.

    :param str xml_file: xml file prepared by meshconvert (command line: dolfin-convert)
    :param logger: logger object to provide log output, if not then print log

    """
    filename, ext = os.path.splitext(xml_file)
    if ext != ".xml":
        print("You must provide a xml file!")
        sys.exit(0)
    mesh_dir = "mesh/"
    if not os.path.exists(mesh_dir):
        try:
            os.makedirs(mesh_dir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    mesh = d.Mesh(filename + ".xml")
    cd = d.MeshFunction('size_t', mesh, filename + "_physical_region.xml")
    fd = d.MeshFunction('size_t', mesh, filename + "_facet_region.xml")
    hdf = d.HDF5File(mesh.mpi_comm(), mesh_dir + filename + ".h5", "w")
    hdf.write(mesh, "/mesh")
    hdf.write(cd, "/subdomains")
    hdf.write(fd, "/facets")
    if logger is None:
        print("Converted mesh to .h5 format")
    else:
        logger.debug("Converted mesh to .h5 format")
    # clean up directory
    os.remove(filename + ".xml")
    os.remove(filename + "_physical_region.xml")
    os.remove(filename + "_facet_region.xml")
    os.remove(filename + ".msh2")
    os.remove(filename + ".med")
    if logger is None:
        print("removed previous binary and text meshfiles that are not needed any longer")
    else:
        logger.debug("removed previous binary and text meshfiles that are not needed any longer")
    return
