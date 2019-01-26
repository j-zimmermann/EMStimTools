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

import dolfin as d
import numpy as np


class Mesh(object):
    '''
    A Mesh instance has the following attributes:

    - :attr:`meshname` : title of the study and '_mesh'. Also indicates the name of the meshfile (.h5).
    - :attr:`mesh` : the mesh data
    - :attr:`cells` : the subdomains. A fenics function.
    - :attr:`facets` : the facets (needed for BCs). A fenics function.
    - :attr:`dimension` : the geometrical dimension.
    - :attr:`subdomaininfo` : dict with info on subdomains. similarly :attr:`facetinfo` exists.
    - :attr:`refinement_info` : info about whether and which regions to refine. refinement is carried out when initialising mesh.

    A mesh instance has the following method:

    - :meth:`refine`: uniform refinement of marked regions.
    '''
    def __init__(self, title, logger, subdomaininfo, facetinfo, refinement_info):
        self.mesh_dir = "mesh/"
        self.meshname = str(title) + '_mesh'
        self.mesh_file = self.mesh_dir + self.meshname + '.h5'
        self.mesh = d.Mesh()
        self.logger = logger
        self.logger.info("Loading mesh")
        hdf = d.HDF5File(self.mesh.mpi_comm(), self.mesh_file, "r")
        hdf.read(self.mesh, "/mesh", False)
        self.cells = d.MeshFunction('size_t', self.mesh, self.mesh.geometric_dimension())
        self.facets = d.MeshFunction('size_t', self.mesh, self.mesh.geometric_dimension() - 1)
        hdf.read(self.cells, "/subdomains")
        hdf.read(self.facets, "/facets")
        self.dimension = self.mesh.topology().dim()
        self.subdomaininfo = subdomaininfo
        self.facetinfo = facetinfo
        self.refinement_info = refinement_info
        self.logger.info("Loaded mesh")
        if self.refinement_info is not None:
            self.refine()

    def refine(self):
        for n in range(self.refinement_info['max_iter']):
            cell_markers = d.MeshFunction('bool', self.mesh, self.dimension)
            cell_markers.set_all(False)
            hlp = np.asarray(self.cells.array(), dtype=np.int32)
            for i in range(hlp.size):
                if any(hlp[i] == self.subdomaininfo[a] for a in self.refinement_info['regions']):
                    cell_markers[i] = True

            self.mesh = d.refine(self.mesh, cell_markers)
            self.facets = d.adapt(self.facets, self.mesh)
            self.cells = d.adapt(self.cells, self.mesh)
        self.logger.info("refined mesh {} times at {}".format(self.refinement_info['max_iter'], self.refinement_info['regions']))
        try:
            meshsave = self.refinement_info['save_mesh']
            self.logger.info('Will save mesh.')
            if meshsave is True:
                hdfout = d.HDF5File(self.mesh.mpi_comm(), self.meshname + "_refined.h5", "w")
                hdfout.write(self.mesh, "/mesh")
                hdfout.write(self.cells, "/subdomains")
                hdfout.write(self.facets, "/facets")
            else:
                self.logger.info("Mesh won't be saved")
        except KeyError:
            self.logger.info("Mesh won't be saved")
