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
import os
import errno
import numpy as np
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt

from emstimtools.utils.loadyml import (return_study, return_solver)
from emstimtools.utils.parseequations import parseEquations
from .Mesh import Mesh
from .Boundary import Boundary


class Fenics(object):
    """
    This class implements the main interface to FEniCS.
    The attributes :attr:`data` and :attr:`logger` are taken from main simulation class.

    Currently, the following modules are avaibable:

    - :class:`.ES` : ElectroStatics
    - :class:`.EQS` : ElectroQuasiStatics
    - :class:`.JouleHeating` : Joule Heating due to Ohmic heat

    We need to have (use yaml inputfile for that):

        - :attr:`mesh` : to be imported from h5 file OR from other object (e.g. in multiphysics)
        - :attr:`element` : element-type
        - :attr:`degree` : degree of polynomial
        - :attr:`solver` : solver type

    For the FEA certain options are set:

        - :attr:`ConstantsSpace` : a FunctionSpace with piecewise-constant Discontinuous Galerkin elements
        - :attr:`cells_array` : array for cells that can be used to save material constants

    For further information on which elements can be used, have a look at the
    `FEniCS website <https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/functions/functionspace/FunctionSpace.html>`_

    For visualization:
    Starting from Fenics 2017.2 and Paraview 5.5, XDMF files can be used for visualization AND checkpointing. See for instance slides `here <https://orbilu.uni.lu/bitstream/10993/35848/1/main.pdf>`_

    To tune the solver, choose something from the following keys.
    The following keys are accessible by just using the keyword, e.g. 'linear_solver' in the YAML-file:

        +---------------------------+---------+---------+
        | linear_variational_solver |  type   |  value  |
        +===========================+=========+=========+
        | linear_solver             |  string | default |
        +---------------------------+---------+---------+
        | preconditioner            |  string | default |
        +---------------------------+---------+---------+
        | print_matrix              |    bool |    0    |
        +---------------------------+---------+---------+
        | print_rhs                 |    bool |    0    |
        +---------------------------+---------+---------+
        | symmetric                 |    bool |    0    |
        +---------------------------+---------+---------+

    The following keys are accessible by using a nested dictionary starting with 'krylov_solver' in the YAML-file (only for Krylov subspace solvers):

        +--------------------------+----------+---------+
        | **krylov_solver**        |   type   |  value  |
        +==========================+==========+=========+
        | absolute_tolerance       |  double  | <unset> |
        +--------------------------+----------+---------+
        | divergence_limit         |  double  | <unset> |
        +--------------------------+----------+---------+
        | error_on_nonconvergence  |    bool  | <unset> |
        +--------------------------+----------+---------+
        | maximum_iterations       |     int  | <unset> |
        +--------------------------+----------+---------+
        | monitor_convergence      |    bool  | <unset> |
        +--------------------------+----------+---------+
        | nonzero_initial_guess    |    bool  | <unset> |
        +--------------------------+----------+---------+
        | relative_tolerance       |  double  | <unset> |
        +--------------------------+----------+---------+
        | report                   |    bool  | <unset> |
        +--------------------------+----------+---------+

    The following keys are accessible by using a nested dictionary starting with 'lu_solver' in the YAML-file (only for direct solvers):

        +----------------+-------+-------+
        | **lu_solver**  |  type | value |
        +================+=======+=======+
        | report         |  bool |   1   |
        +----------------+-------+-------+
        | symmetric      |  bool |   0   |
        +----------------+-------+-------+
        | verbose        |  bool |   0   |
        +----------------+-------+-------+


    .. todo:: make output nicer
    .. todo:: implement more properties, e.g. `capacitance or resistance <http://www.iue.tuwien.ac.at/phd/heinzl/node55.html>`_
    .. todo:: analyse system matrix, idea at bottom of file in comment
    """
    def __init__(self, data, logger):
        self.data = data
        self.study = return_study(self.data)
        self.logger = logger
        self.array_tag = False
        self.array_iter = None
        # needed, otherwise refinement does not work
        d.parameters["refinement_algorithm"] = "plaza_with_parent_facets"
        self.plot = None
        if 'plot' in self.data:
            self.plot = self.data['plot']
        self.load_mesh()
        if 'output' in self.data:
            self._prepare_result_dir()
        self._initalize_point_evaluation()
        self._prepare_fem()

    def load_mesh(self):
        """
        load mesh including subdomains and facets.
        if it is defined, plot the mesh.

        .. todo:: write check if subdomains and facets match mesh.

        """
        subdomains = None
        facets = None
        refinement = None
        if 'subdomains' in self.data:
            subdomains = self.data['subdomains']
        if 'facets' in self.data:
            facets = self.data['facets']
        if 'mesh' in self.data:
            change = self.data['mesh']
        self.mesh = Mesh(self.data['study'], self.logger, subdomains, facets, change)
        self.dx = d.Measure('dx', domain=self.mesh.mesh)
        self.ds = d.Measure('ds', domain=self.mesh.mesh, subdomain_data=self.mesh.facets)
        # try to plot mesh
        try:
            if 'mesh' in self.plot:
                if self.plot['mesh'] is True:
                    self._plot_mesh()
        except TypeError:
            pass

    def _plot_mesh(self):
        if self.mesh.dimension < 3:
            self.logger.info("Plotting mesh")
            d.plot(self.mesh.mesh, title='mesh')
            plt.show()
            p = d.plot(self.mesh.cells, title="subdomains ")
            loc = plticker.MultipleLocator(base=1.0)
            # some hacks from the internet
            plt.colorbar(p, format='%d', ticks=loc, fraction=0.03, pad=0.04)
            plt.show()
            self.logger.debug("Mesh was plotted, subdomains and facets also match the mesh")
        else:
            self.logger.info("Plotting mesh does not make sense in 3D. Output is written to paraview readable files.")
            self.mesh_dir = "mesh/"
            if not os.path.exists(self.mesh_dir):
                try:
                    os.makedirs(self.mesh_dir)
                except OSError as exc:  # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise

            tmp = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.mesh_dir + self.mesh.meshname + str('_plot.xdmf'))
            tmp.write(self.mesh.mesh)
            tmp.close()
            tmp = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.mesh_dir + self.mesh.meshname + str('_subdomains') + str('.xdmf'))
            tmp.write(self.mesh.cells)
            tmp.close()
            tmp = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.mesh_dir + self.mesh.meshname + str('_facets') + str('.xdmf'))
            tmp.write(self.mesh.facets)
            tmp.close()

            self.logger.debug("Mesh was written out.")

    def _prepare_result_dir(self):
        self.result_dir = "results/"
        if not os.path.exists(self.result_dir):
            try:
                os.makedirs(self.result_dir)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

    def _prepare_fem(self):
        self.element = d.FiniteElement(self.data['element'], self.mesh.mesh.ufl_cell(), self.data['degree'])
        if 'properties' in self.data:
            if 'project_element' not in self.data['properties']:
                self.data['properties']['project_element'] = self.data['element']
            if 'project_degree' not in self.data['properties']:
                if self.data['degree'] > 1:
                    deg = self.data['degree'] - 1
                else:
                    deg = self.data['degree']
                self.data['properties']['project_degree'] = deg
            if 'project_solver' not in self.data['properties']:
                self.data['properties']['project_solver'] = 'mumps'
            if 'project_preconditioner' not in self.data['properties']:
                self.data['properties']['project_preconditioner'] = 'none'

                
        # to set constants element-wise, we choose a discontinuous basis of degree 0 (piecewise constant)
        self.ConstantsSpace = d.FunctionSpace(self.mesh.mesh, 'DG', 0)
        self.cells_array = np.asarray(self.mesh.cells.array(), dtype=np.int32)

    def _setup_solver(self):
        self.solver = d.LinearVariationalSolver(self.problem)
        solver_info = return_solver(self.data)
        if solver_info is not None:
            for i in solver_info:
                if isinstance(solver_info[i], dict):
                    for j in solver_info[i]:
                        self.logger.debug("Setting [" + str(i) + "][" + str(j) + "] as " + str(solver_info[i][j]))
                        self.solver.parameters[i][j] = solver_info[i][j]
                else:
                    self.logger.debug("Setting [" + str(i) + "] as " + str(solver_info[i]))
                    self.solver.parameters[i] = solver_info[i]

    def _set_material_constant(self, phys_constant):
        """
        set material constant. pass the name of the constant as a string.

        .. note:: if you want to add your own constant, pay attention that you set the vector for the parameter properly!

        .. todo:: fix ordering in salome or enable computation when SALOME subdomains are not strictly ordereded. Right now it seems that it works if facets are set before subdomains in SALOME.
        .. todo:: check if elif case for resistivity is needed any longer
        """

        # check if attribute exists
        # for array job: just initialise
        if hasattr(self, phys_constant) is False:
            setattr(self, phys_constant, d.Function(self.ConstantsSpace))
            if self.array_tag is True:
                return
        # allocate fenics function
        if phys_constant == 'resistivity':
            self.data['conductivity'] = {}
            for i in self.mesh.subdomaininfo:
                self.data['conductivity'][i] = 1. / self.data[phys_constant][i]
            setattr(self, 'conductivity', d.Function(self.ConstantsSpace))
            # phys_constant = 'conductivity'

        # make temparray to extract everything
        tmparray = np.empty(len(self.mesh.subdomaininfo))
        for i in self.mesh.subdomaininfo:
            # here float is necessary as YAML does not recognize scientific notation
            tmparray[int(self.mesh.subdomaininfo[i] - np.min(self.mesh.cells.array()))] = float(self.data[phys_constant][i])
        if phys_constant == 'conductivity':
            self.conductivity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        elif phys_constant == 'permittivity':
            self.permittivity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        # needed any longer????
        elif phys_constant == 'resistivity':
            # iterate through array and invert resistivity to get conductivity
            for x in np.nditer(tmparray, op_flags=['readwrite']):
                x[...] = 1. / x
            self.conductivity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        elif phys_constant == 'density':
            self.density.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        elif phys_constant == 'heat_capacity':
            self.heat_capacity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        elif phys_constant == 'thermal_conductivity':
            self.thermal_conductivity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)
        elif phys_constant == 'thermal_resistivity':
            self.thermal_resistivity.vector()[:] = np.choose(self.cells_array - np.min(self.mesh.cells.array()), tmparray)

    def _prepare_boundaries(self):
        self.boundaries = Boundary(self.data['boundaries'], self.mesh.facets, self.mesh.facetinfo, self.mesh.cells, self.mesh.subdomaininfo)

    def _project_on_submeshes(self):
        """
        this function provides utilities to write the solution on subdomains of the mesh to a paraview-readable file.
        Enhances the postprocessing process crucially.
        """
        self.logger.debug('Projecting solution on submeshes')

    def _postprocess(self):
        """
        Steps for postprocessing the results.
        Possibilities are:

        #. Projection on SubMeshes in :meth:`_project_on_submeshes`

        """
        self.logger.debug('Entering Postprocessing')
        if 'submesh' in self.data['postprocess']:
            self._project_on_submeshes()

    def _plot_submesh(self, subdomain):
        try:
            if 'mesh' in self.plot:
                if self.plot['mesh'] is True:
                        tmp = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.mesh_dir + self.mesh.meshname + '_plot_' + str(subdomain) + '.xdmf')
                        tmp.write(self.sub_mesh)
                        tmp.close()
        except TypeError:
            pass

    def _set_function_space(self):
        self.V = d.FunctionSpace(self.mesh.mesh, self.element)
        self.logger.info("Number of DOFs : {}".format(self.V.dim()))

    def _set_frequencies(self):
        if isinstance(self.data['frequencies'], float) or isinstance(self.data['frequencies'], list):
            self.frequencies = np.array(self.data['frequencies'], ndmin=1)
        else:
            self.frequencies = np.array(eval(self.data['frequencies']), ndmin=1)
        self.logger.debug("Will compute for the frequencies: {}".format(self.frequencies))

    def _set_timestepping(self):
        if 'timesteps' in self.data:
            timeinfo = self.data['timesteps']
            try:
                self.t_max = timeinfo['t_max']
            except KeyError:
                self.logger.warning("You have to specify the total simulation time!")
            try:
                self.dt = timeinfo['dt']
            except KeyError:
                self.logger.warning("You have to specify the time step!")
            assert self.t_max > self.dt, "time step is larger than total simulation time!"
            self.N = int(self.t_max / self.dt)
            self.t = 0.
        else:
            raise Exception("You must enter information about timesteps using key \" timesteps \"")

        self.logger.debug("Set timestep to {} and will run for total time {}, i.e. {} steps.".format(self.dt, self.t_max, self.N))

    def _initalize_point_evaluation(self):
        self.evaluate_points_solution = None
        try:
            self.evaluate_points_solution = self.data['evaluate_solution']
        except (KeyError, ValueError):
            pass
        # parse equations to values
        if self.evaluate_points_solution is not None:
            self.evaluate_points_solution = parseEquations(self.evaluate_points_solution, self.data['geometryvalues'], self.data['dimension'])
            if self.array_tag is False:
                filestring = self.result_dir + 'solution_' + self.study + '.dat'
            else:
                filestring = self.result_dir + 'solution_' + self.study + '_' + self.array_iter + '.dat'
            self.datafileTXT = open(filestring, 'w')


"""
from dolfin import *
import ufl
import scipy.sparse as sp
import matplotlib.pyplot as plt


def spy_form_matrix(A, l=None, bcs=None, marker_size=2.0, show=True, pattern_only=True, **kwargs):
    assert isinstance(A, ufl.form.Form)
    if l: assert isinstance(l, ufl.form.Form)
    eigen_matrix = EigenMatrix()
    if not l:
        assemble(A, tensor=eigen_matrix)
        if not bcs is None:
            for bc in bcs: bc.apply(eigen_matrix)
    else:
        if not bcs is None:
            SystemAssembler(A, l, bcs).assemble(eigen_matrix)
        else:
            SystemAssembler(A, l).assemble(eigen_matrix)
    A = eigen_matrix

    row, col, data = A.data()
    if pattern_only:
        data[:] = 1.0

    sp_mat = sp.csr_matrix((data, col, row), dtype='float')
    plt.spy(sp_mat, markersize=marker_size, precision=0, **kwargs)
    print(sp_mat)
    print(sp.tril(sp_mat))
    print("\n")
    print(sp.triu(sp_mat))
    if show: plt.show()


mesh = UnitSquareMesh(4, 4)
V = FunctionSpace(mesh, 'CG', 1)
u, v = TrialFunction(V), TestFunction(V)

a = dot(grad(u), grad(v))*dx
L = Constant(1.0)*v*dx
spy_form_matrix(a, L)
"""
