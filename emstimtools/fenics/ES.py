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

"""
compute EM fields in stimulation tools by solving the ES (Laplace or Poisson) equation.

"""
import dolfin as d
import numpy as np
from .fenics import Fenics


def ESLHS(conductivity, u, v, dx):
    a = d.inner(conductivity * d.grad(u), d.grad(v)) * dx
    return a


def ESRHS(source_term, v, dx):
    L = -source_term * v * dx
    return L


class ES(Fenics):
    r"""
    Define a class, where an entire ES run will be done.
    The following PDE is solved:

    .. math::

         \nabla\cdot\left[\sigma(\vec{r},\omega)\nabla\Phi(\vec{r})\right]=0 \enspace,

    where :math:`\sigma` is the conductivity that is defined for each subdomain.

    We need the following attributes:

        - :attr:`conductivity`  : conductivities
        - :attr:`V` : FunctionSpace, e.g. CG 2

    .. todo:: add source term in case of charge sources, i.e. Poisson equation
    .. todo:: implement array jobs
    """
    def __init__(self, data, logger):
        # use __init__ from Fenics class
        super().__init__(data, logger)
        ###
        # Physics part
        ###
        if self.data['conductivity'] is not None:
            self._set_material_constant('conductivity')
        else:
            raise Exception("Conductivity must be known!")
        self._prepare_output()
        ###
        # FEM part
        ###
        self._set_function_space()
        self._prepare_boundaries()
        self.boundaries.set_DirichletBC(self.V, self.logger)
        self._set_problem()
        self._setup_solver()
        self.logger.info("Starting ES simulation.")
        self._solve_problem()
        self._close_output()
        self.logger.debug('checking for postprocessing steps')
        if 'postprocess' in self.data:
            self._postprocess()

    def _set_problem(self):
        # for now hard-coded
        self.source_term = d.Constant(0.0)
        # trial function
        u = d.TrialFunction(self.V)
        # test function
        v = d.TestFunction(self.V)
        a = ESLHS(self.conductivity, u, v, self.dx)
        L = ESRHS(self.source_term, v, self.dx)
        a, L = self.boundaries.set_RobinBCPoisson(a, L, u, v, self.ds, 'R', 'V_ref', self.logger)
        # vector for solution
        self.u = d.Function(self.V)
        self.u.rename('potential', 'potential')
        self.problem = d.LinearVariationalProblem(a, L, self.u, self.boundaries.bc)

    def _solve_problem(self):
        self.solver.solve()
        self._write_results()

    def _write_results(self):
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5.write(self.u, "/potential")
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self.logger.debug("Write to xdmf file")
                    self.datafileXDMF.write(self.u)
        if 'properties' in self.data:
            try:
                if self.data['properties']['E-Field'] is True:
                    self.get_field()
                    self.datafileXDMFE.write(self.Efield)
                if 'Current' in self.data['properties']:
                    self.get_current()
                if 'Impedance' in self.data['properties']:
                    self.get_impedance()
            except KeyError:
                pass

    def get_field(self):
        """
        compute field
        """
        self.Vector = d.VectorFunctionSpace(self.mesh.mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        self.Efield = d.project(-d.grad(self.u), self.Vector, solver_type=self.data['properties']['project_solver'], preconditioner_type=self.data['properties']['project_preconditioner'])
        self.Efield.rename('E-Field', 'E-Field')

    def _prepare_output(self):
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.logger.info("Will store solution to hdf5 format")
                    self.datafileHDF5 = d.HDF5File(self.mesh.mesh.mpi_comm(), self.result_dir + "solution" + self.study + ".h5", "w")

            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self._open_xdmf_output_solution()
        if 'properties' in self.data:
            try:
                if self.data['properties']['E-Field'] is True:
                    self._open_output_efield()
            except KeyError:
                pass

    def _open_xdmf_output_solution(self, add=''):
        self.logger.info("Will store solution to xdmf format")
        filename = self.result_dir + "solution" + self.study + add + '.xdmf'
        self.logger.info("Will store solution to file {}".format(filename))
        self.datafileXDMF = d.XDMFFile(self.mesh.mesh.mpi_comm(), filename)

    def _open_xdmf_output_difference(self, add=''):
        filename = self.result_dir + "solution_difference" + self.study + add + '.xdmf'
        self.logger.info("Will store solution difference to file {}".format(filename))
        self.datafileXDMFdiff = d.XDMFFile(self.mesh.mesh.mpi_comm(), filename)

    def _open_output_efield(self, add=''):
        filename = self.result_dir + "field" + self.study + add + '.xdmf'
        self.logger.info("Will store e-field to file {}".format(filename))
        self.datafileXDMFE = d.XDMFFile(self.mesh.mesh.mpi_comm(), filename)

    def _close_solution_xdmf(self):
        self.logger.debug("Closed XDMF file for solution")
        self.datafileXDMF.close()

    def _close_solution_difference_xdmf(self):
        self.logger.debug("Closed XDMF file for solution difference")
        self.datafileXDMFdiff.close()

    def _close_efield_xdmf(self):
        self.logger.debug("Closed XDMF file for e-field")
        self.datafileXDMFE.close()

    def _close_output(self):
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5.close()

            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self._close_solution_xdmf()
        if 'properties' in self.data:
            try:
                if self.data['properties']['E-Field'] is True:
                    self._close_efield_xdmf()
            except KeyError:
                pass

    def _write_projected_solution(self, i):
        """
        take domain :param str i: and project solution there
        """
        self.sub_mesh = d.SubMesh(self.mesh.mesh, self.mesh.cells, self.mesh.subdomaininfo[i])
        V_sub = d.FunctionSpace(self.sub_mesh, self.element)
        # real part, mumps to avoid memory overflow
        u_sub = d.project(self.u, V_sub, solver_type=self.data['properties']['project_solver'], preconditioner_type=self.data['properties']['project_preconditioner'])
        u_sub.rename('potential', 'potential')
        self.datafileXDMF.write(u_sub)

    def _write_projected_efield(self, i):
        """
        take domain :param str i: and project e-field there
        """
        Vector_sub = d.VectorFunctionSpace(self.sub_mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        efield_sub = d.project(self.Efield, Vector_sub, solver_type=self.data['properties']['project_solver'], preconditioner_type=self.data['properties']['project_preconditioner'])
        efield_sub.rename('E-Field', 'E-Field')
        self.datafileXDMFE.write(efield_sub)
        self.logger.debug("Wrote e-field for domain " + i)

    def _project_on_submeshes(self):
        """
        iterate through all domains that are provided in input file and write solution and if needed field.

        .. todo:: support for array jobs
        """
        self.logger.debug('Projecting ES solutions on submeshes')
        domains = self.data['postprocess']['submesh']  # a list of strings
        for i in domains:
            self._open_xdmf_output_solution(add='_' + i)
            self._write_projected_solution(i)
            self._close_solution_xdmf()
            if 'properties' in self.data:
                try:
                    if self.data['properties']['E-Field'] is True:
                        self._open_output_efield(add='_' + i)
                        self._write_projected_efield(i)
                        self._close_efield_xdmf()
                except KeyError:
                    pass

    def get_current(self):
        self.currents = {}
        self.logger.debug('Computing the electric currents through certain surfaces.')
        self.normal = d.FacetNormal(self.mesh.mesh)
        if not isinstance(self.data['properties']['Current'], list):
            self.data['properties']['Current'] = [self.data['properties']['Current']]
        for facet in self.data['properties']['Current']:
            current = self.compute_current(facet)
            self.currents[facet] = current
        self.logger.info('The currents are:')
        self.logger.info(self.currents)

    def compute_current(self, facet):
        if not hasattr(self, "normal"):
            self.normal = d.FacetNormal(self.mesh.mesh)
        if facet not in self.mesh.facetinfo:
            raise Exception('The facet {} is not in this geometry. Choose one of the following facets {}'.format(facet, self.mesh.facetinfo.keys))
        index = self.boundaries.facet_dict[facet]
        current = d.dot(-self.conductivity * d.grad(self.u), self.normal) * self.ds(index)
        return d.assemble(current)

    def get_impedance(self):
        """
        .. todo:: implement proper impedance on surface
        """
        self.logger.debug('Computing the impedance on a surface.')
        if not isinstance(self.data['properties']['Impedance'], list):
            self.data['properties']['Impedance'] = [self.data['properties']['Impedance']]
        self.impedances = {}
        for facet in self.data['properties']['Impedance']:
            if facet not in self.mesh.facetinfo:
                raise Exception('The facet {} is not in this geometry. Choose one of the following facets {}'.format(facet, self.mesh.facetinfo.keys))
            # get voltage:
            try:
                voltage = self.data['boundaries']['Dirichlet'][facet]
            except KeyError:
                self.logger.error("Currently, there is no possibility to compute the impedance on boundaries that do not have a fixed voltage.")
            if hasattr(self, 'currents'):
                if facet in self.currents:
                    self.impedances[facet] = voltage / self.currents[facet]
                    continue
            else:
                self.impedances[facet] = voltage / self.compute_current(facet)
        self.logger.info('The impedances are:')
        self.logger.info(self.impedances)

    def _interpolate_external_solution(self, ext_solution):
        """
        Takes an externally generated solution (e.g. on a coarser mesh) and interpolates it onto the finer mesh.
        Requires :param:`ext_solution` to be a DOLFIN Expression or a DOLFIN function from a Lagrange space.
        """

        self.u_ext = d.Function(self.V)
        # needed since otherwise errors are introduced on the boundary
        self.u_ext.set_allow_extrapolation(True)
        ext_solution.set_allow_extrapolation(True)

        d.LagrangeInterpolator.interpolate(self.u_ext, ext_solution)

    def compare_solution(self, ext_solution, save_difference=False):
        """
        compare solution to other FE-generated solution
        """
        self._interpolate_external_solution(ext_solution)
        self.solution_difference = d.Function(self.V)
        self.solution_difference.rename('solution_difference', 'solution_difference')
        self.solution_difference.vector()[:] = self.u.vector()[:]
        self.solution_difference.vector()[:] -= self.u_ext.vector()[:]
        if save_difference:
            if not hasattr(self, "result_dir"):
                self._prepare_result_dir()
            self._open_xdmf_output_difference()
            self.datafileXDMFdiff.write(self.solution_difference)
            self._close_solution_difference_xdmf()

    def get_solution_error(self, ext_solution, save_difference=False, norm_type="l2"):
        """
        get the error in a certain norm. possible are l1, l2 and linf norms.
        Note that here only the norm of the difference vector is taken!
        """
        self.compare_solution(ext_solution, save_difference)
        return self.solution_difference.vector().norm(norm_type)

    def get_solution_l2error(self, ext_solution, save_difference=False):
        self.compare_solution(ext_solution, save_difference)
        return np.sqrt(d.assemble(d.inner(self.u - self.u_ext, self.u - self.u_ext) * self.dx))

