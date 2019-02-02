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
from .fenics import Fenics


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
        dx = d.dx
        a = d.inner(self.conductivity * d.grad(u), d.grad(v)) * dx
        L = -self.source_term * v * dx
        # vector for solution
        self.u = d.Function(self.V)
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
            except KeyError:
                pass

    def get_field(self):
        """
        compute field
        """
        self.Vector = d.VectorFunctionSpace(self.mesh.mesh, self.data['element'], self.data['degree'] - 1)
        self.Efield = d.FunctionSpace(self.Vector)
        self.Efield = d.project(-d.grad(self.u), self.Vector, solver_type='mumps')

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
        self.logger.info("Will store solution to file {}".format(self.result_dir + "solution" + self.study + add + '.xdmf'))
        self.datafileXDMF = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "solution" + self.study + add + '.xdmf')

    def _open_output_efield(self, add=''):
        self.logger.info("Will store e-field to file {}".format(self.result_dir + "field" + self.study + add + '.xdmf'))
        self.datafileXDMFE = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "field" + self.study + add + '.xdmf')

    def _close_solution_xdmf(self):
        self.logger.debug("Closed XDMF file for solution")
        self.datafileXDMF.close()

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
        u_sub = d.project(self.u, V_sub, solver_type='mumps')
        self.datafileXDMF.write(u_sub)

    def _write_projected_efield(self, i):
        """
        take domain :param str i: and project e-field there
        """
        Vector_sub = d.VectorFunctionSpace(self.sub_mesh, self.data['element'], self.data['degree'] - 1)
        efield_sub = d.project(self.Efield, Vector_sub, solver_type='mumps')
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
