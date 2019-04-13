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
from emstimtools.fenics.fenics import Fenics
from emstimtools.fenics.Boundary import Boundary
from emstimtools.utils.loadyml import return_solver


def HeatLHSstat(T, v, thermal_conductivity, dx):
    """
    LHS for stationary case
    """
    a = d.inner(thermal_conductivity * d.grad(T), d.grad(v)) * dx
    return a


def HeatLHSstattime(T, v, density, heat_capacity, thermal_conductivity, dx, dt):
    """
    LHS for time-dependent case

    .. todo:: check

    """
    a = (density * heat_capacity * T * v * dx + thermal_conductivity * dt * d.dot(d.grad(T), d.grad(v)) * dx)
    return a


def HeatRHSstat(v, dx):
    """
    RHS for stationary case
    """
    L = d.Constant(0) * v * dx
    return L


def HeatRHStime(T, v, density, heat_capacity, thermal_conductivity, dx, dt):
    """
    RHS for time-dependent case

    .. todo:: check

    """
    L = density * heat_capacity * T
    return L


class Heat(Fenics):
    """
    Solve heat equation.
    """
    def __init__(self, data, logger):
        # use __init__ from Fenics class
        super().__init__(data, logger)
        self.logger.info("Starting Heat simulation.")
        ###
        # FEM part
        ###
        self._set_function_space()
        self._prepare_boundaries()
        self.logger.debug("Set DirichletBC")
        self.boundaries.set_DirichletBC(self.V, self.logger)
        self.datafileXDMFT = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "Tsolution" + self.    study + '.xdmf')
        ###
        # Physics part
        ###
        stationary = 'timesteps' not in self.data
        if stationary:
            self._set_material_constants(stat=True)
            self._set_stationary_problem()
            self.logger.debug("Stationary problem set.")
            self.N = 1
            self.t = 0
        else:
            self._set_timestepping()
            self._set_material_constants()
            self._set_problem()
            self.logger.debug("Problem set.")
        self._setup_solver()
        self._solve_problem()
        self.datafileXDMFT.close()

    '''
    def _prepare_output(self):
        """
        write to HDF5 or XDMF file if specified
        if solution at points is evaluated, also write to txt files
        """
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.logger.info("Will store solution to hdf5 format")
                    filestring = self.result_dir + "solution" + self.study + ".h5"
                    self.datafileHDF5 = d.HDF5File(self.mesh.mesh.mpi_comm(), filestring, "w")
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self.logger.info("Will store solution to xdmf format")
                    filestring = self.result_dir + "solution" + self.study + ".xdmf"
                    self.datafileXDMFT = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)
    '''

    def _prepare_boundaries(self):
        self.boundaries = Boundary(self.data['boundaries'], self.mesh.facets, self.mesh.facetinfo, self.mesh.cells, self.mesh.subdomaininfo)

    def _set_stationary_problem(self):
        T = d.TrialFunction(self.V)
        v = d.TestFunction(self.V)
        a = HeatLHSstat(T, v, self.thermal_conductivity, self.dx)
        L = HeatRHSstat(v, self.dx)
        a, L = self.boundaries.set_RobinBCPoisson(a, L, T, v, self.ds, 'h', 'T_ref', self.logger)

        self.T_n1 = d.Function(self.V)  # solution for T
        self.Tproblem = d.LinearVariationalProblem(a, L, self.T_n1, self.boundaries.bc)

    def _set_problem(self):
        """
        .. todo:: needs implementation
        """
        self.logger.debug("Set problem")
        # for now hard-coded
        self.T_n = d.Function(self.V)
        self.T_n.rename('temperature', 'temperature')
        self.T_n = d.interpolate(self.T_ref, self.V)

        self.logger.debug("test functions")
        v = d.TestFunction(self.V)
        # temperature
        v = d.TestFunction(self.V)
        T_n1 = d.TrialFunction(self.V)

        self.T_n1 = d.Function(self.V)  # solution for T
        self.Tproblem = d.LinearVariationalProblem(a_1, L_1, self.T_n1, self.boundariesT.bc)
        self.logger.debug("set thermal problem")

    def _setup_solver(self):
        self.Tsolver = d.LinearVariationalSolver(self.Tproblem)
        solver_info = return_solver(self.data)
        if solver_info is not None:
            for i in solver_info:
                self.logger.debug("Setting " + str(self.Tsolver.parameters[i]) + " as " + str(solver_info[i]))
                self.Tsolver.parameters[i] = solver_info[i]

    def _solve_problem(self):
        append = False
        for n in range(self.N):
            # Compute solution
            self.Tsolver.solve()
            # if stationary: leave here
            if self.N == 1:
                break

            # write previous solution to file
            if self.data['output']['HDF5'] is True:
                self.datafileHDF5T.write(self.T_n, "/T_{}".format(n))
            if self.data['output']['XDMF'] is True:
                self.datafileXDMFT.write_checkpoint(self.T_n, 'temperature', self.t, append=append)
            # Update previous solution
            self.T_n.assign(self.T_n1)

            # Update current time
            self.t += self.dt
            # to write time series
            append = True
        try:
            if self.data['output']['HDF5'] is True:
                self.datafileHDF5T.write(self.T_n1, "/T_{}".format(self.N))
        except KeyError:
            pass
        try:
            if self.data['output']['XDMF'] is True:
                self.datafileXDMFT.write_checkpoint(self.T_n1, 'temperature', self.t, append=append)
        except KeyError:
            pass

    def _set_material_constants(self, stat=False):
        # first what is needed for stationary case
        if 'thermal_conductivity' in self.data:
            self._set_material_constant('thermal_conductivity')
        else:
            self.logger.warning("thermal conductivity must be known!")
        if stat:
            return

        # global definitions for time-dependent case
        try:
            temp_info = self.data['temperature']
            self.T_0 = d.Constant(temp_info['initial'])
            self.T_ref = d.Constant(temp_info['reference'])
        except KeyError:
            self.logger.warning("You must provide the temperature in a dictionary containing \" initial \" and \" reference \" as keys.")
        # material-specific definitions
        if 'heat_capacity' in self.data:
            self._set_material_constant('heat_capacity')
        else:
            raise Exception("Heat capacity must be known!")
        if 'density' in self.data:
            self._set_material_constant('density')
        else:
            self.logger.warning("density must be known!")

    def compute_flux(self):
        self.hflux = d.project(-self.thermal_conductivity * d.grad(self.T_n), d.VectorFunctionSpace(self.mesh.mesh, self.data['properties']['project_element'], self.data['properties']['project_degree']))
