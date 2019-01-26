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

"""Simulate Joule (Resistive) Heating in FEniCS"""
import dolfin as d
from emstimtools.fenics.fenics import Fenics
from emstimtools.fenics.Boundary import Boundary
from emstimtools.utils.loadyml import return_solver


class JouleHeating(Fenics):
    """
    Solve a Joule heating problem. Currently only for electrostatics.
    Idea:

    #. initializes an object with a mesh and all needed stuff. then two problems are defined and solved.

    One can use a temperature-dependent resistivity as in COMSOL.

    .. todo:: implement dielectric heating, i.e. when there is time-harmonic signal
    .. todo:: write function for write out to files, avoid Key Errors
    """
    def __init__(self, data, logger):
        # use __init__ from Fenics class
        super(JouleHeating, self).__init__(data, logger)
        self.logger.info("Starting Joule Heating simulation.")
        ###
        # Physics part
        ###
        self._set_timestepping()
        self._set_material_constants()
        ###
        # FEM part
        ###
        self._set_function_space()
        self._prepare_boundaries()
        self.logger.debug("Set DirichletBC")
        self.boundariesEM.set_DirichletBC(self.V, self.logger)
        self.boundariesT.set_DirichletBC(self.V, self.logger)
        self._set_problem()
        self.logger.debug("Problem set.")
        self._setup_solver()
        self._solve_problem()

    def _prepare_boundaries(self):
        """ need to redefine function as we have two kinds of boundaries """
        self.boundariesEM = Boundary(self.data['boundaries']['electric'], self.mesh.facets, self.mesh.facetinfo)
        self.boundariesT = Boundary(self.data['boundaries']['thermal'], self.mesh.facets, self.mesh.facetinfo)

    def _set_problem(self):
        self.logger.debug("Set problem")
        # for now hard-coded
        self.source_term = d.Constant(0.0)
        self.logger.debug("trial functions")
        self.psi = d.TrialFunction(self.V)
        # Define initial value, which will later be the solution from the previous (n-th) time step
        self.T_n = d.Function(self.V)
        self.T_n = d.interpolate(self.T_ref, self.V)

        self.logger.debug("test functions")
        v = d.TestFunction(self.V)
        dx = d.dx
        self.logger.debug("vector for solution")
        self.u = d.Function(self.V)
        a = d.inner(self._sigma_T(self.T_n) * d.grad(self.psi), d.grad(v)) * dx
        L = self.source_term * v * dx
        self.psi = d.Function(self.V)  # solution for potential
        self.EMproblem = d.LinearVariationalProblem(a, L, self.psi, self.boundariesEM.bc)
        self.logger.debug("set EM problem")
        # temperature
        vT = d.TestFunction(self.V)
        T_n1 = d.TrialFunction(self.V)

        a_1 = (self.density * self.heat_capacity * T_n1 * vT * dx + self.thermal_conductivity * self.dt * d.dot(d.grad(T_n1), d.grad(vT)) * dx)
        L_1 = (self.density * self.heat_capacity * self.T_n + self.dt * self._sigma_T(self.T_n) * d.dot(d.grad(self.psi), d.grad(self.psi))) * vT * dx
        self.T_n1 = d.Function(self.V)  # solution for T
        self.Tproblem = d.LinearVariationalProblem(a_1, L_1, self.T_n1, self.boundariesT.bc)
        self.logger.debug("set thermal problem")

    def _setup_solver(self):
        self.EMsolver = d.LinearVariationalSolver(self.EMproblem)
        self.Tsolver = d.LinearVariationalSolver(self.Tproblem)
        solver_info = return_solver(self.data)
        if solver_info is not None:
            for i in solver_info:
                self.logger.debug("Setting " + str(self.EMsolver.parameters[i]) + " as " + str(solver_info[i]))
                self.EMsolver.parameters[i] = solver_info[i]
                self.Tsolver.parameters[i] = solver_info[i]

    def _solve_problem(self):
        append = False
        for n in range(self.N):
            if 'thermal_resistivity' in self.data:
                # solve EM problem for every time step
                self.EMsolver.solve()
                if self.data['output']['XDMF'] is True:
                    self.datafileXDMFEM.write_checkpoint(self.psi, 'potential', self.t, append=append)
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5EM.write(self.psi, "/V_{}".format(n))

            else:
                # solve once and forever, only solution for temperature for each time step
                if n == 0:
                    self.EMsolver.solve()
                    if self.data['output']['XDMF'] is True:
                        self.datafileXDMFEM.write_checkpoint(self.psi, 'potential', self.t, append=append)
                    if self.data['output']['HDF5'] is True:
                        self.datafileHDF5EM.write(self.psi, "/V_{}".format(n))

            # Compute solution
            self.Tsolver.solve()
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
        if self.data['output']['HDF5'] is True:
            self.datafileHDF5T.write(self.T_n, "/T_{}".format(self.N))
            if 'thermal_resistivity' in self.data:
                self.datafileHDF5EM.write(self.psi, "/V_{}".format(self.N))
        if self.data['output']['XDMF'] is True:
            self.datafileXDMFT.write_checkpoint(self.T_n, 'temperature', self.t, append=append)
            if 'thermal_resistivity' in self.data:
                self.datafileXDMFEM.write_checkpoint(self.psi, 'potential', self.t, append=append)

    def _set_material_constants(self):
        # global definitions
        try:
            temp_info = self.data['temperature']
            self.T_0 = d.Constant(temp_info['initial'])
            self.T_ref = d.Constant(temp_info['reference'])
        except KeyError:
            self.logger.warning("You must provide the temperature in a dictionary containing \" initial \" and \" reference \" as keys.")
        # material-specific definitions
        if 'resistivity' in self.data:
            self._set_material_constant('resistivity')
        elif 'conductivity' in self.data:
            self._set_material_constant('conductivity')
        else:
            raise Exception("Conductivity or resistivity must be known!")
        if 'heat_capacity' in self.data:
            self._set_material_constant('heat_capacity')
        else:
            raise Exception("Heat capacity must be known!")
        if 'density' in self.data:
            self._set_material_constant('density')
        else:
            self.logger.warning("density must be known!")
        if 'thermal_conductivity' in self.data:
            self._set_material_constant('thermal_conductivity')
        else:
            self.logger.warning("thermal conductivity must be known!")

        if 'thermal_resistivity' in self.data:
            self._set_material_constant('thermal_resistivity')

        # temperature-dependent resistivity
    def _sigma_T(self, Temp):
        if 'thermal_resistivity' in self.data:
            return self.conductivity / (1. + self.thermal_resistivity * (Temp - self.T_0))
        else:
            return self.conductivity

    def compute_flux(self):
        self.hflux = d.project(-self.thermal_conductivity * d.grad(self.T_n), d.VectorFunctionSpace(self.mesh.mesh, 'Lagrange', 1))
