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
import numpy as np
from emstimtools.fenics.fenics import Fenics
from emstimtools.fenics.Boundary import Boundary
from emstimtools.utils.loadyml import return_solver
from emstimtools.fenics.ES import ESLHS, ESRHS
from emstimtools.fenics.EQS import EQSLHS, EQSRHS
from emstimtools.fenics.HeatEquation import HeatLHStime, HeatLHSstat, HeatRHStime, HeatRHSstat


class JouleHeating(Fenics):
    """
    Solve a Joule heating problem.

    The time-dependent or stationary heat equation (think of :math:`\frac{\partial T}{\partial t} \to 0`)  is solved:
    .. math::

        \rho C_p \frac{\partial T}{\partial t} =  \nabla \cdot (k_{\mathrm{iso}} \nabla T) + Q_{\mathrm{e}}

    The heat :math:`Q_{\mathrm{e}}` due to the current is either:
    .. math::

        Q_\mathrm{e} = \sigma |\nabla V|^2

    or in the frequency-dependent case (:math:`V=V_r + i V_i` is a phasor)
    .. math::

        Q_\mathrm{e} =  \frac{1}{2} \mathrm{Re}\left[(\sigma + i \omega \epsilon) (\nabla V)^\ast \nabla V \right] = \frac{1}{2} \sigma (|\nabla V_r|^2 i+ |\nabla V_i|^2)

    One can use a temperature-dependent resistivity as in COMSOL.
    .. note ::

        Technically, this is not completely solved correctly. A better formulation for a truly implicit scheme is needed. Here, the temperature of the previous time-step is used, but it should be the one of the current time-step. Then, one had to solve a non-linear problem.
    """
    def __init__(self, data, logger):
        # use __init__ from Fenics class
        super().__init__(data, logger)
        self.logger.info("Starting Joule Heating simulation.")
        ###
        # Physics part
        ###
        self.stationary = 'timesteps' not in self.data
        self.static = 'frequencies' not in self.data
        self._prepare_output()
        if not self.static:
            self._set_frequencies()
            self.eps0 = 8.854187817e-12
            if len(self.frequencies) > 1:
                raise Exception("Currently not implemented to have more than one frequency!")
        if not self.stationary:
            self._set_timestepping()
        self._set_material_constants()
        ###
        # FEM part
        ###
        self._set_function_space()
        self._set_function_space_EL()
        self._prepare_boundaries()
        self._set_material_constants()
        self.logger.debug("Set DirichletBC")
        self._set_boundaries()

        if self.stationary:
            self.N = 1
            self.t = 0
            self.logger.debug("Stationary problem will be set.")
        else:
            self._set_timestepping()
            self.logger.debug("Time-dependent problem will be set.")
        self._set_problem()

        self.logger.debug("Problem set.")
        self._setup_solver()
        self._solve_problem()
        self._close_output()

    def _set_function_space_EL(self):
        if not self.static:
            Ec = self.element * self.element
        else:
            Ec = self.element
        self.VEL = d.FunctionSpace(self.mesh.mesh, Ec)
        self.logger.info("Number of DOFs for electric part : {}".format(self.V.dim()))

    def _set_boundaries(self):
        if self.static:
            self.boundariesEM.set_DirichletBC(self.VEL, self.logger)
        else:
            self.boundariesEM.set_DirichletBC(self.VEL, self.logger, subspace=0)
            self.boundariesEM.set_DirichletBC(self.VEL, self.logger, subspace=1, value=0.0)
        self.boundariesT.set_DirichletBC(self.V, self.logger)

    def _prepare_boundaries(self):
        """ need to redefine function as we have two kinds of boundaries """
        self.boundariesEM = Boundary(self.data['boundaries']['electric'], self.mesh.facets, self.mesh.facetinfo, self.mesh.cells, self.mesh.subdomaininfo)
        self.boundariesT = Boundary(self.data['boundaries']['thermal'], self.mesh.facets, self.mesh.facetinfo, self.mesh.cells, self.mesh.subdomaininfo)

    def _set_problem(self):
        self.logger.debug("Set problem")
        # for now hard-coded
        self.source_term = d.Constant(0.0)
        self.logger.debug("trial functions")
        if self.stationary:
            self.T_n = None
        else:
            # Define initial value, which will later be the solution from the previous (n-th) time step
            self.T_n = d.Function(self.V)
            self.T_n.rename('temperature', 'temperature')
            self.T_n = d.interpolate(self.T_0, self.V)

        # EM part ###
        if self.static:
            psi = d.TrialFunction(self.VEL)
            v = d.TestFunction(self.VEL)
            a = d.inner(self._sigma_T(self.T_n) * d.grad(psi), d.grad(v)) * self.dx
            L = ESRHS(self.source_term, v, self.dx)
            a, L = self.boundariesEM.set_RobinBCPoisson(a, L, psi, v, self.ds, 'R', 'V_ref', self.logger)
        if self.static is False:
            frequency = self.frequencies[0]
            self.permittivity_f = self.permittivity.copy(deepcopy=True)
            self.permittivity_f.vector()[:] = self.permittivity.vector()[:] * self.eps0 * frequency * 2. * np.pi
            u_r, u_i = d.TrialFunction(self.VEL)
            v_r, v_i = d.TestFunction(self.VEL)
            a = EQSLHS(self.conductivity, self.permittivity_f, u_r, u_i, v_r, v_i, self.dx)
            L = EQSRHS(self.source_term, v_r, v_i, self.dx)
            a, L = self.boundariesEM.set_RobinBCPoissonComplex(a, L, u_r, u_i, v_r, v_i, self.ds, 'Y_real', 'Y_imag', 'V_0r', 'V_0i', self.logger)

        self.psi = d.Function(self.VEL)
        self.EMproblem = d.LinearVariationalProblem(a, L, self.psi, self.boundariesEM.bc)

        # thermal part ###
        v = d.TestFunction(self.V)
        T_n1 = d.TrialFunction(self.V)

        if self.stationary:
            self.N = 1
            self.t = 0
            a = HeatLHSstat(T_n1, v, self.thermal_conductivity, self.dx)
            L = HeatRHSstat(v, self.dx)
            if self.static:
                L += self._sigma_T(self.T_n) * d.dot(d.grad(self.psi), d.grad(self.psi)) * v * self.dx
            else:
                L += self._sigma_T(self.T_n) * (d.dot(d.grad(self.psi.sub(0)), d.grad(self.psi.sub(0))) + d.dot(d.grad(self.psi.sub(1)), d.grad(self.psi.sub(1)))) * v * self.dx

        else:
            a = HeatLHStime(T_n1, v, self.density, self.heat_capacity, self.thermal_conductivity, self.dx, self.dt)
            L = HeatRHStime(self.T_n, v, self.density, self.heat_capacity, self.thermal_conductivity, self.dx)
            # coupling term
            if self.static:
                L += self.dt * self._sigma_T(self.T_n) * d.dot(d.grad(self.psi), d.grad(self.psi)) * v * self.dx
            else:
                L += self.dt * self._sigma_T(self.T_n) * (d.dot(d.grad(self.psi.sub(0)), d.grad(self.psi.sub(0))) + d.dot(d.grad(self.psi.sub(1)), d.grad(self.psi.sub(1)))) * v * self.dx

        a, L = self.boundariesT.set_RobinBCPoisson(a, L, T_n1, v, self.ds, 'h', 'T_ref', self.logger)
        self.T_n1 = d.Function(self.V)  # solution for T
        self.T_n1.rename('temperature', 'temperature')
        self.Tproblem = d.LinearVariationalProblem(a, L, self.T_n1, self.boundariesT.bc)
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

    def _prepare_output(self):
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.logger.info("Will store solution to hdf5 format")
                    self.datafileHDF5EM = d.HDF5File(self.mesh.mesh.mpi_comm(), self.result_dir + "EMsolution" + self.study + ".h5", "w")
                    self.datafileHDF5T = d.HDF5File(self.mesh.mesh.mpi_comm(), self.result_dir + "Tsolution" + self.study + ".h5", "w")
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self.logger.info("Will store solution to xdmf format")
                    self.datafileXDMFT = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "Tsolution" + self.study + '.xdmf')
                    if not self.static:
                        self.datafileXDMFEM_real = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "EMsolution_real" + self.study + '.xdmf')
                        self.datafileXDMFEM_imag = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "EMsolution_imag" + self.study + '.xdmf')
                    else:
                        self.datafileXDMFEM = d.XDMFFile(self.mesh.mesh.mpi_comm(), self.result_dir + "EMsolution" + self.study + '.xdmf')

    def _close_output(self):
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5EM.close()
                    self.datafileHDF5T.close()
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    if not self.static:
                        self.datafileXDMFEM_real.close()
                        self.datafileXDMFEM_imag.close()
                    else:
                        self.datafileXDMFEM.close()
                    self.datafileXDMFT.close()

    def _write_EM_output(self, n, append=False):
        if 'XDMF' in self.data['output']:
            if self.data['output']['XDMF']:
                if not self.static:
                    self.datafileXDMFEM_real.write_checkpoint(self.psi.sub(0), 'potential real', self.t, append=append)
                    self.datafileXDMFEM_imag.write_checkpoint(self.psi.sub(1), 'potential imag', self.t, append=append)
                else:
                    self.datafileXDMFEM.write_checkpoint(self.psi, 'potential', self.t, append=append)
        if 'HDF5' in self.data['output']:
            if self.data['output']['HDF5'] is True:
                self.datafileHDF5EM.write(self.psi, "/V_{}".format(n))

    def _write_T_output(self, T, n, append=False):
        if 'XDMF' in self.data['output']:
            if self.data['output']['XDMF']:
                self.datafileXDMFT.write_checkpoint(T, 'temperature', self.t, append=append)
        if 'HDF5' in self.data['output']:
            if self.data['output']['HDF5'] is True:
                self.datafileHDF5T.write(T, "/T_{}".format(n))

    def _solve_problem(self):
        append = False
        for n in range(self.N):
            if 'thermal_resistivity' in self.data:
                # solve EM problem for every time step
                self.EMsolver.solve()
                self.logger.debug("Solved EM problem at n=" + str(n))

            else:
                # solve once and forever, only solution for temperature for each time step
                if n == 0:
                    self.EMsolver.solve()
                    self.logger.debug("Solved EM problem at n=" + str(n))
                    self._write_EM_output(n, append=append)

            # Compute solution
            self.Tsolver.solve()
            self.logger.debug("Solved thermal problem at n=" + str(n))
            # if stationary: leave here
            if self.N == 1:
                break

            # write previous solution to file
            self._write_T_output(self.T_n, n, append=append)
            # Update previous solution
            self.T_n.assign(self.T_n1)

            # Update current time
            self.t += self.dt
            # to write time series
            append = True
        self._write_T_output(self.T_n1, self.N, append=append)
        if 'thermal_resistivity' in self.data:
            self._write_EM_output(n, append=append)

    def _set_material_constants(self):
        # global definitions
        try:
            temp_info = self.data['temperature']
            self.T_0 = d.Constant(temp_info['initial'])
            if 'thermal_resistivity' in self.data:
                self.T_ref = d.Constant(temp_info['reference'])
        except KeyError:
            self.logger.warning("You must provide the temperature in a dictionary containing \" initial \" and \" reference \" (time-dependent case) as keys.")
        # material-specific definitions
        if 'resistivity' in self.data:
            self._set_material_constant('resistivity')
        elif 'conductivity' in self.data:
            self._set_material_constant('conductivity')
        else:
            raise Exception("Conductivity or resistivity must be known!")
        if self.static is False:
            if 'permittivity' in self.data:
                self._set_material_constant('permittivity')
            else:
                raise Exception("Permittivity must be known!")

        if self.stationary is False:
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
