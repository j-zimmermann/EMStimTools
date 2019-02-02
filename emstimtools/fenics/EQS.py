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
compute EM fields in stimulation tools by solving the EQS equation.

"""
import numpy as np
import dolfin as d
from .fenics import Fenics
from itertools import product


def dict_product(d):
    keys = d.keys()
    for element in product(*d.values()):
        yield dict(zip(keys, element))


class EQS(Fenics):
    """
    Define a class, where an entire EQS run will be done.
    We need the following attributes:

    - :attr:`frequencies` : frequencies, either array or single frequency
    - :attr:`conductivity_info`  : conductivities
    - :attr:`permittivity_info` : permittivities
    - :attr:`V` : FunctionSpace as a mixed function space for real and imaginary part


    General Information can be obtained by using:
    - :meth:`get_info`

    We use the following steps:

    #. set the frequencies for which the model will be evaluated
    #. set the material constants, i.e. conductivity and permittivity
    #. set function space, here a mixed function space for real and imaginary part
    #. set BCs. The imaginary part is always set to 0.
    #. set problem, i.e. bilinear form
    #. set the solver as defined in YAML input file
    #. solve the problem. for an array job: evaluate all possible parameter combinations.
    #. during the solve, the solution may be evaluated at certain points
    #. prepare the output and print it

    .. todo:: integrate complex number support that is partly available in dolfin-x
    .. todo:: add source term in case of charge sources
    .. todo:: add support for storing the field, see e.g. :meth:`get_fieldnorm`

    """
    def __init__(self, data, logger):
        # use __init__ from Fenics class
        super(EQS, self).__init__(data, logger)
        self.logger.info("Starting EQS simulation.")
        ###
        # Physics part
        ###
        self._set_frequencies()
        self._set_material_constants()
        ###
        # FEM part
        ###
        self._set_function_space()
        self._prepare_boundaries()
        self.boundaries.set_DirichletBC(self.V, self.logger, subspace=0)
        # boundary for imaginary part always 0
        self.boundaries.set_DirichletBC(self.V, self.logger, subspace=1, value=0.0)
        self._set_problem()
        self.logger.debug("Problem set.")
        self._setup_solver()
        ###
        # solution part
        ###
        # if there is an array job, iterate through
        if self.array_tag is False:
            self._prepare_output()
            self._solve_problem()
            self._close_output()
        else:
            self.combinations = len(self.parameter_combinations)
            for _ in range(self.combinations):  # run for all possible parameter combinations
                self._update_material_constant()
                self._prepare_output()
                self._solve_problem()
                self._close_output()
        self.logger.debug('checking for postprocessing steps')
        if 'postprocess' in self.data:
            self._postprocess()

    def _update_material_constant(self):
        """ also update output files"""
        tmp = self.parameter_combinations.pop()
        self.data['conductivity'] = tmp[0]
        self.data['permittivity'] = tmp[1]
        self._set_material_constant('conductivity')
        self._set_material_constant('permittivity')
        self.logger.info('Running for conductivity')
        self.logger.info(tmp[0])
        self.logger.info('Running for permittivity')
        self.logger.info(tmp[1])

    def _set_material_constants(self):
        """
        when the material constants are set, it is first determined if an array job has to run.
        """
        self.eps0 = 8.854187817e-12
        if self.data['conductivity'] is not None:
            self.conductivity_info = self.data['conductivity']
        else:
            raise Exception("Conductivity must be known!")

        for key in self.conductivity_info:
            if isinstance(self.conductivity_info[key], list):
                self.array_tag = True
                self.array_iter = 0

        if self.data['permittivity'] is not None:
            self.permittivity_info = self.data['permittivity']
            if self.array_tag is False:
                for key in self.permittivity_info:
                    if isinstance(self.permittivity_info[key], list):
                        self.array_tag = True
        else:
            raise Exception("Permittivity must be known!")
        self._set_material_constant('conductivity')
        # beware that eps0 and frequency need to be considered!!!
        self._set_material_constant('permittivity')
        # have second array for this
        self.permittivity_f = self.permittivity.copy(deepcopy=True)

        if self.array_tag is True:
            # delete conductivity array, will be written later
            self.data['conductivity_array'] = self.data.pop('conductivity')
            self.conductivity_array = []
            self.conductivity_array_values = [[]]
            # append all arrays
            for key in self.conductivity_info:
                self.conductivity_array.append(key)
                if isinstance(self.conductivity_info[key], list):
                    self.conductivity_array_values.append(self.conductivity_info[key])
                else:
                    self.conductivity_array_values.append([self.conductivity_info[key]])
            self.conductivity_array_values.pop(0)
            d = dict(zip(self.conductivity_array, self.conductivity_array_values))
            self.conductivity_arrays = list(dict_product(d))
            # delete permittivity dict, will be written later
            self.data['permittivity_array'] = self.data.pop('permittivity')
            permittivity_array = []
            permittivity_array_values = [[]]
            # append all arrays
            for key in self.permittivity_info:
                permittivity_array.append(key)
                if isinstance(self.permittivity_info[key], list):
                    permittivity_array_values.append(self.permittivity_info[key])
                else:
                    permittivity_array_values.append([self.permittivity_info[key]])

            permittivity_array_values.pop(0)
            d = dict(zip(permittivity_array, permittivity_array_values))
            self.permittivity_arrays = list(dict_product(d))
            # all possible combinations:
            self.parameter_combinations = list(product(self.conductivity_arrays, self.permittivity_arrays))

    def _set_function_space(self):
        Ec = self.element * self.element
        self.V = d.FunctionSpace(self.mesh.mesh, Ec)
        self.logger.info("Number of DOFs : {}".format(self.V.dim()))

    def _set_problem(self):
        # for now hard-coded
        self.logger.debug("Set problem")
        self.source_term = d.Constant(0.0)
        self.logger.debug("trial functions for real and imaginary part")
        u_r, u_i = d.TrialFunction(self.V)
        self.logger.debug("test functions for real and imaginary part")
        v_r, v_i = d.TestFunction(self.V)
        dx = d.dx

        a = (d.inner(d.grad(self.conductivity * u_r), d.grad(v_r)) * dx
             - d.inner(d.grad(self.permittivity_f * u_i), d.grad(v_r)) * dx
             - d.inner(d.grad(self.permittivity_f * u_r), d.grad(v_i)) * dx
             - d.inner(d.grad(self.conductivity * u_i), d.grad(v_i)) * dx
             + d.inner(d.grad(self.conductivity * u_r), d.grad(v_i)) * dx
             - d.inner(d.grad(self.permittivity_f * u_i), d.grad(v_i)) * dx
             + d.inner(d.grad(self.permittivity_f * u_r), d.grad(v_r)) * dx
             + d.inner(d.grad(self.conductivity * u_i), d.grad(v_r)) * dx)
        L = -(self.source_term * v_r + self.source_term * v_i) * dx
        self.logger.debug("vector for solution")
        self.u = d.Function(self.V)
        self.problem = d.LinearVariationalProblem(a, L, self.u, self.boundaries.bc)

    def _solve_problem(self):
        append = False
        for frequency in self.frequencies:
            self.logger.info('Solving for {}.'.format(frequency))
            # use eps = frequency * eps_r * eps0
            self.permittivity_f.vector()[:] = self.permittivity.vector()[:] * self.eps0 * frequency * 2. * np.pi
            self.solver.solve()
            # self.u_r, self.u_i = self.u.split(deepcopy=True)
            try:
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5.write(self.u.sub(0), "/real part {}".format(frequency))
                    # self.datafileHDF5.write(self.u_r, "/real part {}".format(frequency))
                    self.datafileHDF5.write(self.u.sub(1), "/imaginary part {}".format(frequency))
                    # self.datafileHDF5.write(self.u_i, "/imaginary part {}".format(frequency))
            except KeyError:
                pass
            try:
                if self.data['output']['XDMF'] is True:
                    self.datafileXDMFreal.write_checkpoint(self.u.sub(0), 'potential real part', frequency, append=append)
                    self.datafileXDMFimag.write_checkpoint(self.u.sub(1), 'potential imag part', frequency, append=append)
            except KeyError:
                pass
            self._write_derived_quantities(frequency, append)
            self._point_evaluation(frequency)
            append = True

    def _write_derived_quantities(self, frequency, append):
        self.logger.debug("Computing derived quantities")
        if 'properties' in self.data:
            try:
                if self.data['properties']['E-Field'] is True:
                    self.get_field()
                    self.dataXDMFEreal.write_checkpoint(self.Efield_real, 'E-Field real part', frequency, append=append)
                    self.dataXDMFEimag.write_checkpoint(self.Efield_imag, 'E-Field imag part', frequency, append=append)
            except KeyError:
                pass
            try:
                if self.data['properties']['E-Field-norm'] is True:
                    self.get_fieldnorm()
                    self.dataXDMFEnorm.write_checkpoint(self.normEr, 'E-Field norm', frequency, append=append)
            except KeyError:
                pass

    def get_info(self):
        info_str = """\
        # The simulation is run for:
        # - frequencies: {}
        """.format(str(self.frequencies))
        info_str += """\
        # - conductivities:
        """
        str1 = ''
        str2 = ''
        for i in self.mesh.subdomaininfo:
            str1 += "# {} : {}\n".format(i, self.data['conductivity'][i])
            str2 += "# {} : {}\n".format(i, self.data['permittivity'][i])
        info_str += str1
        info_str += """\
        # - permittivities
        """
        info_str += str2
        info_str += """\
        # - boundaries {}
        """.format(self.boundaries.boundary_dict)
        return info_str

    def get_fieldnorm(self):
        """
        compute field norm.
        tried to solve it by means of projection => inaccurate

        .. todo:: decide whether field should be of same order as solution
        .. todo:: decide whether solver in `d.project` call should always be `mumps`
        """
        self.Vnorm = d.FunctionSpace(self.mesh.mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        self.normEr = d.Function(self.Vnorm)
        if not hasattr(self, 'Efield_real'):
            self.get_field()
        Ex_r, Ey_r, Ez_r = self.Efield_real.split(deepcopy=True)
        Ex_i, Ey_i, Ez_i = self.Efield_imag.split(deepcopy=True)
        norm_r = Ex_r.vector() * Ex_r.vector() + Ey_r.vector() * Ey_r.vector() + Ez_r.vector() * Ez_r.vector()
        norm_i = Ex_i.vector() * Ex_i.vector() + Ey_i.vector() * Ey_i.vector() + Ez_i.vector() * Ez_i.vector()
        normsquared = norm_r + norm_i
        norm = np.sqrt(normsquared.get_local())
        self.normEr.vector().set_local(norm)
        self.normEr.vector().apply('')
        self.logger.debug("Computed norm of E-Field")

    def get_field(self):
        """
        compute field

        .. todo:: find better representation for real and imaginary part
        .. todo:: decide whether field should be of same order as solution
        """
        self.logger.debug("Entering Field computation")
        if self.data['degree'] > 1:
            self.Vector = d.VectorFunctionSpace(self.mesh.mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        else:
            raise Exception("Use for this computation a 2nd or higher order element")

        self.Efield_real = d.project(-d.grad(self.u.sub(0)), self.Vector, solver_type='mumps')
        self.Efield_imag = d.project(-d.grad(self.u.sub(1)), self.Vector, solver_type='mumps')
        self.Efield_real.rename('E-field real', 'E-field real')
        self.Efield_imag.rename('E-field imag', 'E-field imag')
        self.logger.debug("Computed E-Field")

    def _prepare_output(self):
        """
        write to HDF5 or XDMF file if specified
        if solution at points is evaluated, also write to txt files
        """
        if 'output' in self.data:
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.logger.info("Will store solution to hdf5 format")
                    if self.array_tag is False:
                        filestring = self.result_dir + "solution" + self.study + ".h5"
                    else:
                        filestring = self.result_dir + "solution" + self.study + "_" + str(self.array_iter) + ".h5"
                    self.datafileHDF5 = d.HDF5File(self.mesh.mesh.mpi_comm(), filestring, "w")
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self.logger.info("Will store solution to xdmf format")
                    self._open_xdmf_output_solution()
            if 'properties' in self.data:
                try:
                    if self.data['properties']['E-Field'] is True:
                        self._open_output_efield()
                except KeyError:
                    pass
                try:
                    if self.data['properties']['E-Field-norm'] is True:
                        self._open_output_efieldnorm()
                except KeyError:
                    pass

        # txt files section
        if self.evaluate_points_solution is not None:
            if self.array_tag is False:
                filestring = self.result_dir + 'solution_' + self.study + '.dat'
            else:
                filestring = self.result_dir + 'solution_' + self.study + '_' + str(self.array_iter) + '.dat'
            self.datafileTXT = open(filestring, 'w')
            self.datafileTXT.write(self.get_info())
            self.datafileTXT.write('\n')
            self.datafileTXT.write("# output at certain points, for each point real and imaginary part plus absolute value\n")
            self.datafileTXT.write("# frequency")
            for i in self.evaluate_points_solution:
                self.datafileTXT.write(" " + i + " " + str(self.evaluate_points_solution[i]))
            self.datafileTXT.write("\n")

        # increment counter
        if self.array_iter is not None:
            self.array_iter += 1

    def _open_xdmf_output_solution(self, add=''):
        if self.array_tag is False:
            filestring = self.result_dir + "solution_real" + self.study + add + ".xdmf"
        else:
            filestring = self.result_dir + "solution_real" + self.study + "_" + str(self.array_iter) + add + ".xdmf"
        self.logger.debug('opening xdmffile ' + filestring)
        self.datafileXDMFreal = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)
        if self.array_tag is False:
            filestring = self.result_dir + "solution_imag" + self.study + add + ".xdmf"
        else:
            filestring = self.result_dir + "solution_imag" + self.study + "_" + str(self.array_iter) + add + ".xdmf"

        self.logger.debug('opening xdmffile ' + filestring)
        self.datafileXDMFimag = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)

    def _open_output_efield(self, add=''):
        """
        for the output we can add information, e.g. about the respective subdomain
        """
        if self.array_tag is False:
            filestring = self.result_dir + "efield_real" + self.study + add + ".xdmf"
        else:
            filestring = self.result_dir + "efield_real" + self.study + "_" + str(self.array_iter) + add + ".xdmf"
        self.logger.debug('opening xdmffile ' + filestring)
        self.dataXDMFEreal = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)
        if self.array_tag is False:
            filestring = self.result_dir + "efield_imag" + self.study + add + ".xdmf"
        else:
            filestring = self.result_dir + "efield_imag" + self.study + "_" + str(self.array_iter) + add + ".xdmf"
        self.logger.debug('opening xdmffile ' + filestring)
        self.dataXDMFEimag = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)

    def _open_output_efieldnorm(self, add=''):
        """
        for the output we can add information, e.g. about the respective subdomain
        """
        if self.array_tag is False:
            filestring = self.result_dir + "efieldnorm" + self.study + add + ".xdmf"
        else:
            filestring = self.result_dir + "efieldnorm" + self.study + "_" + str(self.array_iter) + add + ".xdmf"
        self.logger.debug('opening xdmffile ' + filestring)
        self.dataXDMFEnorm = d.XDMFFile(self.mesh.mesh.mpi_comm(), filestring)

    def _close_output(self):
        """
        close all possible output after successfull solution and before postprocessing starts.
        """
        self.logger.debug('cleaning up')
        if 'output' in self.data:
            self.logger.debug("going to close output")
            if 'HDF5' in self.data['output']:
                if self.data['output']['HDF5'] is True:
                    self.datafileHDF5.close()
                    self.logger.debug("Closed HDF5 file(s)")
            if 'XDMF' in self.data['output']:
                if self.data['output']['XDMF'] is True:
                    self._close_solution_xdmf()
        if self.evaluate_points_solution is not None:
            self.datafileTXT.close()
            self.logger.debug('Closed TXT file')

        if 'properties' in self.data:
                try:
                    if self.data['properties']['E-Field'] is True:
                        self._close_efield_xdmf()
                except KeyError:
                    pass
                try:
                    if self.data['properties']['E-Field-norm'] is True:
                        self._close_efieldnorm_xdmf()
                except KeyError:
                    pass

    def _close_solution_xdmf(self):
        self.datafileXDMFreal.close()
        self.datafileXDMFimag.close()
        self.logger.debug("Closed XDMF files")

    def _close_efield_xdmf(self):
        self.dataXDMFEreal.close()
        self.dataXDMFEimag.close()
        self.logger.debug('Closed E-Field files')

    def _close_efieldnorm_xdmf(self):
        self.dataXDMFEnorm.close()
        self.logger.debug('Closed E-Field norm files')

    def _write_projected_solution(self, i):
        """
        take domain :param str i: and project solution there
        """
        self.sub_mesh = d.SubMesh(self.mesh.mesh, self.mesh.cells, self.mesh.subdomaininfo[i])
        V_sub = d.FunctionSpace(self.sub_mesh, self.element)
        # real part, mumps to avoid memory overflow
        # u_sub = d.project(self.u_r, V_sub, solver_type='mumps')
        u_sub = d.project(self.u.sub(0), V_sub, solver_type='mumps')
        u_sub.rename('potential real part', 'potential real part sub')
        self.datafileXDMFreal.write(u_sub)
        # imag part
        u_sub.rename('potential imag part', 'potential imag part sub')
        # u_sub = d.project(self.u_i, V_sub, solver_type='mumps')
        u_sub = d.project(self.u.sub(1), V_sub, solver_type='mumps')
        self.datafileXDMFimag.write(u_sub)

    def _write_projected_efield(self, i):
        """
        take domain :param str i: and project e-field there
        """
        Vector_sub = d.VectorFunctionSpace(self.sub_mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        efield_sub = d.project(self.Efield_real, Vector_sub, solver_type='mumps')
        efield_sub.rename('E-Field real part', 'E-Field real part sub')
        self.dataXDMFEreal.write(efield_sub)
        efield_sub = d.project(self.Efield_imag, Vector_sub, solver_type='mumps')
        efield_sub.rename('E-Field imag part', 'E-Field imag part sub')
        self.dataXDMFEimag.write(efield_sub)
        self.logger.debug("Wrote e-field for domain " + i)

    def _write_projected_efieldnorm(self, i):
        """
        take domain :param str i: and project e-field norm there
        """
        Vnorm_sub = d.FunctionSpace(self.sub_mesh, self.data['properties']['project_element'], self.data['properties']['project_degree'])
        efieldnorm_sub = d.project(self.normEr, Vnorm_sub, solver_type='mumps')
        efieldnorm_sub.rename('E-Field norm', 'E-Field norm sub')
        self.dataXDMFEnorm.write(efieldnorm_sub)
        self.logger.debug("Wrote e-field norm for domain " + i)

    def _project_on_submeshes(self):
        """
        iterate through all domains that are provided in input file and write solution and if needed field.

        .. todo:: support for array jobs
        """
        self.logger.debug('Projecting EQS solutions on submeshes')
        domains = self.data['postprocess']['submesh']  # a list of strings
        for i in domains:
            self._open_xdmf_output_solution(add='_' + i)
            self._write_projected_solution(i)
            self._close_solution_xdmf()
            self._plot_submesh(i)
            if 'properties' in self.data:
                try:
                    if self.data['properties']['E-Field'] is True:
                        self._open_output_efield(add='_' + i)
                        self._write_projected_efield(i)
                        self._close_efield_xdmf()
                except KeyError:
                    pass
                try:
                    if self.data['properties']['E-Field-norm'] is True:
                        self._open_output_efieldnorm(add='_' + i)
                        self._write_projected_efieldnorm(i)
                        self._close_efieldnorm_xdmf()
                except KeyError:
                    pass

    def _point_evaluation(self, frequency):
        if self.evaluate_points_solution is not None:
            self.datafileTXT.write("%f" % frequency)
            for i in self.evaluate_points_solution:
                if(self.mesh.dimension == 2):
                    Phi_real = self.u_r(d.Point(self.evaluate_points_solution[i][0], self.evaluate_points_solution[i][1]))
                    Phi_imag = self.u_i(d.Point(self.evaluate_points_solution[i][0], self.evaluate_points_solution[i][1]))
                else:
                    Phi_real = self.u_r(self.evaluate_points_solution[i][0], self.evaluate_points_solution[i][1], self.evaluate_points_solution[i][2])
                    Phi_imag = self.u_i(self.evaluate_points_solution[i][0], self.evaluate_points_solution[i][1], self.evaluate_points_solution[i][2])
                self.datafileTXT.write(" %f %f %f" % (Phi_real, Phi_imag, d.sqrt(Phi_real * Phi_real + Phi_imag * Phi_imag)))
            self.datafileTXT.write("\n")
