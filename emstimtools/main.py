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

import logging
from shutil import copyfile
from .utils.loadyml import (load_yaml_file, return_study, return_geometry, return_dimension)
from .utils.changesalomefile import changeSALOMEfile
from .utils.runsalomefile import callSalome
from .utils.convertmesh import (convertMEDfile, convertmshfile, createh5file)
from .fenics.EQS import EQS
from .fenics.JouleHeating import JouleHeating
from .fenics.ES import ES
from .fenics.HeatEquation import Heat
from .dataanalysis.data import SensitivityAnalysis


class Simulation(object):
    """
    This is the main simulation class.
    It contains all necessary information.
    The most important parameter is the yamlfile.
    Information on how to prepare the yamlfile and how it is processed can be found in :ref:`yaml`.
    It will be used to generate a dictionary called :attr:`data`.
    After initializing an instance of 'Simulation', the simulation can be run by
    calling the method :meth:`run`.

    Logging events is possible. Logs are written to standard output and/or a file called 'EMStimTools_:attr:`study`'.
    The default level of the logger is 'WARNING'.
    Another level can be chosen according to the `official documentation <https://docs.python.org/2/library/logging.html#levels>`_.
    The method :meth:`set_log_level` takes the level either as a string or an int value.
    Setting 'INFO' increases the output, and setting it to 'DEBUG' increases it even further.

    Possible keys for :attr:`data` are listed in the following. They are structured according to their role:

    #. *General information*

        .. csv-table:: General header
            :header: "**Key**",   "**Role**",   "**Format**"
            :widths: 15, 25, 10

            "study",  "name of the study", "str"
            "SALOMEfile", "python file from SALOME", "str"
            "physics",  "PDE(s) to solve, see :class:`.Fenics`", "str"
            "dimension", "geometrical dimension", "int"

    #. *Geometry*

        .. csv-table:: Information about geometry
            :header: "**Key**",   "**Role**",   "**Format**"
            :widths: 15, 25, 10

            "geometry", "list with names of SALOME notebook parameters", "list of str"
            "geometryvalues", "values for notebook parameters", "dict"
            "materials", "subdomains as specified in SALOME", "list of str"

    #. *General FEM*

        General FEM keys are given here. For more detailed documentation for the respective physics, please refer to the :class:`.Fenics` section.

        .. csv-table:: General settings for FEniCS
            :header: "**Key**",   "**Role**",   "**Format**"
            :widths: 15, 25, 10

            "element", "Finite Element, choose from `here <http://femtable.org/>`_", "str"
            "degree", "degree of finite element", "int"
            "solver", "direct access to `linear <https://fenicsproject.org/docs/dolfin/2018.1.0/cpp/d7/d34/LinearVariationalSolver_8h_source.html>`_ and `non-linear <https://fenicsproject.org/docs/dolfin/2018.1.0/cpp/dc/de4/NonlinearVariationalSolver_8h_source.html>`_ solver parameters in FEniCS", "dict"

    #. Output

        .. csv-table:: Output possibilities of EMStimTools
            :header: "**Key**",   "**Role**",   "**Format**"
            :widths: 15, 25, 10

            "output", "choose option for writing solution. HDF5 or XDMF.", "dict, e.g. XDMF : yes"
            "properties", "compute derived quantities, output format as chosen above", "dict"
            "plot", "plot e.g. mesh", "dict"
            "postprocess", "choose options for postprocessing, see :meth:`.Fenics._postprocess`", "dict"

    .. todo:: work on MPI possibility (at least for FEniCS). The idea is: one run only for creating the hdf5 file and then a second run for the simulation and possible postprocessing in a 3rd call (nasty but inevitable). maybe fenicstools could help here.
    """

    def __init__(self, yamlfile):
        if not str(yamlfile).endswith('.yml'):
            raise Exception("You must provide a YAML input file.")
        self.yamlfile = yamlfile
        # to check whether input file is there
        self.mesh_tag = None
        # get dictionary
        self._get_information()

        # add a log file
        self._add_logger()

    def run(self):
        """
        If everything is alright, the simulation can be started from this command.
        The work flow is:

        #. take SALOME template file and make it general.
        #. determine loops, create mesh accordingly
        """
        self._prepare_geometry()
        self._create_meshfile()
        self._run_study()
        self._data_analysis()
        self.logger.info("Run successful.")

    def _get_information(self):
        """ Load in yaml file.
            Stores dictionary 'data'
        """
        self.data = load_yaml_file(self.yamlfile)
        self.study = return_study(self.data)
        self.dimension = return_dimension(self.data)
        if 'salome_parameters' not in self.data:
            self.data['salome_parameters'] = None

    def _add_logger(self):
        self.logger = logging.getLogger('EMStimTools_' + str(self.study))
        self.fh = logging.FileHandler(str(self.study) + '.log', encoding='utf-8')
        self.logger.addHandler(self.fh)
        self.ch = logging.StreamHandler()
        self.ch.setLevel('INFO')
        self.logger.addHandler(self.ch)
        self.logger.info("Loaded YAML file successfully.")
        self.logger.info("Going to run study\" " + str(self.study) + "\"")

    def _salomefile_processing(self):
        self.SALOMEfile = self.data['SALOMEfile']
        self.geometryparameter, self.geometryvalues = return_geometry(self.data)

        self.logger.debug("Read in SALOMEfile " + str(self.SALOMEfile))
        changeSALOMEfile(self.SALOMEfile, self.geometryparameter, self.geometryvalues, meshname=str(self.study) + '_mesh', salome_parameters=self.data['salome_parameters'])
        self.logger.info("SALOME file " + str(self.SALOMEfile) + " prepared")

    def _prepare_geometry(self):
        """
        Take a template SALOME file and make it general.
        """
        if 'SALOMEfile' in self.data:
            self._salomefile_processing()
        elif 'meshfile' in self.data:
            self.logger.info("You choose to run the study from a meshfile. SALOME will be skipped.")
            self.mesh_tag = True
        else:
            raise Exception("You must provide a SALOMEfile or meshfile")

    def _create_meshfile(self):
        if self.mesh_tag is None:
            self.logger.info('Run SALOME, now follows SALOME output:')
            callSalome(self.SALOMEfile)
        else:
            self.logger.info("Use MED file " + self.data['meshfile'])
            self.logger.info("It will not be destroyed but copied and saved.")
            if self.data['meshfile'] == str(self.study) + '_mesh.med':
                copyfile(self.data['meshfile'], self.data['meshfile'] + '.tmp')
                copytag = True
            else:
                copyfile(self.data['meshfile'], str(self.study) + '_mesh.med')
                copytag = False
        self.logger.info('now convert to HDF5 file:')
        self.logger.debug('MED file to MSH')
        convertMEDfile(str(self.study) + '_mesh.med')
        self.logger.debug('MSH file to XML')

        subdomains, facets = convertmshfile(str(self.study) + '_mesh.msh2', self.dimension)
        self._check_consistency(subdomains, facets)
        self.logger.debug("Subdomains")
        self.logger.debug(subdomains)
        self.logger.debug("Facets")
        self.logger.debug(facets)
        self.data.update({'subdomains': subdomains})
        self.data.update({'facets': facets})
        self.logger.debug('Updated dictionary')
        self.logger.debug(' XML file to HDF5')
        createh5file(str(self.study) + '_mesh.xml')
        if self.mesh_tag is not None:
            if copytag is True:
                copyfile(self.data['meshfile'] + '.tmp', self.data['meshfile'])
        self.logger.info('done')

    def _check_consistency(self, subdomains, facets):
        if 'materials' not in self.data:
            raise Exception("You must provide material information!\n Use a list of strings for that.")
        assert isinstance(subdomains, dict) is True
        assert isinstance(facets, dict) is True
        if not all(key in self.data['materials'] for key in list(subdomains.keys())):
            raise Exception("Materials do not match subdomains! Check your spelling and consistency of yaml and salome-file.\n Subdomains are: " + str(list(subdomains.keys())) + ", materials: " + str(self.data['materials']))

    def _run_study(self):
        self.logger.info("Running study for physics" + self.data['physics'])
        if self.data['physics'] == 'EQS':
            self.fenics_study = EQS(self.data, self.logger)
        elif self.data['physics'] == 'JouleHeating':
            self.fenics_study = JouleHeating(self.data, self.logger)
        elif self.data['physics'] == 'ES':
            self.fenics_study = ES(self.data, self.logger)
        elif self.data['physics'] == 'Heat':
            self.fenics_study = Heat(self.data, self.logger)

        else:
            raise Exception("Specify a study that exists!!!")
        self.logger.info("FEniCS run successful")

    def set_log_level(self, level, stream=False, filelog=False):
        """
        Choose an appropriate level for the logger,
        examples can be found `here <https://docs.python.org/2/library/logging.html#levels>`_.
        By default, the level for the entire logger is changed.

        :param bool stream: change level only for standard output if True
        :param bool filelog: change level for file log only if True
        """
        try:
            if not stream or filelog:
                self.ch.setLevel(level)
                self.fh.setLevel(level)
                self.logger.setLevel(level)
                print("Changed global log level to " + str(level))
            elif stream:
                self.ch.setLevel(level)
                print("Changed stream log level to " + str(level))
            elif filelog:
                self.fh.setLevel(level)
                print("Changed file log level to " + str(level))

        except ValueError:
            print("Use a correct logging level!")
            print("Valid values can be found here: https://docs.python.org/2/library/logging.html#levels")
            print("Continuing with previous value")
            pass

    def _data_analysis(self):
        """
        currently only sensitivity analysis for EQS available!

        .. todo:: implement more interfaces such as UQ

        """
        if self.data['physics'] != 'EQS':
            return
        if 'dataanalysis' in self.data:
            if 'sensitivityanalysis' in self.data['dataanalysis']:
                SensitivityAnalysis(self.data)

    def _prepare_dataoutput(self):
        """

        .. todo:: this is very preliminary

        """
        if 'sensitivityanalysis' in self.data['dataanalysis']:
            if 'solution' in self.data['dataanalysis']['sensitivityanalysis']:
                if 'evaluate_solution' not in self.data:
                    self.data['evaluate_solution'] = {}
                self.data['evaluate_solution']['sensitivity_point'] = self.data['dataanalysis']['sensitivityanalysis']['solution']
            else:
                self.logger.warning('Not implemented')
        return
