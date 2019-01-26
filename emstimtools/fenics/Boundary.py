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
import random

class Boundary(object):
    """
    Define Boundary conditions.
    This class takes the :attr:`boundary_dict`, where
    all information about the boundaries is stored.
    Also the :attr:`facets` and their indices stored in :attr:`facet_dict` have to be known.
    """

    def __init__(self, boundary_dict, facets, facet_dict, subdomains, subdomain_dict):
        self.boundary_dict = boundary_dict
        self.facets = facets
        self.facet_dict = facet_dict
        self.subdomains = subdomains
        self.subdomain_dict = subdomain_dict
        self.bc = []

    def set_DirichletBC(self, V, logger, subspace=None, value=None):
        """
        Set Dirichlet Boundary Conditions.

        :param V: FEniCS Function Space
        :param logger: global logger
        :param subspace: for instance for EQS we need values for subspaces
        :param value: if there is no value in the input file, we can set it (needed for the imaginary part in EQS for instance)

        """
        if self.boundary_dict['Dirichlet'] is not None:
            self.boundary_info = self.boundary_dict['Dirichlet']
            if value is None:
                valuetag = None
                logger.debug("Set Dirichlet boundaries with info: {}".format(self.boundary_info))
            else:
                valuetag = True
                logger.debug("Set Dirichlet boundaries for all facets with value {}.".format(value))
            for boundary_key in self.boundary_info:
                if valuetag is None:
                    value = self.boundary_info[boundary_key]
                if subspace is None:
                    if boundary_key in self.facet_dict:
                        self.bc.append(d.DirichletBC(V, d.Constant(value), self.facets, self.facet_dict[boundary_key]))
                        logger.info("Set boundary for: {}, value {}, facet number {}".format(boundary_key, value, self.facet_dict[boundary_key]))
                    elif boundary_key in self.subdomain_dict:
                        raise Exception('BC must be facet!')
                        # self.bc.append(d.DirichletBC(V, d.Constant(value), self.subdomains, self.subdomain_dict[boundary_key]))
                        # logger.info("Set boundary for: {}, value {}, subdomain number {}".format(boundary_key, value, self.subdomain_dict[boundary_key]))
                    else:
                        raise Exception("BC neither match facets nor subdomains!!!!")
                else:
                    if boundary_key in self.facet_dict:
                        self.bc.append(d.DirichletBC(V.sub(subspace), value, self.facets, self.facet_dict[boundary_key]))
                        logger.info("Set boundary for: {} and subspace {}, value {}, facet number {}".format(boundary_key, subspace, value, self.facet_dict[boundary_key]))
                    elif boundary_key in self.subdomain_dict:
                        raise Exception('BC must be facet!')
                        # self.bc.append(d.DirichletBC(V.sub(subspace), value, self.subdomains, self.subdomain_dict[boundary_key]))
                        # logger.info("Set boundary for: {} and subspace {}, value {}, subdomain number {}".format(boundary_key, subspace, value, self.subdomain_dict[boundary_key]))

                    else:
                        raise Exception("BC neither match facets nor subdomains!!!!")

    def set_NeumannBC(self, logger):
        """
        Function to set Neumann Boundaries.
        .. todo:: need to implement Neumann Boundaries. Idea at: https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1005.html
        """
        return
