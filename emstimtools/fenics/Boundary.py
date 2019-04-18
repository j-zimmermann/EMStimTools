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


class Boundary(object):
    """
    Define Boundary conditions.
    This class takes the :attr:`boundary_dict`, where
    all information about the boundaries is stored.
    The data is structured as:
    - :attr:`facets` : MeshFunction containing the indices corresponding to the facets
    - :attr:`facet_dict` : Mapping from name of facet to index
    - :attr:`subdomains` : MeshFunction containing the indices of the subdomains
    - :attr:`subdomain_dict` : Mapping from name of subdomain to index
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
        if 'Dirichlet' in self.boundary_dict:
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

    def set_RobinBCPoisson(self, a, L, u, v, ds, p_name, q_name, logger, dt=False):
        r"""
        Function to set Robin BC for Poisson like equation.
        Find more information in the FEnics tutorial chapter 1.5.3
        A Robin BC for means that we have an expression like:

        .. math::

            - k \frac{\partial u}{\partial n} = p (u - q) \enspace ,

        where q is a constant value, u is the unknown and p is a constant factor as well. k might be material property.
        Essentially, it enters the weak formulation by changing the LHS like

        .. math::

            a = a + \int_\Gamma p u v \mathrm{d}s

        and the RHS like:

        .. math::

            L = L + \int_\Gamma p q v \mathrm{d}s

        """
        if 'Robin' in self.boundary_dict:
            self.boundary_info = self.boundary_dict['Robin']
            for boundary_key in self.boundary_info:
                self.logger.debug('Set Robin BC for ' + boundary_key + ' with info ' + self.boundary_info[boundary_key])
                index = self.facet_dict[boundary_key]
                try:
                    p = self.boundary_info[boundary_key][p_name]
                    q = self.boundary_info[boundary_key][q_name]
                    if dt is not False:
                        p = p * dt
                except KeyError:
                    raise("You have not specified all parameters needed for the Robin BC!")
                a += p * u * v * ds(index)
                L += p * q * v * ds(index)
        return a, L

    def set_RobinBCPoissonComplex(self, a, L, u_r, u_i, v_r, v_i, ds, p_r_name, p_i_name, q_r_name, q_i_name, logger):
        """
        use Robin BC for complex-valued problems (i.e. with mixed function spaces)
        """
        if 'Robin' in self.boundary_dict:
            self.boundary_info = self.boundary_dict['Robin']
            for boundary_key in self.boundary_info:
                self.logger.debug('Set Robin BC for ' + boundary_key + ' with info ' + self.boundary_info[boundary_key])
                index = self.facet_dict[boundary_key]
                try:
                    p_r = self.boundary_info[boundary_key][p_r_name]
                    p_i = self.boundary_info[boundary_key][p_i_name]
                    q_r = self.boundary_info[boundary_key][q_r_name]
                    q_i = self.boundary_info[boundary_key][q_r_name]
                except KeyError:
                    raise("You have not specified all parameters needed for the Robin BC!")
                a += p_r * u_r * v_r * ds(index)
                a += p_r * u_r * v_i * ds(index)
                a -= p_i * u_i * v_r * ds(index)
                a -= p_i * u_i * v_i * ds(index)
                a += p_r * u_i * v_r * ds(index)
                a -= p_r * u_i * v_i * ds(index)
                a += p_i * u_r * v_r * ds(index)
                a -= p_i * u_r * v_i * ds(index)

                L -= p_r * q_r * v_r * ds(index)
                L -= p_r * q_r * v_i * ds(index)
                L += p_i * q_i * v_r * ds(index)
                L += p_i * q_i * v_i * ds(index)
                L -= p_r * q_i * v_r * ds(index)
                L += p_r * q_i * v_i * ds(index)
                L -= p_i * q_r * v_r * ds(index)
                L += p_i * q_r * v_i * ds(index)
        return a, L
