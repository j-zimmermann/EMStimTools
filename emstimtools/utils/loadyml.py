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

import yaml


def load_yaml_file(filename):
    stream = open(filename, "r")
    data = yaml.load(stream)
    return data


def return_study(data):
    study = None
    try:
        study = data['study']
    except (KeyError, ValueError):
        print("Please name your study!")
    return study


def return_dimension(data):
    dimension = None
    try:
        dimension = data['dimension']
    except (KeyError, ValueError):
        print("Please specify the dimension!")
    return dimension


def return_geometry(data):
    """
    .. todo:: write a ordered dict, that accepts arrays and allow geometry parameters to be written as array. then use 'salome -t xy.py args:a,b,c'

    :param list geometry: Contains the names of all parts of the geometry are given as specified in SALOME.
    :param dict geometryvalues: Contains values and takes the values of geometry as keys.
    :return: geometry and geometryvalues.
    """
    geometry = None
    geometryvalues = None
    try:
        geometryvalues = data['geometryvalues']
    except (KeyError, ValueError):
        print("No geometry info")

    try:
        geometry = data['geometry']
    except (KeyError, ValueError):
        print("No geometry info")
    return geometry, geometryvalues


def return_subdomains(data):
    subdomains = None
    facets = None
    try:
        subdomains = data['subdomains']
    except (KeyError, ValueError):
        print("no subdomain info")
    try:
        facets = data['facets']
    except (KeyError, ValueError):
        print("no facet info")
    return subdomains, facets


def return_boundaries(data):
    boundaries = None
    try:
        boundaries = data['boundaries']
    except (KeyError, ValueError):
        print("no boundary info")
    return boundaries


def return_points_field(data):
    evaluatepoints = None
    try:
        evaluatepoints = data['evaluatefield']
    except (KeyError, ValueError):
        print("no points given")
    return evaluatepoints


def return_points_pot(data):
    evaluatepoints = None
    try:
        evaluatepoints = data['evaluatepot']
    except (KeyError, ValueError):
        print("no points given")
    return evaluatepoints


def return_refinement(data):
    refinement = None
    try:
        refinement = data['refinement']
    except (KeyError, ValueError):
        print("no refinement chosen")
    return refinement


def return_conductivity(data):
    conductivity = None
    try:
        conductivity = data['conductivity']
    except (KeyError, ValueError):
        print("no conductivity data")
    return conductivity


def return_permittivity(data):
    permittivity = None
    try:
        permittivity = data['permittivity']
    except (KeyError, ValueError):
        print("no permittivity data")
    return permittivity


def return_solver(data):
    solver = None
    try:
        solver = data['solver']
    except (KeyError, ValueError):
        print("no solver data")
    return solver
