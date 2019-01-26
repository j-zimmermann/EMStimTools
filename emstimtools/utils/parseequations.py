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

from collections import OrderedDict


def parseEquations(dictwithstrings, dictwithvalues, dimension):
    """
    take a dict with equations written as strings and turn them into values
    """
    points_eval = dict()
    for i in range(dimension):
        for p in dictwithstrings:
            coords = dictwithstrings[p]
            assert(len(coords) == dimension)
            values = []
            for j in range(dimension):
                iter_str = coords[j].split()
                for idx, element in enumerate(iter_str):
                    if element in dictwithvalues:
                        iter_str[idx] = str(dictwithvalues[element])
                values.append(eval(''.join(iter_str)))

            points_eval[p] = values
    points_eval = OrderedDict(sorted(points_eval.items(), key=lambda t: t[0]))
    return points_eval
