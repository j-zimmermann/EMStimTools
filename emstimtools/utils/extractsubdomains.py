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


def extractSubdomainsandFacets(filename, model_dim, exportyml=None):
    """
    Take msh file and extract subdomains and facets.
    """
    mshfile = open(filename, "r")
    subdomain_data = {}
    facet_data = {}
    while True:
        line = mshfile.readline()
        if not line:
            # EOF
            break
        if "PhysicalNames" in line:
            line = mshfile.readline()
            physicalnames = int(line)
            # go through section and fill dictionaries
            for _ in range(physicalnames):
                line = mshfile.readline()
                key = line.split(" ")[2].replace('"', "").replace("\n", "")
                # get group and dimension
                group = int(line.split(" ")[1])
                dim = int(line.split(" ")[0])
                if dim == model_dim:
                    subdomain_data[key] = group
                elif dim == model_dim - 1:
                    facet_data[key] = group
            line = mshfile.readline()
            assert line.strip() == "$EndPhysicalNames"
            break
    if exportyml is None:
        return subdomain_data, facet_data
    else:
        stream = open(exportyml, 'w')
        yaml.dump({"subdomains": subdomain_data}, stream, default_flow_style=False)
        yaml.dump({"facets": facet_data}, stream, default_flow_style=False)
