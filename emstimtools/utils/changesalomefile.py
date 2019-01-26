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

import sys
import os
import re


def _substitute_variables(f, out, items, values):
    while True:
        line = f.readline()
        regex_note = r"notebook\.set\(\"([A-Za-z_]*)\",.[0-9]*\.?[0-9]*\)"
        match = re.match(regex_note, line)
        if match is not None:
            regex_search = r",.[0-9]*\.?[0-9]*\)"
            line = re.sub(regex_search, ", " + match.group(1) + ")", line)
        out.write(line)
        if "End of NoteBook variables section" in line:
            return


def changeSALOMEfile(SALOME_PY_FILE, geometry, geometryvalues, meshname=None):
        """
        Take a SALOME file, which contains an absolute path and make it more general.
        Substitute numbers by variables and make file thus flexible.

        :param str SALOME_PY_FILE: name of template SALOME file that will be changed
        :param list geometry: list with name of geometrical entities
        :param dict geometryvalues: dictionary with geometrical values
        :param str meshname: usually name of study and '_mesh'

        .. note:: whenever something does not work, check the following points:

            #. is the error due to a missing declaration of geometryvalues? => most likely "sys.path.insert( 0, r'/path/to/file')" has been changed. then the list of declarations is not added! set the "sys.path.insert" to the aforementioned option.

        .. warning:: if you use a STEP file as input file, it must be located in the same directory as your meshfile and the study itself!
        """
        # take meshname from salome file, maybe better solution needed
        filename, ext = os.path.splitext(SALOME_PY_FILE)
        if ext != '.py':
            print("must be python file")
            sys.exit(0)
        if meshname is None:
            meshname = filename
        SALOME_OUT = SALOME_PY_FILE + '.tmp'
        stream = open(SALOME_PY_FILE, "r")
        out = open(SALOME_OUT, "w")

        while True:
            line = stream.readline()
            if not line:
                # EOF
                break

            # import os
            if line.strip() == "import sys":
                out.write("import os\n")
                out.write(line.strip() + " #\n")
                continue
            # get current path
            if "notebook =" in line.strip() and "#" not in line.strip():
                out.write("pwd = os.getcwd()\n")
                out.write("print(\"Current path is {}\".format(pwd))\n")
                out.write(line.strip() + " #\n")
                continue
            match = re.match("sys.path.insert\(.0,.r\'[/A-Za-z0-9]*\'\)", line.strip())
            # change system path
            # and write geometry plus meshname if not done previously
            if match is not None:
                out.write("sys.path.insert(0, pwd)\n")
                # read next line to write meshname
                line = stream.readline()
                # write meshname
                out.write("meshname = \'" + meshname + "\'\n")
                if geometry is not None:
                    out.write("\n### geometry parameters###\n")
                    for item in geometry:
                        out.write(item + " = " + str(geometryvalues[item]) + "\n")
                continue
            # write geometry parameters if file has been edited before
            if '### geometry parameters###' in line:
                out.write(line)
                for item in geometry:
                    # read line, do not do anything else
                    line = stream.readline()
                    out.write(item + " = " + str(geometryvalues[item]) + "\n")
                continue

            # write meshname
            if "meshname = " in line:
                out.write("meshname = \'" + meshname + "\'\n")
                continue

            if "Begin of NoteBook variables section" in line and geometry is not None:
                out.write(line)
                _substitute_variables(stream, out, geometry, geometryvalues)
                continue
            elif geometry is None:
                print("You have not specified a list of gemeotry objects!")
            # change path in ExportMED file
            if "ExportMED" in line.strip():
                line = re.sub(r'\(.r.[/[A-Za-z_0-9]*.med', '(pwd + \'/\' + str(meshname) + \'.med', line)
            if "ImportSTEP" in line.strip():
                line = re.sub(r'\(.[/[A-Za-z_0-9]*/(.*).step', r'''(pwd + '/' + '\1.step''', line)
                line = line.replace('\"', '\'')
            # always write the current line
            out.write(line)
        # move file
        os.rename(SALOME_PY_FILE + '.tmp', SALOME_PY_FILE)
