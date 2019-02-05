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

import subprocess
import os


def callSalome(pythonfile):
    """
    calls salome for pythonfile and kills afterwards the opened session

    .. todo:: find solution for killSalomeWithPort.py, i.e. add it to the path automatically or add it to the sources of this package

    .. warning:: If the error "Can't find a free port to launch omniNames.    Try to kill the running servers and then launch SALOME again." pops up, use the trick described on http://academic.bancey.com/successfully-clearing-ports-in-salome-code-aster/ . Briefly this means to run:

        .. code-block:: bash

                rm /tmp/.salome_PortManager.*
                rm /tmp/.omniORB_*

    """
    subprocess.call('salome -t --ns-port-log=$PWD/portlog.txt ' + pythonfile, shell=True)
    # requires SALOME-8.3.0-UB16.04/BINARIES-UB16.04/KERNEL/bin/salome/ in path
    subprocess.call('tail -1 portlog.txt | xargs killSalomeWithPort.py', shell=True)
    os.remove('portlog.txt')
    return
