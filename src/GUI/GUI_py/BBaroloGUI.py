#!/usr/bin/env python3
"""
Main file for BBaroloGUI
"""

########################################################################
# Copyright (C) 2023 Enrico Di Teodoro
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

#import argparse
from imports import *

from bbarolowindow import MainWindow

def read_command_line(argv):
    """Read options from the command line."""
    ...
    #par = argparse.ArgumentParser(prog=__proname__, epilog=f'{__proname__} {__version__}')
    #par.add_argument("-f", "--rcfile", type=str, help="File with rotation curve.")
    #par.add_argument("-c", "--columns", type=int, nargs='+', help="Column numbers for radius and rotation Velocity")
    
    #return par.parse_args()
    

if __name__ == "__main__":

    args = read_command_line(sys.argv)
    
    app = QApplication(sys.argv)
    frame = MainWindow(args=args)
    frame.show()
    app.exec()
 
