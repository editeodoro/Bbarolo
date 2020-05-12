#!/usr/bin/env python
# ---------------------------------------------------------------------
#
#-----------------------------------------------------------------------

import sys, argparse

def read_command_line(argv):
    """Read options from the command line."""
    par = argparse.ArgumentParser(prog='BBaroloGUI', epilog='BBaroloGUI 1.0')
    par.add_argument("--command-line", default=False, action="store_true", help="Run in the command line.")
    par.add_argument("-style", type=str, help="Style of windows.")
    return par.parse_args()
    

if __name__ == '__main__':
    
    args = read_command_line(sys.argv)
    
    if args.command_line:
        # Using as a command line utility
        ...
    else:
        # Using the QT GUI
        from bbarolowindow import MainWindow
        from PyQt5.QtWidgets import QApplication
        app = QApplication(sys.argv)
        frame = MainWindow(args=args)
        frame.show()
        app.exec_()

