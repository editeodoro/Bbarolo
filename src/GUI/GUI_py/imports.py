__proname__ = "BBaroloGUI"
__version__ = "0.1.0"
__author__  = "Enrico Di Teodoro (enrico.diteodoro@gmail.com)"

import os, sys

try: 
    from PySide6.QtWidgets import *
    from PySide6.QtCore import *
    from PySide6.QtGui import *
#    from PySide6.QtCharts import QLineSeries, QChart, QChartView
except ImportError:
    try:
        from PyQt6.QtWidgets import *
        from PyQt6.QtCore import *
        from PyQt6.QtGui import *
        from PyQt6.QtCore import pyqtSignal as Signal, pyqtSlot as Slot
    except ImportError:
        raise ImportError(f" ERROR : {__proname__} requires PySide6 or PyQt6") 