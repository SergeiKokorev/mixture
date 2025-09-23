import os
import sys


from PySide6.QtWidgets import QMainWindow, QWidget, QHBoxLayout
from PySide6.QtCore import Qt


DIR = os.path.abspath(os.pardir)
LIB = os.path.abspath(os.path.join("lib", "fluid_lib.json"))
print(LIB)
sys.path.append(DIR)


from fluid.properties import *
from fluid.datalib_handle import load_lib
from gui.widgets import *



class View(QMainWindow):
    
    def __init__(self, /, parent=None, flag=Qt.WindowType.Window):
        super().__init__(parent, flag)

        fluids = load_lib(LIB)
        self.fluid_list = FluidList(fluids)
        widget = QWidget()
        layout = QHBoxLayout()
        layout.addWidget(self.fluid_list)
        widget.setLayout(layout)


        self.setCentralWidget(widget)

