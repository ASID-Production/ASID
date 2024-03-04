from Source import MainWindow
import sys
import os

module_dir = os.path.normpath(os.path.join(os.path.split(__file__)[0], '../module'))
sys.path.append(module_dir)

MainWindow.show()
