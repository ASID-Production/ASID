import os.path

DIALOG = None

def execute():
    from PySide6 import QtWidgets

    class PathDialog(QtWidgets.QDialog):

        def __init__(self):
            QtWidgets.QDialog.__init__(self)
            self.line_edit = QtWidgets.QLineEdit(parent=self)
            layout = QtWidgets.QVBoxLayout()
            layout.addWidget(self.line_edit)
            self.setLayout(layout)
            self.line_edit.returnPressed.connect(self.accept)
            self.accepted.connect(self.foo)
            self.path = None

        def foo(self):
            import os
            text = self.line_edit.text()
            self.line_edit.clear()
            if os.path.exists(text):
                self.path = text

        def show(self):
            self.path = None
            QtWidgets.QDialog.show(self)

    global DIALOG
    if DIALOG is None:
        DIALOG = PathDialog()
        func = lambda: foo(DIALOG.path)
        DIALOG.accepted.connect(func)
    DIALOG.show()

    def foo(path):
        if path is not None:
            from .ChemPackSource.parsers import PARSER
            import numpy as np
            cpProp_path = os.path.normpath(os.path.join(path, 'aim/CPprop.txt'))
            Cps_sys, _ = PARSER.parsCpProp(cpProp_path, False, None)

            pdb_path = os.path.normpath(os.path.join(path, f'aim/CPs.pdb'))
            CPs_pdb_sys, _ = PARSER.parsPdb(pdb_path, False, None)

            pdb_path = os.path.normpath(os.path.join(path, f'xtbopt.pdb'))
            mol_sys, _ = PARSER.parsPdb(pdb_path, True, None)

            cps = Cps_sys.children[0].children.copy()
            i = 0
            i2 = 0
            for cp in CPs_pdb_sys.children[0].children:
                i += 1
                i2 = 0
                for cp2 in cps:
                    i2 += 1
                    if np.all(np.abs(cp.coord - cp2.coord) < 0.01):
                        cp.data_cp = cp2
                        cps.remove(cp2)
                        break
            a = 0


def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Pars CPs')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions