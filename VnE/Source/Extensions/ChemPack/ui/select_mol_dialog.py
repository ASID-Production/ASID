from PySide6 import QtWidgets


class SelectMolDialog(QtWidgets.QDialog):

    def __init__(self, MOLECULE_SYSTEMS, call_method):
        self.call_method = call_method
        self.mol_sys = MOLECULE_SYSTEMS
        QtWidgets.QDialog.__init__(self)
        self.setWindowTitle('Select molecule')
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.cb_frame = QtWidgets.QFrame(self)
        self.combo_box_layout = QtWidgets.QFormLayout()
        self.cb_label = QtWidgets.QLabel('System:')
        self.combo_box = QtWidgets.QComboBox(self.cb_frame)
        self.combo_box_layout.addRow(self.cb_label, self.combo_box)
        self.cb_frame.setLayout(self.combo_box_layout)

        self.main_layout.addWidget(self.cb_frame)

        buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.main_layout.addWidget(buttonBox)

        self.curr_sys = None
        self.mols = []
        for mol_l in self.mol_sys:
            if mol_l.isValid():
                self.mols.append((mol_l, MOLECULE_SYSTEMS[mol_l]))
                self.combo_box.addItem(mol_l.name)

        self.changeSystem(0)

        self.combo_box.currentIndexChanged.connect(self.changeSystem)

        self.accepted.connect(self.apply)

    def changeSystem(self, ind):
        if ind != -1:
            try:
                self.curr_sys = self.mols[ind]
            except IndexError:
                return

    def apply(self):
        if self.curr_sys is not None:
            self.call_method(self.curr_sys)
        else:
            return
