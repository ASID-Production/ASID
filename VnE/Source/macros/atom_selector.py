SELECTOR = True
if SELECTOR:
    from PySide6 import QtWidgets
    class Dialog(QtWidgets.QDialog):
        def __init__(self, *args, regex=None, select=None, **kwargs):
            from PySide6 import QtWidgets
            super().__init__(*args, **kwargs)
            self.regex = regex
            self.select = select
            self.mols = []
            self.curr_sys = None

            self.main_layout = QtWidgets.QFormLayout(self)
            self.prop_field = QtWidgets.QLineEdit()
            self.prop_field.setPlaceholderText('label')
            self.regex_field = QtWidgets.QLineEdit()
            self.regex_field.setPlaceholderText(r'C\d{2}')

            self.type_combo = QtWidgets.QComboBox()
            self.type_combo.addItems(['Mask', 'List'])
            self.syst_combo = QtWidgets.QComboBox()

            self.select_button = QtWidgets.QPushButton('Select')

            self.main_layout.addRow('System', self.syst_combo)
            self.main_layout.addRow('Property', self.prop_field)
            self.main_layout.addRow(self.type_combo, self.regex_field)
            self.main_layout.addRow(self.select_button)

            from .Extensions.ChemPack import MOLECULE_SYSTEMS

            for mol_l in MOLECULE_SYSTEMS:
                if mol_l.isValid():
                    self.mols.append(mol_l)
                    self.syst_combo.addItem(mol_l.name)

            self.changeSystem(0)

            self.syst_combo.currentIndexChanged.connect(self.changeSystem)
            self.type_combo.currentIndexChanged.connect(self.typeChange)

            self.select_button.pressed.connect(self.apply)

        def changeSystem(self, ind):
            if ind != -1:
                try:
                    self.curr_sys = self.mols[ind]
                except IndexError:
                    return

        def typeChange(self, ind):
            if ind == 0:
                self.regex_field.setPlaceholderText(r'C\d{2}')
            if ind == 1:
                self.regex_field.setPlaceholderText(r'C1;C2;S5')
            self.regex_field.setText('')

        def apply(self):
            if self.curr_sys:
                atoms = self.curr_sys.children[0].children
            else:
                QtWidgets.QMessageBox.warning(self, 'Error', 'No molecule system selected')
                return
            prop = self.prop_field.text()
            if self.type_combo.currentIndex() == 1:
                l = self.regex_field.text().split(';')
                self.execute(atoms, prop, l=l)
            elif self.type_combo.currentIndex() == 0:
                self.execute(atoms, prop, l=None, reg=self.regex_field.text())

        def execute(self, atoms, property, l=None, reg=None):
            if reg:
                sel = []
                for p in atoms:
                    val = getattr(p, property)
                    if val:
                        val = str(val)
                        if self.regex(reg, val):
                            sel.append(p)
                self.select(sel)
            elif l:
                sel = []
                for p in atoms:
                    val = getattr(p, property)
                    if val:
                        val = str(val)
                        if val in l:
                            sel.append(p)
                self.select(sel)

    DIALOG = Dialog(self.parent, regex=regex, select=select)
    DIALOG.show()
