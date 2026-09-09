from VnE.Source.Extensions.ChemPack import attachCPprop
#asd
SELECTOR = True
if SELECTOR:
    from PySide6 import QtWidgets
    class Dialog(QtWidgets.QDialog):
        def __init__(self, *args, regex=None, select=None, model=None, attach=None, detach=None, **kwargs):
            from PySide6 import QtWidgets
            super().__init__(*args, **kwargs)
            self.regex = regex
            self.select = select
            self.mols = []
            self.curr_sys = None
            self.model = model
            self.attach = attach
            self.detach = detach

            self.main_layout = QtWidgets.QFormLayout(self)
            self.prop_field = QtWidgets.QLineEdit()
            self.prop_field.setText('label')
            self.prop_field.setPlaceholderText('label')
            self.regex_field = QtWidgets.QLineEdit()
            self.regex_field.setText(r'.*')
            self.regex_field.setPlaceholderText(r'C\d{2}')
            dash_label = QtWidgets.QLabel('-')
            self.regex_field2 = QtWidgets.QLineEdit()
            self.regex_field2.setPlaceholderText(r'C\d{2}')
            self.regex_field2.setText(r'.*')
            self.regex_frame = QtWidgets.QFrame()
            frame_layout = QtWidgets.QHBoxLayout()
            frame_layout.setContentsMargins(0, 0, 0, 0)
            self.regex_frame.setLayout(frame_layout)
            frame_layout.addWidget(self.regex_field)
            frame_layout.addWidget(dash_label)
            frame_layout.addWidget(self.regex_field2)

            self.dist_field = QtWidgets.QLineEdit()
            self.dist_field.setPlaceholderText(r'0.0')
            dash_label = QtWidgets.QLabel('-')
            self.dist_field2 = QtWidgets.QLineEdit()
            self.dist_field2.setPlaceholderText(r'-1.0')
            self.dist_frame = QtWidgets.QFrame()
            frame_layout = QtWidgets.QHBoxLayout()
            frame_layout.setContentsMargins(0, 0, 0, 0)
            self.dist_frame.setLayout(frame_layout)
            frame_layout.addWidget(self.dist_field)
            frame_layout.addWidget(dash_label)
            frame_layout.addWidget(self.dist_field2)

            self.radio_frame = QtWidgets.QFrame()
            self.radio_bond = QtWidgets.QRadioButton('Bond')
            self.radio_line = QtWidgets.QRadioButton('Line')
            frame_layout = QtWidgets.QHBoxLayout()
            frame_layout.setContentsMargins(0, 0, 0, 0)
            frame_layout.addWidget(self.radio_bond)
            frame_layout.addWidget(self.radio_line)
            self.radio_line.setAutoExclusive(True)
            self.radio_bond.setAutoExclusive(True)
            self.radio_bond.setChecked(True)
            self.radio_frame.setLayout(frame_layout)

            self.rad_field = QtWidgets.QLineEdit()
            self.rad_field.setPlaceholderText('0.05')
            self.freq_field = QtWidgets.QLineEdit()
            self.freq_field.setPlaceholderText('3')
            self.freq_field.setEnabled(False)
            self.hfreq_field = QtWidgets.QLineEdit()
            self.hfreq_field.setPlaceholderText('2')
            self.hfreq_field.setEnabled(False)

            self.type_combo = QtWidgets.QComboBox()
            self.type_combo.addItems(['Mask', 'List'])
            self.syst_combo = QtWidgets.QComboBox()

            self.select_button = QtWidgets.QPushButton('Select')
            self.delete_button = QtWidgets.QPushButton('Delete')
            self.modify_button = QtWidgets.QPushButton('Modify')
            self.create_button = QtWidgets.QPushButton('Create')

            self.button_frame = QtWidgets.QFrame()
            frame_layout = QtWidgets.QVBoxLayout()
            self.button_frame.setLayout(frame_layout)
            frame_layout.addWidget(self.select_button)
            frame_layout.addWidget(self.modify_button)
            frame_layout.addWidget(self.create_button)
            frame_layout.addWidget(self.delete_button)


            self.main_layout.addRow('System', self.syst_combo)
            self.main_layout.addRow('Property', self.prop_field)
            self.main_layout.addRow(self.type_combo, self.regex_frame)
            self.main_layout.addRow('Distance', self.dist_frame)
            self.main_layout.addRow(self.radio_frame)
            self.main_layout.addRow('Radius', self.rad_field)
            self.main_layout.addRow('Freq', self.freq_field)
            self.main_layout.addRow('HFreq', self.hfreq_field)
            self.main_layout.addRow(self.button_frame)

            from .Extensions.ChemPack import MOLECULE_SYSTEMS

            for mol_l in MOLECULE_SYSTEMS:
                if mol_l.isValid():
                    self.mols.append(mol_l)
                    self.syst_combo.addItem(mol_l.name)

            self.changeSystem(0)

            self.syst_combo.currentIndexChanged.connect(self.changeSystem)
            self.type_combo.currentIndexChanged.connect(self.typeChange)

            self.radio_bond.toggled.connect(self.setBond)
            self.radio_line.toggled.connect(self.setLine)

            self.select_button.pressed.connect(self.selectBonds)
            self.delete_button.pressed.connect(self.deleteBonds)
            self.create_button.pressed.connect(self.createBonds)
            self.modify_button.pressed.connect(self.modifyBonds)

        def changeSystem(self, ind):
            if ind != -1:
                try:
                    self.curr_sys = self.mols[ind]
                except IndexError:
                    return

        def typeChange(self, ind):
            if ind == 0:
                self.regex_field.setPlaceholderText(r'C\d{2}')
                self.regex_field2.setPlaceholderText(r'C\d{2}')
            if ind == 1:
                self.regex_field.setPlaceholderText(r'C1;C2;S5')
                self.regex_field2.setPlaceholderText(r'C1;C2;S5')
            self.regex_field.setText('')
            self.regex_field2.setText('')

        def setBond(self, b):
            if b:
                self.freq_field.setEnabled(False)
                self.hfreq_field.setEnabled(False)

        def setLine(self, b):
            if b:
                self.freq_field.setEnabled(True)
                self.hfreq_field.setEnabled(True)

        def getEnv(self):
            if self.curr_sys:
                atoms = self.curr_sys.children[0].children
                bonds_l = self.curr_sys.children[1]
            else:
                QtWidgets.QMessageBox.warning(self, 'Error', 'No molecule system selected')
                return
            prop = self.prop_field.text()
            l, regex = None, None
            if self.type_combo.currentIndex() == 1:
                l = (self.regex_field.text().split(';'), self.regex_field2.text().split(';'))
                if not (l[0] and l[1]):
                    l = None
            elif self.type_combo.currentIndex() == 0:
                regex = (self.regex_field.text(), self.regex_field2.text())
                if not (regex[0] and regex[1]):
                    regex = None
            try:
                dist = (float(self.dist_field.text()), float(self.dist_field2.text()))
            except ValueError:
                dist = None
            try:
                rad = float(self.rad_field.text())
            except ValueError:
                rad = None
            try:
                freq = float(self.freq_field.text())
            except ValueError:
                freq = None
            try:
                hfreq = float(self.hfreq_field.text())
            except ValueError:
                hfreq = None

            return atoms, bonds_l, prop, l, regex, dist, rad, freq, hfreq

        def selectBonds(self, select_bonds=True):
            atoms, bonds_l, prop, l, regex, dist, *_ = self.getEnv()
            if regex:
                sel1 = self.selectAtoms(atoms, prop, reg=regex[0])
                sel2 = self.selectAtoms(atoms, prop, reg=regex[1])
            if l:
                sel1 = self.selectAtoms(atoms, prop, reg=l[0])
                sel2 = self.selectAtoms(atoms, prop, reg=l[1])
            atoms2 = [x._atom for x in sel2]
            b_sel = []
            for p1 in sel1:
                a1 = p1._atom
                if a1 is None:
                    continue
                for b in a1.bonds():
                    a2 = b.get(a1)
                    if a2 in atoms2:
                        if b not in b_sel:
                            b_sel.append(b)
            if dist:
                import numpy as np
                for b in b_sel:
                    d = np.linalg.norm(b.parents()[0].coord - b.parents()[1].coord)
                    if dist[0] <= d <= dist[1]:
                        continue
                    else:
                        b_sel.remove(b)
            points = []
            for b in b_sel:
                points += b.point()
            if points and select_bonds:
                self.select(points)
                return
            return b_sel

        def deleteBonds(self):
            bonds = self.selectBonds(False)
            for b in bonds:
                p_list = b.point()[0].parent
                ind = self.model.index(by_point=p_list)
                self.model.removeRow(ind.row(), ind.parent())
                b.destroy()

        def createBonds(self):
            atoms, bonds_l, prop, l, regex, dist, *_ = self.getEnv()
            if regex:
                sel1 = self.selectAtoms(atoms, prop, reg=regex[0])
                sel2 = self.selectAtoms(atoms, prop, reg=regex[1])
            if l:
                sel1 = self.selectAtoms(atoms, prop, reg=l[0])
                sel2 = self.selectAtoms(atoms, prop, reg=l[1])
            atoms2 = [x._atom for x in sel2]
            if sel1 and sel2:
                from .Extensions.ChemPack.MoleculeClass import Bond
            if dist:
                import numpy as np
                for p1 in sel1:
                    a1 = p1._atom
                    if a1 is None:
                        continue
                    for p2 in sel2:
                        a2 = p2._atom
                        if a2 is None or a1 is a2:
                            continue
                        d = np.linalg.norm(a1.coord - a2.coord)
                        if dist[0] <= d <= dist[1]:
                            bond = Bond(a1, a2)
                            r = a1.addBond(bond), a2.addBond(bond)
                            if all(r):
                                bond_l = PointsList(parent=bonds_l,
                                                    name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}',
                                                    rad=bonds_l, freq=bonds_l, hfreq=bonds_l)
                                b1 = Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                                           name=bond.parents()[0].point(),
                                           rad=bond_l,
                                           parent=bond_l,
                                           freq=bond_l, hfreq=bond_l)
                                b2 = Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                                           name=bond.parents()[1].point(),
                                           rad=bond_l,
                                           parent=bond_l,
                                           freq=bond_l, hfreq=bond_l)
                                bond.assignPoint((b1, b2))
            else:
                for p1 in sel1:
                    a1 = p1._atom
                    if a1 is None:
                        continue
                    for p2 in sel2:
                        a2 = p2._atom
                        if a2 is None or a1 is a2:
                            continue
                        bond = Bond(a1, a2)
                        r = a1.addBond(bond), a2.addBond(bond)
                        if all(r):
                            bond_l = PointsList(parent=bonds_l,
                                                name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}',
                                                rad=bonds_l, freq=bonds_l, hfreq=bonds_l)
                            b1 = Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                                       name=bond.parents()[0].point(),
                                       rad=bond_l,
                                       parent=bond_l,
                                       freq=bond_l, hfreq=bond_l)
                            b2 = Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                                       name=bond.parents()[1].point(),
                                       rad=bond_l,
                                       parent=bond_l,
                                       freq=bond_l, hfreq=bond_l)
                            bond.assignPoint((b1, b2))
            self.model.update()

        def modifyBonds(self):
            atoms, bonds_l, prop, l, regex, dist, rad, freq, hfreq = self.getEnv()

            bonds = self.selectBonds(False)
            for b in bonds:
                p_list = b.point()[0].parent
                if self.radio_bond.isChecked():
                    for o in p_list.observers:
                        if o.NAME != 'Bond':
                            self.detach(p_list, o.NAME)
                    self.attach(p_list, 'Bond')
                elif self.radio_line.isChecked():
                    for o in p_list.observers:
                        if o.NAME != 'Line':
                            self.detach(p_list, o.NAME)
                    self.attach(p_list, 'Line')
                    if freq:
                        p_list.freq = freq
                    else:
                        p_list.freq = p_list.parent
                        b.point()[0].freq = p_list
                        b.point()[1].freq = p_list
                    if hfreq:
                        p_list.hfreq = hfreq
                    else:
                        p_list.hfreq = p_list.parent
                        b.point()[0].hfreq = p_list
                        b.point()[1].hfreq = p_list
                if rad:
                    p_list.rad = rad
                else:
                    p_list.rad = p_list.parent
                    b.point()[0].rad = p_list
                    b.point()[1].rad = p_list

        def selectAtoms(self, atoms, property, l=None, reg=None):
            sel = []
            if reg:
                for p in atoms:
                    val = getattr(p, property)
                    if val:
                        val = str(val)
                        if self.regex(reg, val):
                            sel.append(p)
            elif l:
                for p in atoms:
                    val = getattr(p, property)
                    if val:
                        val = str(val)
                        if val in l:
                            sel.append(p)
            return sel

    DIALOG = Dialog(self.parent, regex=regex, select=select, model=self.model, attach=attach, detach=detach)
    DIALOG.show()
