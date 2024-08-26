# Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
# This file is part of ASID - Atomistic Simulation Instruments and Database
# For more information see <https://github.com/ASID-Production/ASID>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ******************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# ******************************************************************************************

from PySide6 import QtWidgets
import numpy as np

import debug

DIALOG = None


class SymOpDialog(QtWidgets.QDialog):

    def __init__(self, tree_model):

        from ..ChemPack import MOLECULE_SYSTEMS

        QtWidgets.QDialog.__init__(self)
        self.setWindowTitle('Symmetry operations')
        self.tree_model = tree_model
        self.main_layout = QtWidgets.QVBoxLayout(self)

        self.cb_frame = QtWidgets.QFrame(self)
        self.combo_box_layout = QtWidgets.QFormLayout()
        self.cb_label = QtWidgets.QLabel('System:')
        self.combo_box = QtWidgets.QComboBox(self.cb_frame)
        self.combo_box_layout.addRow(self.cb_label, self.combo_box)
        self.cb_frame.setLayout(self.combo_box_layout)
        self.combo_box.currentIndexChanged.connect(self.changeSystem)

        self.main_layout.addWidget(self.cb_frame)

        self.check_box = QtWidgets.QCheckBox('Center cell', self)
        self.main_layout.addWidget(self.check_box)

        self.sop_frame = QtWidgets.QFrame(self)
        self.main_layout.addWidget(self.sop_frame)

        self.reset_button = QtWidgets.QPushButton('Reset')
        self.reset_button.pressed.connect(self.resetSystem)
        self.main_layout.addWidget(self.reset_button)

        self.curr_sys = None
        self.mols = []
        self.original_systems = []

        for mol_l in MOLECULE_SYSTEMS:
            if mol_l.isValid():
                self.mols.append((mol_l, MOLECULE_SYSTEMS[mol_l]))
                self.combo_box.addItem(mol_l.name)
                self.original_systems.append((mol_l.copy(obs_clone=False), MOLECULE_SYSTEMS[mol_l]))

        self.changeSystem(0)

    def updateSymOps(self):
        if self.curr_sys is None:
            sop_layout = QtWidgets.QFormLayout(self.sop_frame)
            self.sop_frame.setLayout(sop_layout)
            return
        else:
            sop_layout = QtWidgets.QFormLayout(self.sop_frame)
            atom = None
            try:
                atom = self.curr_sys[1].children[0].children[0]
            except IndexError:
                pass
            widget = self.sop_frame.layout().takeAt(0)
            while widget:
                widget.widget().setParent(None)
                widget = self.sop_frame.layout().takeAt(0)
            if atom is not None:
                for i, symop in enumerate(atom.cif_sym_codes):
                    button = QtWidgets.QPushButton('Apply')
                    func = self.applySymOpFunc(i)
                    button.pressed.connect(func)
                    self.sop_frame.layout().addRow(symop[1], button)

    def changeSystem(self, ind):
        if ind != -1:
            try:
                self.curr_sys = self.mols[ind]
            except IndexError:
                return
            self.updateSymOps()
            self.update()

    def resetSystem(self):
        from ..ChemPack import TREE_MODEL, PARSER, loadMolSys, UNIFORM_MODEL
        i = self.mols.index(self.curr_sys)
        orig_sys = list(self.original_systems[i])
        orig_sys[0] = orig_sys[0].copy(obs_clone=False)
        index = TREE_MODEL.index(0, 0, by_point=self.curr_sys[0])
        parent = index.parent()
        TREE_MODEL.removeRow(index.row(), parent=parent)
        i = self.mols.index(self.curr_sys)
        self.mols.pop(i)
        self.mols.insert(i, orig_sys)
        self.curr_sys = self.mols[i]
        TREE_MODEL.getRoot().addChild(orig_sys[0])
        TREE_MODEL.insertRow(0, parent=parent)
        atoms_index = TREE_MODEL.index(0, 0, by_point=orig_sys[0].children[0])
        bonds_index = TREE_MODEL.index(0, 0, by_point=orig_sys[0].children[1])
        TREE_MODEL.attachObserver(atoms_index, 'Sphere')
        TREE_MODEL.attachObserver(bonds_index, 'Bond')
        coords = np.array([a.coord for a in orig_sys[0].children[0].children])
        cent = np.sum(coords, axis=0) / len(coords)
        if UNIFORM_MODEL is not None:
            root = UNIFORM_MODEL.root()
            root.rotation_point = cent.copy()

    def applySymOpFunc(self, ind):
        from . import MoleculeClass
        from ..ChemPack import PALETTE, TREE_MODEL, PARSER, loadMolSys
        import numpy as np

        def fracToDec(a, b, c, al, be, ga, coords):
            al = (al/180)*np.pi
            be = (be/180)*np.pi
            ga = (ga/180)*np.pi

            sin = np.sin
            cos = np.cos
            cot = lambda x: np.tan(x)**-1
            csc = lambda x: np.sin(x)**-1

            mat = np.array([[a*sin(be)*np.sqrt(1-(cot(al)*cot(be) - csc(al)*csc(be)*cos(ga))**2), 0, 0],
                            [a*csc(al)*cos(ga) - a*cot(al)*cos(be), b*sin(al), 0],
                            [a*cos(be), b*cos(al), c]])
            mat = mat.transpose()
            for i in range(len(coords)):
                coords[i] = (coords[i] @ mat).astype(dtype=np.float32)
            return coords

        def applySymOps():
            import cpplib

            try:
                cell_list_copy = self.curr_sys[0].children[2].copy(obs_clone=False)
            except IndexError:
                cell_list_copy = None

            atoms = self.curr_sys[1].children[0].children
            sym_codes = atoms[0].cif_sym_codes
            sym_code = sym_codes[ind][1]
            data = [[x.atom_type, *list(x.cif_frac_coords)] for x in atoms]
            flag = 0b0 | (int(self.check_box.isChecked()) << 1)
            new_atoms = cpplib.GenSymm(data, flag, [sym_code])
            new_coords = [x[1:] for x in new_atoms]
            cell = [atoms[0].cif_cell_a,
                    atoms[0].cif_cell_b,
                    atoms[0].cif_cell_c,
                    atoms[0].cif_cell_al,
                    atoms[0].cif_cell_be,
                    atoms[0].cif_cell_ga]
            args = []
            args += cell
            args.append(new_coords)
            new_dec_coords = fracToDec(*args)
            new_mol_sys = MoleculeClass.MoleculeSystem()
            new_mol_sys.name = self.curr_sys[1].name
            new_mol_sys.file_name = self.curr_sys[1].file_name
            mol = MoleculeClass.Molecule(new_mol_sys)
            for i, atom in enumerate(new_atoms):
                cif_data = {'cif_sym_codes': sym_codes,
                            'cif_cell_a': cell[0],
                            'cif_cell_b': cell[1],
                            'cif_cell_c': cell[2],
                            'cif_cell_al': cell[3],
                            'cif_cell_be': cell[4],
                            'cif_cell_ga': cell[5],
                            'cif_frac_coords': np.array(atom[1:], dtype=np.float32)}
                coord = np.array(new_dec_coords[i], dtype=np.float32)
                if i < len(atoms):
                    atom = MoleculeClass.Atom(coord.copy(), atom[0], parent=mol, name=atoms[i].name, **cif_data)
                else:
                    atom = MoleculeClass.Atom(coord.copy(), atom[0], parent=mol, name=PALETTE.getName(atom[0]), **cif_data)
            args = PARSER.parsMolSys(new_mol_sys, True, TREE_MODEL.getRoot())
            loadMolSys(args[0], args[1])
            index = TREE_MODEL.index(0, 0, by_point=self.curr_sys[0])
            TREE_MODEL.removeRow(index.row(), parent=index.parent())
            i = self.mols.index(self.curr_sys)
            self.mols.pop(i)
            self.mols.insert(i, (args[1][0], args[0]))
            self.curr_sys = self.mols[i]

            if cell_list_copy:
                self.curr_sys[0].addChild(cell_list_copy)
                index = TREE_MODEL.index(0, 0, by_point=self.curr_sys[0])
                TREE_MODEL.insertRow(TREE_MODEL.rowCount(parent=index), parent=index)

        return applySymOps


def show():
    from ..ChemPack import TREE_MODEL
    dialog = SymOpDialog(TREE_MODEL)
    global DIALOG
    DIALOG = dialog
    dialog.show()