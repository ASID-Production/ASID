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

import debug

TREE_MODEL = None

RES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HOH']
POS_RES = ['LYS', 'ARG']
SEMI_RES = ['HIS']
NEG_RES = ['ASP', 'GLU']

DIALOG = None


def execute():
    from PySide6 import QtWidgets
    from .ChemPack import MOLECULE_SYSTEMS

    class ResDialog(QtWidgets.QDialog):
        def __init__(self):
            QtWidgets.QDialog.__init__(self)

            button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
            button_box.accepted.connect(self.execute)
            button_box.rejected.connect(self.reject)

            res_label = QtWidgets.QLabel('Res id')
            self.res_line = QtWidgets.QLineEdit()
            range_label = QtWidgets.QLabel('Range')
            self.range_line = QtWidgets.QLineEdit()
            chain_label = QtWidgets.QLabel('Target chain')
            self.chain_line = QtWidgets.QLineEdit()

            form_frame = QtWidgets.QFrame()

            main_layout = QtWidgets.QVBoxLayout()
            form_layout = QtWidgets.QFormLayout()

            form_layout.addRow(res_label, self.res_line)
            form_layout.addRow(range_label, self.range_line)
            form_layout.addRow(chain_label, self.chain_line)

            form_frame.setLayout(form_layout)

            main_layout.addWidget(form_frame)
            main_layout.addWidget(button_box)

            self.setLayout(main_layout)

        def execute(self):
            from .ChemPackSource import MoleculeClass
            try:
                range = float(self.range_line.text())
                id = int(self.res_line.text())
                chain = self.chain_line.text()
            except ValueError:
                return
            mol_sys = None
            for mol_l in MOLECULE_SYSTEMS:
                if mol_l.pick == 1 and mol_l.isValid():
                    mol_sys = MOLECULE_SYSTEMS[mol_l]
                    break
            if mol_sys is None:
                return
            target = []
            other = []

            def rec(mol):
                for atom in mol.children:
                    if not isinstance(atom, MoleculeClass.Atom):
                        rec(atom)
                    else:
                        if atom.pdb_res_seq == id:
                            if chain:
                                if atom.pdb_chain == chain:
                                    if atom not in target:
                                        target.append(atom)
                            elif atom not in target:
                                target.append(atom)
                        elif atom not in other:
                            other.append(atom)

            rec(mol_sys)

            sel = []

            from .ChemPackSource.contacts import dist

            def rec_add(atom, prev_atom=None, mem=None):
                if mem is None:
                    mem = {}
                if prev_atom is None:
                    prev_atom = atom
                if atom in mem:
                    return
                if atom.pdb_res_seq == prev_atom.pdb_res_seq:
                    if atom not in sel:
                        sel.append(atom)
                    else:
                        return
                else:
                    return
                mem[atom] = True
                for bond in atom.bonds():
                    atom2 = bond.parents()[bond.parents().index(atom)-1]
                    rec_add(atom2, atom, mem)

            for tar_atom in target:
                for other_atom in other:
                    if dist(tar_atom, other_atom) <= range:
                        rec_add(other_atom)
            atoms = sel + target
            molecule = MoleculeClass.Molecule()
            for atom in atoms:
                molecule.addChild(atom)
            if not molecule.children:
                return
            from .ChemPackSource import save_file
            import os

            filters = f'*.{" *.".join(save_file.SAVE_FILE.getFormats())}'
            filename = QtWidgets.QFileDialog.getSaveFileName(filter=filters)
            if type(filename) == tuple:
                filename = filename[0]

            if filename != '':
                format = os.path.basename(filename).split('.')[1]
                save_file.SAVE_FILE.save(molecule, filename, format)

    global DIALOG
    DIALOG = ResDialog()
    DIALOG.show()


def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Cut prot in range')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions