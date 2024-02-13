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


TREE_MODEL = None

RES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HOH']
POS_RES = ['LYS', 'ARG']
SEMI_RES = ['HIS']
NEG_RES = ['ASP', 'GLU']

DIALOG = None
LAYOUT = None


def execute():
    from .ChemPack import MOLECULE_SYSTEMS
    from PySide6 import QtWidgets
    mol_sys = None
    for mol_l in MOLECULE_SYSTEMS:
        if mol_l.pick == 1.0 and mol_l.isValid():
            mol_sys = MOLECULE_SYSTEMS[mol_l]
            break
    if mol_sys is not None:
        text = []
        charge = 0
        mem = []
        pos_res = []
        neg_res = []
        try:
            a = mol_sys.children[0].children[0].pdb_res_seq
        except AttributeError:
            return
        for atom in mol_sys.children[0].children:
            if atom.pdb_res_seq in mem:
                continue
            else:
                if atom.pdb_res_name not in RES:
                    line = f'WARNING: Unknown res name {atom.pdb_res_name}, res seq: {atom.pdb_res_seq}'
                    text.append(line)
                elif atom.pdb_res_name in SEMI_RES:
                    line = f'WARNING: Unknown charge for res name {atom.pdb_res_name}, res seq: {atom.pdb_res_seq}'
                    text.append(line)
                elif atom.pdb_res_name in POS_RES:
                    charge += 1
                    pos_res.append(atom)
                elif atom.pdb_res_name in NEG_RES:
                    charge -= 1
                    neg_res.append(atom)
                mem.append(atom.pdb_res_seq)
        pos_line = f'Positive res: {", ".join([f"{x.pdb_res_seq}-{x.pdb_res_name}" for x in pos_res])}'
        text.append(pos_line)
        neg_line = f'Negative res: {", ".join([f"{x.pdb_res_seq}-{x.pdb_res_name}" for x in neg_res])}'
        text.append(neg_line)
        line = f'FINAL: Charge == {charge}'
        text.append(line)

        class TextViewer(QtWidgets.QWidget):
            def __init__(self, text):
                QtWidgets.QWidget.__init__(self)

                self.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred,
                                                         QtWidgets.QSizePolicy.Policy.Preferred))

                layout = QtWidgets.QVBoxLayout()
                text_edit = QtWidgets.QTextEdit(text=text)
                text_edit.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                                              QtWidgets.QSizePolicy.Policy.Expanding))
                text_edit.setReadOnly(True)
                layout.addWidget(text_edit)
                self.setLayout(layout)

        dialog = TextViewer('\n'.join(text))
        global DIALOG
        DIALOG = dialog
        dialog.show()


def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Calc PDB charge')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions
