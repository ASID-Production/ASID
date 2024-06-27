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


from .ChemPackSource.ui.select_mol_dialog import SelectMolDialog
import os.path as opath

TREE_MODEL = None
DIALOG = None


def execute():
    from .ChemPack import MOLECULE_SYSTEMS
    import rdkit
    from rdkit.Chem import Descriptors3D
    from PySide6 import QtWidgets

    DESCRIPTORS = Descriptors3D.descList

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

    def process(molsys):
        from .ChemPackSource.save_file import SAVE_FILE
        list_obj, molsys = molsys
        name = f'{opath.dirname(__file__)}/../../temp/{id(molsys)}.xyz'
        SAVE_FILE.save(molsys, name, 'xyz')
        rdkit_mol = rdkit.Chem.MolFromXYZFile(name)
        params = Descriptors3D.CalcMolDescriptors3D(rdkit_mol)
        text = ['For detailed description check: https://www.rdkit.org/new_docs/source/rdkit.Chem.Descriptors3D.html\n']
        a = [f'{x}: {y}' for x, y in params.items()]
        text = text + a
        dialog = TextViewer('\n'.join(text))
        global DIALOG
        DIALOG = dialog
        dialog.show()
        return

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()
    pass


def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('RDKit_descriptors.py')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions