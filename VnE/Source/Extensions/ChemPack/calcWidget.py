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

from ..ChemPack import MAIN_WIDGET, MOLECULE_SYSTEMS, TREE_MODEL
from . import contacts
from PySide6.QtWidgets import QLabel, QFrame, QSizePolicy
from PySide6.QtCore import Qt


class CalcWidget(QLabel):

    def __init__(self, *args, **kwargs):
        QLabel.__init__(self, *args, **kwargs)
        self.sel = {}
        TREE_MODEL.item_selected.connect(self.add)
        TREE_MODEL.item_deselected.connect(self.remove)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

    def calc(self):
        atoms = []
        for atms in self.sel.values():
            atoms += atms
        if len(atoms) > 4:
            angle = contacts.angle(atoms, False)
            line = f'{atoms[0].name}--{atoms[1].name}--{atoms[2].name}  {atoms[-1].name}--{atoms[-2].name}--{atoms[-3].name}: {angle: .1f}'
        elif len(atoms) == 4:
            angle = contacts.angle(atoms, False)
            line = f'{atoms[0].name}--{atoms[1].name}--{atoms[2].name}--{atoms[3].name}: {angle: .1f}'
        elif len(atoms) == 3:
            angle = contacts.angle(atoms, False)
            line = f'{atoms[0].name}--{atoms[1].name}--{atoms[2].name}: {angle: .1f}'
        elif len(atoms) == 2:
            dist = contacts.dist(*atoms)
            line = f'{atoms[0].name}--{atoms[1].name}: {dist: .3f}'
        else:
            line = [str(x.name) for x in atoms]
            line = '--'.join(line)
        self.setText(line)

    def add(self, index):
        atoms = []
        point = index.internalPointer()
        for mol_sys in MOLECULE_SYSTEMS.values():
            atoms += mol_sys.findProp('_point', point)
        atoms_with_coords = []
        for atom in atoms:
            if atom.coord is not None:
                atoms_with_coords.append(atom)
        if atoms_with_coords:
            self.sel[point] = atoms_with_coords
            self.calc()

    def remove(self, index):
        point = index.internalPointer()
        res = self.sel.pop(point, None)
        if res is not None:
            self.calc()


opengl_frame = MAIN_WIDGET.findChild(QFrame, 'opengl_frame')
label = CalcWidget()
opengl_frame.layout().addWidget(label)
a = 0