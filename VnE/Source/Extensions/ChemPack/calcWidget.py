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
import numpy as np

from ..ChemPack import MAIN_WIDGET, MOLECULE_SYSTEMS, TREE_MODEL
from . import contacts
from PySide6.QtWidgets import QLabel, QFrame, QSizePolicy
from PySide6.QtCore import Qt, QItemSelectionModel


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
        coords = []
        for atms in list(self.sel.values()):
            if atms.isValid():
                atoms.append(atms)
                coords.append(atms.coord)
            else:
                self.sel.pop(atms)
        if len(atoms) > 6:
            atoms = [*atoms[:3], *atoms[-3:]]
        for i, a in enumerate(atoms):
            b = [all(a.coord == x) for x in coords[i+1:]]
            if any(b):
                ind = [x + i + 1 for x in range(len(b)) if b[x]]
                for i1 in ind:
                    try:
                        coords.pop(i1)
                        atoms.pop(i1)
                    except IndexError:
                        pass
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
        point = index.internalPointer()
        if isinstance(point.coord, np.ndarray) and point.coord.shape == (3,):
            self.sel[point] = point
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