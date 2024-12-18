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
import numpy as np

def execute():
    from ..ChemPack import MOLECULE_SYSTEMS
    from PySide6.QtWidgets import QFileDialog
    import cpplib
    import csv

    mol_sys = None

    for mol_l in MOLECULE_SYSTEMS:
        if mol_l.isValid() and mol_l.pick:
            mol_sys = MOLECULE_SYSTEMS[mol_l]
            break

    if mol_sys is not None:
        filepath = QFileDialog.getSaveFileName(filter='*.csv')
        if not filepath[0]:
            return
        atoms = mol_sys.children[0].children
        arg = [(x.atom_type, *x.coord) for x in atoms]
        ret_dict = cpplib.FindDAT_WC(arg)
        dist = ret_dict['bonds']
        angles = ret_dict['angles']
        tors = ret_dict['tors']
        with open(filepath[0], 'w', newline='') as out:
            writer = csv.writer(out, delimiter=';')
            lines = []
            for pair in dist:
                line = [f'{atoms[pair[0]].name}--{atoms[pair[1]].name}', str(round(pair[2], 4))]
                lines.append(line)
            writer.writerow(['Pair', 'Dist'])
            writer.writerows(lines)
            writer.writerow([''])

            lines = []
            for pair in angles:
                line = [f'{atoms[pair[0]].name}--{atoms[pair[1]].name}--{atoms[pair[2]].name}', str(round(pair[3]/np.pi*180, 2))]
                lines.append(line)
            writer.writerow(['Atoms', 'Angle'])
            writer.writerows(lines)
            writer.writerow([''])

            lines = []
            for pair in tors:
                line = [f'{atoms[pair[0]].name}--{atoms[pair[1]].name}--{atoms[pair[2]].name}--{atoms[pair[3]].name}', str(round(pair[4]/np.pi*180, 2))]
                lines.append(line)
            writer.writerow(['Atoms', 'Angle'])
            writer.writerows(lines)
