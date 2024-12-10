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

DIALOG = None


def execute(mol=None):
    from .ui import select_mol_dialog
    from ..ChemPack import MOLECULE_SYSTEMS
    global DIALOG
    if mol is None:
        DIALOG = select_mol_dialog.SelectMolDialog(MOLECULE_SYSTEMS, attachCPprop)
        DIALOG.show()
    else:
        attachCPprop(mol)


def attachCPprop(molsys):
    from PySide6.QtWidgets import QFileDialog
    from .parsers import PARSER
    molsys = molsys[1]
    file_path = QFileDialog.getOpenFileName(filter='CPprop.txt')
    if file_path:
        cp_props = PARSER.parsCpProp(file_path[0])[0]
        cp_props = cp_props.children[0].children
        vis_cps = molsys.children[0].children.copy()
        vis_cps.sort(key=lambda x: x.coord[0])
        cp_props.sort(key=lambda x: x.coord[0])
        if len(vis_cps) != len(cp_props):
            return
        else:
            for i in range(len(cp_props)):
                vis_cp = vis_cps[i]
                cp_prop = cp_props[i]
                for key in cp_prop.__dict__:
                    if key[0] != '_' and key not in vis_cp.__dict__:
                        vis_cp.__dict__[key] = cp_prop.__dict__[key]
                    if vis_cp.point() is not None and key[0] != '_' and key not in vis_cp.point().__dict__:
                        vis_cp.point().addProperty(key, cp_prop.__dict__[key], False)
        a = 0
    else:
        return