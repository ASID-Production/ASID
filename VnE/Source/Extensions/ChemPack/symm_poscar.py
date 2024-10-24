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

from .ui.select_mol_dialog import SelectMolDialog


def execute(pars):
    from . import MOLECULE_SYSTEMS
    from PySide6 import QtWidgets

    def process(molsys):

        list_obj, molsys = molsys

        from pymatgen.core import Structure

        try:
            struct = Structure.from_file(molsys.file_name)
        except ValueError:
            return None, None

        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        sga = SpacegroupAnalyzer(struct, symprec=0.02)
        struct_symm_conv = sga.get_refined_structure()
        struct_symm_prim = sga.find_primitive()

        from pymatgen.io.vasp import Poscar

        filename_conv, _ = QtWidgets.QFileDialog.getSaveFileName(None, 'POSCAR-conv')
        if filename_conv:
            poscar = Poscar(struct_symm_conv)
            poscar.write_file(filename=filename_conv, significant_figures=16)
            pars(filename_conv, True)
        filename_prim, _ = QtWidgets.QFileDialog.getSaveFileName(None, 'POSCAR-prim')
        if filename_prim:
            poscar = Poscar(struct_symm_prim)
            poscar.write_file(filename=filename_prim, significant_figures=16)
            pars(filename_prim, True)

        return filename_conv, filename_prim

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()
    pass