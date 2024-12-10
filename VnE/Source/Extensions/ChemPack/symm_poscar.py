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
        symm_struct = sga.get_symmetrized_structure()

        atoms = []
        coords = []
        for atm in symm_struct.equivalent_sites:
            atm = atm[0]
            atoms.append(atm.specie)
            coords.append(atm.frac_coords)
        struct = Structure(symm_struct.lattice, atoms, coords)

        from pymatgen.io.vasp import Poscar

        filename_conv, _ = QtWidgets.QFileDialog.getSaveFileName(None, 'POSCAR symm')
        sym_codes = []
        for i, op in enumerate(symm_struct.spacegroup):
            sym_codes.append([i, op.as_xyz_str()])
        added_sym_codes = ['x+1,y,z',
                           'x-1,y,z',
                           'x,y+1,z',
                           'x,y-1,z',
                           'x,y,z+1',
                           'x,y,z-1']
        a = i + 1
        for i, op in enumerate(added_sym_codes):
            sym_codes.append([a+i, op])
        cif_data = {'cif_sym_codes': sym_codes,
                    'cif_space_group': symm_struct.spacegroup.int_symbol}

        if filename_conv:
            poscar = Poscar(struct)
            poscar.write_file(filename=filename_conv, significant_figures=16)
            mol_sys, lists_tuple = pars(filename_conv, True)
            mol = mol_sys.children[0]
            for atom in mol:
                for key, val in cif_data.items():
                    atom.__setattr__(key, val)


        return filename_conv

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()
    pass