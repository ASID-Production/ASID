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


from . import MoleculeClass
from ..ChemPack import PALETTE


class SaveFile:

    def getFormats(self):
        return SaveFile.__formats

    def save(self, molecule, filename, format):
        meth = SaveFile.__formats.get(format, None)
        if meth is not None:
            meth(self, molecule, filename)
        pass

    def save_xyz(self, molecule, filename):
        mol = molecule
        while not isinstance(mol.children[0], MoleculeClass.Atom):
            mol = molecule.children[0]

        out = open(filename, 'w')
        for atom in mol.children:
            line = f'{atom.atom_type}{atom.coord[0]:> 19.9f}{atom.coord[1]:> 17.9f}{atom.coord[2]:> 17.9f}\n'
            out.write(line)
        out.close()
        return

    def save_pdb(self, molecule, filename):
        mol = molecule
        while not isinstance(mol.children[0], MoleculeClass.Atom):
            mol = molecule.children[0]

        out = open(filename, 'w')
        i = 1
        prev_atom = mol.children[0]
        for atom in mol.children:
            if prev_atom.pdb_chain != atom.pdb_chain or prev_atom.pdb_flag != atom.pdb_flag:
                ter_line = f'TER   {i:>5}{prev_atom.pdb_res_name:>9} {prev_atom.pdb_chain}{prev_atom.pdb_res_seq:>4}{prev_atom.pdb_iCode}\n'
                out.write(ter_line)
                i += 1
            line = f'{atom.pdb_flag:<6}{i:>5} {atom.pdb_name:<4}{atom.pdb_alt_loc}{atom.pdb_res_name} {atom.pdb_chain}{atom.pdb_res_seq:>4}{atom.pdb_iCode}{atom.coord[0]:>11.3f}{atom.coord[1]:>8.3f}{atom.coord[2]:>8.3f}{atom.pdb_occupancy:>6.2f}{atom.pdb_tempFactor:>6.2f}{PALETTE.getName(atom.atom_type).upper():>12}{atom.pdb_charge}\n'
            out.write(line)
            prev_atom = atom
            i += 1
        out.write('END\n')
        out.close()
        return

    __formats = {'pdb': save_pdb,
                 'xyz': save_xyz}


SAVE_FILE = SaveFile()
