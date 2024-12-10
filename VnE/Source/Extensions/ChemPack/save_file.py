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

import debug


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
        out.write(f'{len(mol.children)}\n')
        out.write('\n')
        for atom in mol.children:
            line = f'{PALETTE.getName(atom.atom_type)}{atom.coord[0]:> 19.9f}{atom.coord[1]:> 17.9f}{atom.coord[2]:> 17.9f}\n'
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

    def save_cif(self, molecule, filename):
        from pymatgen.core import Structure
        from pymatgen.io import cif
        import numpy as np
        from .MoleculeClass import DefaultData

        def cell_to_vec(cell):

            a, b, c, al, be, ga = cell
            cos = lambda x: np.cos(x/180*np.pi)
            sin = lambda x: np.sin(x/180*np.pi)
            tg = lambda x: np.tan(x/180*np.pi)
            cot = lambda x: 1/np.tan(x/180*np.pi)
            csc = lambda x: 1/np.sin(x/180*np.pi)

            mat = np.array([[a*sin(be)*np.sqrt(1-(cot(al)*cot(be)-csc(al)*csc(be)*cos(ga))**2), 0, 0],
                            [a*csc(al)*cos(ga)-a*cot(al)*cos(be), b*sin(al), 0],
                            [a*cos(be), b*cos(al), c]])
            return mat

        mol = molecule
        while not isinstance(mol.children[0], MoleculeClass.Atom):
            mol = molecule.children[0]

        atm = mol.children[0]
        cell_st = [DefaultData.__getattr__(DefaultData, 'cif_cell_a'),
                   DefaultData.__getattr__(DefaultData, 'cif_cell_b'),
                   DefaultData.__getattr__(DefaultData, 'cif_cell_c'),
                   DefaultData.__getattr__(DefaultData, 'cif_cell_al'),
                   DefaultData.__getattr__(DefaultData, 'cif_cell_be'),
                   DefaultData.__getattr__(DefaultData, 'cif_cell_ga')]
        cell = [atm.cif_cell_a, atm.cif_cell_b, atm.cif_cell_c,
                atm.cif_cell_al, atm.cif_cell_be, atm.cif_cell_ga]
        coords = np.array([x.coord.copy() for x in mol])
        sym_codes =atm.cif_sym_codes
        exclude = ['x+1,y,z',
                   'x-1,y,z',
                   'x,y+1,z',
                   'x,y-1,z',
                   'x,y,z+1',
                   'x,y,z-1']
        sym_codes = [x for x in atm.cif_sym_codes if x[1] not in exclude]
        space_group = atm.cif_space_group
        if cell == cell_st:
            minimum = np.array([min([x[i] for x in coords]) for i in range(3)])
            coords -= minimum
            maximum = np.array([max([x[i] for x in coords]) for i in range(3)])
            cell[0] = maximum[0] + 6
            cell[1] = maximum[1] + 6
            cell[2] = maximum[2] + 6
            coords += 3

        lattice = cell_to_vec(cell)
        species = [x.atom_type for x in mol]
        coords = coords@np.linalg.inv(lattice)
        labels = [x.name for x in mol]
        struct = Structure(lattice, species, coords, to_unit_cell=False, coords_are_cartesian=False, labels=labels)
        writer = cif.CifWriter(struct)
        cif_file = writer.cif_file.data
        f_mol = list(cif_file.keys())[0]
        data = cif_file[f_mol].data
        data['_symmetry_space_group_name_H-M'] = space_group
        data['_symmetry_equiv_pos_site_id'] = [str(x[0]) for x in sym_codes]
        data['_symmetry_equiv_pos_as_xyz'] = [x[1] for x in sym_codes]
        writer.write_file(filename)

    __formats = {'pdb': save_pdb,
                 'xyz': save_xyz,
                 'cif': save_cif,
                 }


SAVE_FILE = SaveFile()
