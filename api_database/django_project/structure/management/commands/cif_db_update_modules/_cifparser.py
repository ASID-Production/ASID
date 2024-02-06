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
# *****************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# *****************************************************************************************

from structure.models import StructureCode, CoordinatesBlock, Other, Cell
from ._element_numbers import element_numbers
import re


def add_coords(cif_block: dict, structure: classmethod):
    # if there are coordinates in the file, then add them to the database
    if {'_atom_site_label', '_atom_site_type_symbol',
            '_atom_site_fract_x', '_atom_site_fract_y',
            '_atom_site_fract_z'}.issubset(cif_block[1].keys()):
        coords, atoms = get_coords(cif_block)
        cb_obj, created = CoordinatesBlock.objects.get_or_create(refcode=structure)
        cb_obj.coordinates = coords
        cb_obj.save()
        return atoms
    return 0


def add_other_info(atoms, structure):
    other_obj, created = Other.objects.get_or_create(refcode=structure)
    max_atom_num = 0
    if atoms:
        for atom in atoms:
            atom_type = re.findall(r'[A-Za-z]{1,3}', atom)
            if atom_type and atom_type[0] in element_numbers.keys():
                num = element_numbers[atom_type[0]]
                if num > max_atom_num:
                    max_atom_num = num
        other_obj.number_atoms_with_sites = len(atoms)
        other_obj.maximum_atomic_number = max_atom_num
    else:
        other_obj.has_3d_structure = False
    other_obj.save()
    return 0


def add_cell_parms_with_error(cif_block, structure: classmethod):
    if {'_cell_length_a', '_cell_length_b', '_cell_length_c',
            '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'}.issubset(cif_block[1].keys()):
        a, b, c, al, be, ga = (cif_block[1]['_cell_length_a'], cif_block[1]['_cell_length_b'],
                               cif_block[1]['_cell_length_c'], cif_block[1]['_cell_angle_alpha'],
                               cif_block[1]['_cell_angle_beta'], cif_block[1]['_cell_angle_gamma'])
        if Cell.objects.filter(refcode=structure).exists():
            cell = Cell.objects.get(refcode=structure)
            cell.a_err = a
            cell.b_err = b
            cell.c_err = c
            cell.al_err = al
            cell.be_err = be
            cell.ga_err = ga
            cell.save()
    return 0


def get_coords(cif_block) -> str:
    atomic_sites = ''
    atoms = []
    coords = cif_block[1].GetLoop('_atom_site_label')
    order: list = coords.GetItemOrder()
    idxs = {
        'lable': order.index('_atom_site_label'),
        'atom_type_idx': order.index('_atom_site_type_symbol'),
        'x_idx': order.index('_atom_site_fract_x'),
        'y_idx': order.index('_atom_site_fract_y'),
        'z_idx': order.index('_atom_site_fract_z')
    }
    if '_atom_site_occupancy' in order:
        idxs['occup'] = order.index('_atom_site_occupancy')
    if '_atom_site_b_iso_or_equiv' in order:
        if 'occup' not in idxs.keys():
            idxs['occup'] = '1'
        idxs['b_iso'] = order.index('_atom_site_b_iso_or_equiv')
    for site in coords:
        atoms.append(site[idxs['lable']])
        temp = list()
        for value in idxs.values():
            temp.append(site[value])
        for j, element in enumerate(temp, start=0):
            # remove question marks
            if j == 0:
                temp[j] = element.replace('?', '')
            # remove the parentheses
            if '(' in element:
                idx = element.index('(')
                # rewrite the atomic element without parentheses
                temp[j] = element[:idx]
        atomic_sites += ' '.join(temp)
        atomic_sites += '\n'
    return atomic_sites, atoms
