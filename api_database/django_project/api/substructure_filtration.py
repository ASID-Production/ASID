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

from structure.management.commands.cif_db_update_modules._add_substructure_filtration import TEMPLATES, SET_ELEMENTS
from django.conf import settings
import cpplib

NUM_ELEM_DICT = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Me', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt'
}


def set_substructure(template_graph, analyse_mol):
    '''Return is the template_graph in analyse_mol.'''
    output = cpplib.SearchMain(template_graph, [analyse_mol, ], 1, False)
    return len(output)


def set_elements(analyse_mol):
    elem_found = set()
    mol_data = analyse_mol.split()
    num_nodes = int(mol_data[1])
    for idx, el in enumerate(mol_data[3:], start=1):
        if idx < num_nodes * 2 and idx % 2 == 1:
            if int(el) > 0:
                element = NUM_ELEM_DICT[int(el)]
                elem_found.add(element)
            else:
                return 0
    return elem_found


def set_filter(analyse_mol):
    result = dict()  # {attr_name: [1, obj_name]}

    for attr_name, data in TEMPLATES.items():
        template_graph, obj_name = data
        output = set_substructure(template_graph, analyse_mol)
        if output:
            result[attr_name] = [True, obj_name]
        else:
            result[attr_name] = [False, obj_name]

    elem_found = set_elements(analyse_mol)
    # only if no multitypes
    if elem_found:
        for attr_name, elements in SET_ELEMENTS.items():
            if elem_found.intersection(set(elements)):
                result[attr_name] = [True, 'Substructure1']
            else:
                result[attr_name] = [False, 'Substructure1']
        if {'C', 'N', 'O'}.issubset(elem_found):
            result['only_CHNO'] = [True, 'Substructure1']
        else:
            result['only_CHNO'] = [False, 'Substructure1']
        if 'C' not in elem_found:
            result['no_C'] = [True, 'Substructure1']
        else:
            result['no_C'] = [False, 'Substructure1']
    return result
