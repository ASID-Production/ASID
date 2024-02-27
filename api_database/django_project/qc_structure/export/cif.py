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

from math import cos, sqrt, radians

CIF_HEAD = '''
#######################################################################
#
#             Atomistic Simulation Instruments and Database
#                                ASID
#
#######################################################################
#
# This CIF file has been generated from the Atomistic Simulation 
# Instruments and Database and include bibliographic, chemical, 
# crystal, experimental, refinement and atomic coordinate data.
#
#######################################################################

'''


def qc_get_cif_content(qc_structure: classmethod) -> str:

    def check_value_exist(cif_parameter: str, value, in_commas: bool) -> str:
        if value:
            if in_commas:
                return f"{cif_parameter} '{value}'\n"
            return f"{cif_parameter} {value}\n"
        return f"{cif_parameter} ?\n"

    text = CIF_HEAD
    text += f"data_{qc_structure.refcode}\n"
    text += check_value_exist('_database_code', qc_structure.refcode, False)
    text += check_value_exist('_chemical_formula_sum', qc_structure.qc_formula.formula_sum, True)
    text += check_value_exist('_chemical_formula_moiety', qc_structure.qc_formula.formula_moiety, True)

    text += check_value_exist('_chemical_name_systematic', qc_structure.qc_name.systematic_name, True)
    text += check_value_exist('_chemical_name_common', qc_structure.qc_name.trivial_name, True)

    volume = qc_structure.qc_cell.a * qc_structure.qc_cell.b * qc_structure.qc_cell.c * sqrt(
        1 + 2 * cos(radians(qc_structure.qc_cell.al)) * cos(radians(qc_structure.qc_cell.be)) * cos(radians(qc_structure.qc_cell.ga)) -
        cos(radians(qc_structure.qc_cell.al)) ** 2 - cos(radians(qc_structure.qc_cell.be)) ** 2 - cos(radians(qc_structure.qc_cell.ga))
    )

    text += check_value_exist('_cell_volume', round(volume, 3), False)
    text += check_value_exist('_exptl_crystal_density_diffrn', qc_structure.qc_properties.calculated_density, False)
    text += '_exptl_special_details ?\n'

    text += check_value_exist('_symmetry_cell_setting', qc_structure.qc_cell.spacegroup.get_system_display(), False)
    text += check_value_exist('_symmetry_space_group_name_H-M', qc_structure.qc_cell.spacegroup.name, True)
    text += check_value_exist('_symmetry_space_group_name_Hall', qc_structure.qc_cell.spacegroup.hall_name, True)
    text += check_value_exist('_symmetry_Int_Tables_number', qc_structure.qc_cell.spacegroup.number, False)

    text += 'loop_\n_symmetry_equiv_pos_site_id\n_symmetry_equiv_pos_as_xyz\n'
    for i, symop in enumerate(qc_structure.qc_cell.spacegroup.symops.split(';')):
        text += f'{i} {symop}\n'

    text += check_value_exist('_cell_length_a', qc_structure.qc_cell.a, False)
    text += check_value_exist('_cell_length_b', qc_structure.qc_cell.b, False)
    text += check_value_exist('_cell_length_c', qc_structure.qc_cell.c, False)
    text += check_value_exist('_cell_angle_alpha', qc_structure.qc_cell.al, False)
    text += check_value_exist('_cell_angle_beta', qc_structure.qc_cell.be, False)
    text += check_value_exist('_cell_angle_gamma', qc_structure.qc_cell.ga, False)
    text += check_value_exist('_cell_formula_units_Z', qc_structure.qc_cell.zvalue, False)

    if qc_structure.qc_coordinates.coordinates:
        text += 'loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n'
        text += qc_structure.qc_coordinates.coordinates

    text += '#END\n'
    return text
