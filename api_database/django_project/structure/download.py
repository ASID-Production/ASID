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
#                                                                     #
#             Atomistic Simulation Instruments and Database           #
#                                ASID                                 #
#                                                                     #
#######################################################################
#                                                                     #
# This CIF file has been generated from the Atomistic Simulation      # 
# Instruments and Database and include bibliographic, chemical,       #
# crystal, experimental, refinement and atomic coordinate data.       #
#                                                                     #
#######################################################################

'''


def create_cif_text(structure: classmethod) -> str:

    def check_value_exist(cif_parameter: str, value, in_commas: bool, percent: bool = False) -> str:
        # check commas need
        if not in_commas and value and type(value) is str and len(value.split()) > 1:
            try:
                float(value)
            except ValueError:
                in_commas = True
        # check exists
        if percent and value:
            value = float(value) / 100
        if value and in_commas:
            return f"{cif_parameter} '{value}'\n"
        elif value and not in_commas:
            return f"{cif_parameter} {value}\n"
        elif cif_parameter == '_journal_name_full' and structure.publication.publication.journal.name:
            return f"'_journal_name_full' {structure.publication.publication.journal.name}\n"
        else:
            return f"{cif_parameter} ?\n"

    def check_cell_params_exist(cif_parameter: str, value, in_commas: bool, cell):
        unit_cell_params = {'_cell_length_a': cell.a, '_cell_length_b': cell.b,
                            '_cell_length_c': cell.c, '_cell_angle_alpha': cell.al,
                            '_cell_angle_beta': cell.be, '_cell_angle_gamma': cell.ga}
        # if there is no data of paramaters' diviations, write  parameters as is
        if value:
            return check_value_exist(cif_parameter, value, in_commas)
        return check_value_exist(cif_parameter, unit_cell_params[cif_parameter], in_commas)

    text = CIF_HEAD
    text += f"data_{structure.refcode}\n"
    text += check_value_exist('_database_code_CSD', structure.refcode, False)
    text += check_value_exist('_database_code_depnum_ccdc_archive', structure.CCDC_number, True)
    text += check_value_exist('_chemical_formula_sum', structure.formula.formula_sum, True)
    text += check_value_exist('_chemical_formula_moiety', structure.formula.formula_moiety, True)
    text += check_value_exist('_chemical_melting_point', structure.crystal_and_structure_info.melting_point, False)
    if hasattr(structure, "publication"):
        # text += check_value_exist('_journal_coden_Cambridge', structure.publication.publication.journal.international_coden, False)
        text += check_value_exist('_journal_volume', structure.publication.publication.volume, False)
        text += check_value_exist('_journal_year', structure.publication.publication.year, False)
        text += check_value_exist('_journal_page_first', structure.publication.publication.page, False)
        if structure.publication.publication.journal:
            text += check_value_exist('_journal_name_full', structure.publication.publication.journal.fullname, True)
        else:
            text += f"_journal_name_full ?\n"
        text += check_value_exist('_journal_DOI', structure.publication.publication.doi, True)

    if structure.authors.count() > 0:
        text += 'loop_\n_publ_author_name\n'
        for author in structure.authors.all():
            if author.initials:
                text += f'"{author.initials + author.family_name}"\n'
            else:
                text += f'"{author.family_name}"\n'

    text += check_value_exist('_chemical_name_systematic', structure.name.systematic_name, True)
    text += check_value_exist('_chemical_name_common', structure.name.trivial_name, True)

    volume = structure.cell.a * structure.cell.b * structure.cell.c * sqrt(
        1 + 2 * cos(radians(structure.cell.al)) * cos(radians(structure.cell.be)) * cos(radians(structure.cell.ga)) -
        cos(radians(structure.cell.al)) ** 2 - cos(radians(structure.cell.be)) ** 2 - cos(radians(structure.cell.ga))
    )

    text += check_value_exist('_cell_volume', round(volume, 3), False)
    text += check_value_exist('_exptl_crystal_density_diffrn', structure.experimental_info.calculated_density_value, False)
    text += check_value_exist('_exptl_crystal_colour', structure.crystal_and_structure_info.color, True)
    text += check_value_exist('_exptl_crystal_description', structure.crystal_and_structure_info.crystal_shape, True)

    if (structure.crystal_and_structure_info.bioactivity or
            structure.crystal_and_structure_info.phase_transitions or
            structure.crystal_and_structure_info.polymorph or
            structure.crystal_and_structure_info.sensitivity or
            structure.crystal_and_structure_info.pressure or
            structure.crystal_and_structure_info.disorder or
            structure.crystal_and_structure_info.recrystallisation_solvent):
        text += '_exptl_special_details\n;\n'
        for param, value in {
            'bioactivity': structure.crystal_and_structure_info.bioactivity,
            'phase transitions': structure.crystal_and_structure_info.phase_transitions,
            'polymorph': structure.crystal_and_structure_info.polymorph,
            'sensitivity': structure.crystal_and_structure_info.sensitivity,
            'pressure': structure.crystal_and_structure_info.pressure,
            'disorder': structure.crystal_and_structure_info.disorder,
            'recrystallisation solvent': structure.crystal_and_structure_info.recrystallisation_solvent,
        }.items():
            if (
                    (param == 'bioactivity' and structure.crystal_and_structure_info.bioactivity) or
                    (param == 'phase transitions' and structure.crystal_and_structure_info.phase_transitions) or
                    (param == 'polymorph' and structure.crystal_and_structure_info.polymorph) or
                    (param == 'sensitivity' and structure.crystal_and_structure_info.sensitivity) or
                    (param == 'pressure' and structure.crystal_and_structure_info.pressure) or
                    (param == 'disorder' and structure.crystal_and_structure_info.disorder) or
                    (param == 'recrystallisation solvent' and structure.crystal_and_structure_info.recrystallisation_solvent)
            ):
                text += f'{param}: {value}\n'
        text += ';\n'
    else:
        text += '_exptl_special_details ?\n'

    text += check_value_exist('_diffrn_ambient_temperature', structure.experimental_info.structure_determination_temperature, False)
    text += check_value_exist('_refine_ls_R_factor_gt', structure.refinement_info.r_factor, False, True)
    text += check_value_exist('_refine_ls_wR_factor_gt', structure.refinement_info.wR_factor, False, True)
    text += check_value_exist('_refine_ls_goodness_of_fit_ref', structure.refinement_info.gof, False, True)
    text += check_value_exist('_symmetry_cell_setting', structure.cell.spacegroup.get_system_display(), False)
    text += check_value_exist('_symmetry_space_group_name_H-M', structure.cell.spacegroup.name, True)
    text += check_value_exist('_symmetry_space_group_name_Hall', structure.cell.spacegroup.hall_name, True)
    text += check_value_exist('_symmetry_Int_Tables_number', structure.cell.spacegroup.number, False)

    text += 'loop_\n_symmetry_equiv_pos_site_id\n_symmetry_equiv_pos_as_xyz\n'
    for i, symop in enumerate(structure.cell.spacegroup.symops.split(';')):
        text += f'{i} {symop}\n'

    text += check_cell_params_exist('_cell_length_a', structure.cell.a_err, False, structure.cell)
    text += check_cell_params_exist('_cell_length_b', structure.cell.b_err, False, structure.cell)
    text += check_cell_params_exist('_cell_length_c', structure.cell.c_err, False, structure.cell)
    text += check_cell_params_exist('_cell_angle_alpha', structure.cell.al_err, False, structure.cell)
    text += check_cell_params_exist('_cell_angle_beta', structure.cell.be_err, False, structure.cell)
    text += check_cell_params_exist('_cell_angle_gamma', structure.cell.ga_err, False, structure.cell)
    text += check_value_exist('_cell_formula_units_Z', structure.cell.zvalue, False)
    text += check_value_exist('_cell_formula_units_Z_prime', structure.cell.zprime, False)

    if hasattr(structure, "coordinates"):
        if structure.characteristics.has_3d_structure and structure.coordinates.coordinates:
            # TODO: after dump loading leave only 5 values in coordinates strings!!!
            parms_num = len(structure.coordinates.coordinates.split('\n')[0].split())
            text += 'loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n'
            if parms_num >= 6:
                text += '_atom_site_occupancy\n'
            if parms_num >= 7:
                text += '_atom_site_B_iso_or_equiv\n'
            text += structure.coordinates.coordinates

    text += '\n#END\n'
    return text
