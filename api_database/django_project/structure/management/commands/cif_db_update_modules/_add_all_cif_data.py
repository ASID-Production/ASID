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

from structure.models import (StructureCode, Author, Spacegroup, SYSTEMS,
                              Cell, CENTRINGS, ReducedCell, Journal,
                              Publication, RefcodePublicationConnection, Formula,
                              get_elements_list, ElementsManager,
                              CompoundName, ExperimentalInfo, RefinementInfo,
                              CrystalAndStructureInfo)
from django_project.loggers import all_cif_data_logger as logger_1
import re
from api.filters import get_reduced_cell
from math import cos, sqrt, radians
from ._cifparser import add_cell_parms_with_error
import json
import os

compound_names = {
    'systematic_name': ['_chemical_name_systematic'],
    'trivial_name': ['_chemical_name_common', '_chemical_name_mineral']
}
experimental_info = {
    'structure_determination_temperature': [
        '_diffrn_ambient_temperature',
        '_cell_measurement_temperature',
    ],
    'measurement_reflns': [
        '_diffrn_reflns_number',
        '_refine_ls_number_reflns',
        '_reflns_number_total',
    ],
    'measurement_theta_max': [
        '_cell_measurement_theta_max',
        '_diffrn_reflns_theta_max'
    ],
    'measurement_theta_min': [
        '_cell_measurement_theta_min',
        '_diffrn_reflns_theta_min'
    ],
    'calculated_density_value': ['_exptl_crystal_density_diffrn'],
    'absorpt_correction_type': ['_exptl_absorpt_correction_type'],
    'measurement_device': ['_diffrn_measurement_device'],
    'measurement_method': ['_diffrn_measurement_method'],
    'monochromator': ['_diffrn_radiation_monochromator'],
    'radiation_type': ['_diffrn_radiation_type'],
    'wavelength': [
        '_diffrn_radiation_wavelength',
        '_cell_measurement_wavelength',
        '_diffrn_refln_wavelength',
        '_refln_wavelength'
    ],
    'radiation_source': ['_diffrn_radiation_source'],
}
refinement_info = {
    'r_factor': ['_refine_ls_r_factor_gt'],
    'wR_factor': ['_refine_ls_wr_factor_ref'],
    'gof': ['_refine_ls_goodness_of_fit_ref'],
    'diff_density_max': ['_refine_diff_density_max'],
    'diff_density_min': ['_refine_diff_density_min'],
    'extinction_coef': ['_refine_ls_extinction_coef'],
    'flack': ['_refine_ls_abs_structure_flack'],
}
crystal_and_structure_info = {
    'color': [
        '_exptl_crystal_colour',
        '_exptl_crystal_colour_primary',
    ],
    'bioactivity': ['_exptl_special_details'],
    'crystal_shape': ['_exptl_crystal_description'],
    'phase_transitions': ['_exptl_special_details'],
    'polymorph': ['_exptl_special_details'],
    'melting_point': [
        '_chemical_melting_point',
        '_exptl_special_details'
    ],
    'sensitivity': ['_exptl_special_details'],
    'pressure': [
        '_cell_measurement_pressure',
        '_diffrn_ambient_pressure',
        '_exptl_special_details'
    ],
    'disorder': ['_exptl_special_details'],
    'recrystallisation_solvent': [
        '_exptl_crystal_recrystallization_method',
        '_exptl_crystal_preparation'],
    'size_max': ['_exptl_crystal_size_max'],
    'size_min': ['_exptl_crystal_size_min'],
    'size_mid': ['_exptl_crystal_size_mid'],
}
formula = {
    'formula_moiety': [
        '_chemical_formula_moiety',
        '_chemical_formula_structural',
        '_chemical_formula_iupac'
    ],
    'formula_sum': ['_chemical_formula_sum']
}
journal = {
    'international_coden': ['_journal_coden_cambridge'],
    'name': ['_citation_journal_abbrev'],
    'fullname': ['_journal_name_full'],
}


def add_refcode(cif_block, refcode):
    data_in_cif = cif_block.keys()
    struct_obj, created = StructureCode.objects.get_or_create(refcode=refcode)
    logger_1.info(f'Refcode: {refcode}. Created new: {created}')
    if '_database_code_depnum_ccdc_archive' in data_in_cif:
        ccdc = re.findall(r'\d+', cif_block['_database_code_depnum_ccdc_archive'])[0]
        struct_obj.CCDC_number = str(ccdc)
    if '_database_code_icsd' in data_in_cif:
        struct_obj.ICSD = True
    elif '_cod_database_code' in data_in_cif:
        struct_obj.COD = True
    struct_obj.save()
    return struct_obj


def add_author(cif_block, struct_obj, authors=None):
    if not authors:
        try:
            authors = []
            authors_temp = cif_block.GetLoop('_publ_author_name')
            for author in authors_temp:
                authors.append(author[0])
        except:
            if '_publ_author_name' in cif_block.keys():
                authors = cif_block['_publ_author_name'].split(';')
            else:
                return 0
    # add to database
    for author in authors:
        author = str(author).replace('\n', '')
        author = author.replace('\r', '')
        if author[0] in ['\'', '"']:
            author = author[1:]
        if author[-1] in ['\'', '"']:
            author = author[:-1]
        author_split = re.findall('[^, .]+', author)
        if len(author) == 1:
            author_obj, created = Author.objects.get_or_create(family_name=author_split[0])
        elif len(author) == 2:
            if ',' in author:
                family, initials = author_split
            else:
                family = author_split[-1]
                initials = author.replace(family, '')
            author_obj, created = Author.objects.get_or_create(family_name=family, initials=initials)
        else:
            if ',' in author:
                family = author.split(',')[0]
                initials = ' '.join(re.findall('[^, ]+', author)[1:])
            else:
                family = author_split[-1]
                initials = author.replace(family, '')
            author_obj, created = Author.objects.get_or_create(family_name=family, initials=initials)
        logger_1.info(f'Author: {author_obj}')
        struct_obj.authors.add(author_obj)
        return author_obj
    return 0


def get_or_create_space_group(cif_block):

    def get_sg_number_by_h_m_name(h_m_name: str):
        file = open(os.path.join(os.path.dirname(__file__), '../symops.json'))
        symops_list = json.load(file)
        for i, symops in enumerate(symops_list):
            if i > 0:
                if symops['hermann_mauguin'] == h_m_name:
                    return symops['number']

    def get_system_from_numb(sg_number: int):
        file = open(os.path.join(os.path.dirname(__file__), '../symops.json'))
        symops_list = json.load(file)
        for i, symops in enumerate(symops_list):
            if i > 0:
                if symops['number'] == sg_number:
                    return symops['crystal_class']

    data_in_cif = cif_block.keys()
    # sg number
    if '_space_group_it_number' in data_in_cif:
        sg_number = int(cif_block['_space_group_it_number'])
    elif '_symmetry_int_tables_number' in data_in_cif:
        sg_number = int(cif_block['_symmetry_int_tables_number'])
    elif '_symmetry_space_group_name_h-m' in data_in_cif:
        sg_number = get_sg_number_by_h_m_name(cif_block['_symmetry_space_group_name_h-m'])
    else:
        raise Exception('No space group number found in cif!')
    logger_1.info(f'Space group number: {sg_number}')
    # system
    if '_space_group_crystal_system' in data_in_cif:
        system = cif_block['_space_group_crystal_system']
    elif '_symmetry_cell_setting' in data_in_cif:
        system = cif_block['_symmetry_cell_setting']
    else:
        system = get_system_from_numb(sg_number)
    logger_1.info(f'Crystal system: {system}')
    # get system number in db choices
    if not system.isalpha():
        system = get_system_from_numb(sg_number)
    if system:
        system_id = 0
        for syst in SYSTEMS:
            if system in syst:
                system_id = syst[0]
        if system_id:
            logger_1.info(f'Space group system number: {system_id}')
        else:
            raise Exception('No space group system found in cif!')
    else:
        raise Exception('No space group system found in cif!')
    # check if space group has already existed in db
    if '_symmetry_space_group_name_h-m' in data_in_cif or '_space_group_name_h-m_alt' in data_in_cif:
        if '_symmetry_space_group_name_h-m' in data_in_cif:
            space_group_name = cif_block['_symmetry_space_group_name_h-m'].split()
        else:
            space_group_name = cif_block['_space_group_name_h-m_alt'].split()
        if system == 'monoclinic':
            if space_group_name[1] == '1' and space_group_name[3] == '1':
                space_group_name = [space_group_name[0], space_group_name[2]]
        space_group_name = ''.join(space_group_name)
        logger_1.info(f'Space group name: {space_group_name}')
    else:
        raise Exception('No space group found in cif!')
    hall = ''
    if '_symmetry_space_group_name_hall' in data_in_cif or '_space_group_name_hall' in data_in_cif:
        if '_symmetry_space_group_name_hall' in data_in_cif:
            hall = cif_block['_symmetry_space_group_name_hall']
        else:
            hall = cif_block['_space_group_name_hall']
        logger_1.info(f'Hall name: {hall}')
        if Spacegroup.objects.filter(name=space_group_name, hall_name=hall).exists():
            space_group = Spacegroup.objects.get(name=space_group_name, hall_name=hall)
            return space_group
    else:
        space_group = Spacegroup.objects.filter(name=space_group_name)
        if space_group.count() == 1:
            return space_group[0]
        elif space_group.count() > 1:
            raise Exception('No space group found in cif!')
        elif space_group.count() == 0:
            # if such a group is not in the database, then add it
            logger_1.info(f'Creating new space group object...')
    # create a space group object
    space_group, create = Spacegroup.objects.get_or_create(
        name=space_group_name,
        number=int(sg_number),
        system=system_id
    )
    # hall name
    if hall:
        space_group.hall_name = hall
    # symops
    try:
        symops = cif_block.GetLoop('_symmetry_equiv_pos_as_xyz')
        key = '_symmetry_equiv_pos_as_xyz'
    except:
        try:
            symops = cif_block.GetLoop('_space_group_symop_operation_xyz')
            key = '_space_group_symop_operation_xyz'
        except:
            space_group.delete()
            raise Exception('No symops found in cif!')
    keys: list = symops.GetItemOrder()
    symops_list = []
    for item in symops:
        symop = item[keys.index(key)].replace(' ', '')
        symops_list.append(symop)
    operations = ";".join(symops_list)
    space_group.symops = operations
    space_group.save()
    return space_group


def add_cell(cif_block, space_group, struct_obj):
    data_in_cif = cif_block.keys()
    if {
        '_cell_length_a', '_cell_length_b', '_cell_length_c',
        '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'
    }.issubset(data_in_cif):
        a, b, c, al, be, ga = (cif_block['_cell_length_a'].split('(')[0],
                               cif_block['_cell_length_b'].split('(')[0],
                               cif_block['_cell_length_c'].split('(')[0],
                               cif_block['_cell_angle_alpha'].split('(')[0],
                               cif_block['_cell_angle_beta'].split('(')[0],
                               cif_block['_cell_angle_gamma'].split('(')[0])
    else:
        raise Exception('No unit cell parameters were found!')
    if space_group.name[0].upper() in ['P', 'I', 'A', 'B', 'C', 'F', 'R']:
        centring = space_group.name[0]
    elif space_group.name[1].upper() in ['P', 'I', 'A', 'B', 'C', 'F', 'R']:
        centring = space_group.name[1]
    else:
        raise Exception('Invalid space group name: could not find centring in space group name!')
    centring_id = 0
    for item in CENTRINGS:
        if centring.replace('-', '') in item:
            centring_id = item[0]
    if not centring_id:
        raise Exception('Invalid unit cell centring!')
    if '_cell_formula_units_z' in data_in_cif:
        z_val = float(cif_block['_cell_formula_units_z'])
    else:
        logger_1.error('No Z value was found!')
        z_val = len(space_group.symops.split(';'))
    cell, created = Cell.objects.get_or_create(
        a=float(a), b=float(b), c=float(c),
        al=float(al), be=float(be), ga=float(ga),
        spacegroup=space_group, refcode=struct_obj,
        centring=centring_id, zvalue=z_val
    )
    if '_cell_formula_units_Z_prime' in data_in_cif:
        z_prime = float(cif_block['_cell_formula_units_Z_prime'])
        cell.zprime = z_prime
        cell.save()
    if created:
        add_cell_parms_with_error([1, cif_block], struct_obj)


def add_reduced_cell(struct_obj):
    centrings = dict((v, k) for v, k in CENTRINGS)
    centring = centrings[struct_obj.cell.centring]
    params = [struct_obj.cell.a, struct_obj.cell.b, struct_obj.cell.c,
              struct_obj.cell.al, struct_obj.cell.be, struct_obj.cell.ga]
    reduced_params = get_reduced_cell(params, centring)
    a, b, c, al, be, ga = reduced_params
    volume = (
            a * b * c * sqrt(1 + 2 * cos(radians(al)) * cos(radians(be)) *
            cos(radians(ga)) - cos(radians(al)) ** 2 -
            cos(radians(be)) ** 2 - cos(radians(ga)) ** 2)
    )
    rc, created = ReducedCell.objects.get_or_create(
        refcode=struct_obj,
        a=round(a, 3), b=round(b, 3), c=round(c, 3),
        al=round(al, 3), be=round(be, 3), ga=round(ga, 3),
        volume=round(volume, 3),
    )


def add_journal_and_publication(cif_block, struct_obj):

    def add_authors_and_get_flag():
        for author in authors:
            author_obj = add_author(cif_block[1], struct_obj, authors=[author, ])
            if created:
                pub_obj.authors.add(author_obj)
            if author_obj in pub_obj.authors.all():
                continue
            else:
                return False
        return True

    filtr_j = dict()
    for key, values in journal.items():
        for data_key in values:
            if data_key in cif_block[1].keys():
                if data_key == '_citation_journal_abbrev':
                    citation = cif_block[1].GetLoop('_citation_journal_abbrev')
                    value = citation['_citation_journal_abbrev'][0]
                else:
                    value = cif_block[1][data_key]
                if value:
                    if key == 'name':
                        filtr_j['name'] = value
                    elif key == 'fullname':
                        filtr_j['fullname'] = value
    j_obj = Journal.objects.filter(**filtr_j)
    filtr = dict()
    if j_obj.count() == 1:
        journal_obj = j_obj[0]
        filtr['journal'] = journal_obj
    elif j_obj.count() == 0:
        journal_obj = Journal.objects.create(**filtr_j)
        filtr['journal'] = journal_obj
    # publication
    year = 0
    if '_journal_year' in cif_block[1].keys():
        year = int(cif_block[1]['_journal_year'])
        filtr['year'] = year
    if '_journal_page_first' in cif_block[1].keys():
        page = cif_block[1]['_journal_page_first']
        filtr['page'] = page
    if '_journal_volume' in cif_block[1].keys():
        volume = cif_block[1]['_journal_volume']
        filtr['volume'] = volume
    if '_journal_doi' in cif_block[1].keys():
        doi = cif_block[1]['_journal_doi']
        filtr['doi'] = doi
    if '_citation_year' in cif_block[1].keys():
        citation = cif_block[1].GetLoop('_citation_year')
        citation_keys = citation.GetItemOrder()
        year = int(citation['_citation_year'][0])
        filtr['year'] = year
        if '_citation_doi' in citation_keys:
            doi = citation['_citation_doi'][0]
            filtr['doi'] = doi
        if '_citation_journal_volume' in citation_keys:
            volume = citation['_citation_journal_volume'][0]
            filtr['volume'] = volume
        if '_citation_page_first' in citation_keys:
            page = citation['_citation_page_first'][0]
            filtr['page'] = page
    # filter is_null
    filtr_is_null = dict()
    for key in ['journal', 'page', 'volume', 'doi']:
        if key not in filtr.keys():
            filtr_is_null[f'{key}__isnull'] = True
    if year:
        pub_obj, created = Publication.objects.filter(**filtr_is_null).get_or_create(**filtr)
    else:
        return 0
    # authors
    flag = False
    authors = []
    if '_publ_author_name' in cif_block[1].keys():
        try:
            authors_temp = cif_block[1].GetLoop('_publ_author_name')
            for author in authors_temp:
                authors.append(author[0])
        except:
            authors = cif_block[1]['_publ_author_name'].split(';')
        flag = add_authors_and_get_flag()
    if '_citation_author_name' in cif_block[1].keys():
        try:
            authors = cif_block[1].GetLoop('_citation_author_name')
        except:
            authors = cif_block[1]['_citation_author_name'].split(';')
        flag = add_authors_and_get_flag()
    # refcode connection
    if flag:
        ref_pub_obj, created = RefcodePublicationConnection.objects.get_or_create(refcode=struct_obj)
        ref_pub_obj.publication = pub_obj
        ref_pub_obj.save()


def add_element_composition(cif_block, struct_obj):
    '''Must be call after add_formula function!!!'''
    if Formula.objects.filter(refcode=struct_obj).exists():
        formula_sum = struct_obj.formula.formula_sum
        formula_moiety = struct_obj.formula.formula_moiety
        elements_from_formula = dict()
        if formula_sum:
            elements = formula_sum.split()
            for element in elements:
                atom_type = re.findall(r'[A-Za-z]{1,3}', element)[0]
                count = re.findall(r'\d+', element)
                if count:
                    count = float(count[0])
                else:
                    count = 1
                elements_from_formula[atom_type] = count
            el_manager, created = ElementsManager.objects.get_or_create(refcode=struct_obj)
            el_manager.save_elements(elements_from_formula)
            return 0
        elif formula_moiety:
            moiety_formula = dict()
            mols = formula_moiety.split(',')
            for mol in mols:
                not_splited = False
                if '(' in mol:
                    multipl = re.findall(r'^\d+', mol)
                    if multipl and mol.startswith(multipl[0] + '('):
                        multipl_coof = float(multipl[0])
                    elif mol.startswith('('):
                        multipl = re.findall(r'[)]\d+', mol)
                        if multipl and mol.endswith(multipl[0]):
                            multipl_coof = float(multipl[0].replace(')', ''))
                        elif mol.startswith('(') and mol.endswith(')'):
                            multipl_coof = 1
                        else:
                            not_splited = True
                    else:
                        not_splited = True
                    if not not_splited:
                        mol = mol.split('(')[1].split(')')[0].split()
                        for item in mol:
                            if re.search(r'[A-Za-z]', item):
                                atom_type = re.findall(r'[A-Za-z]{1,3}', item)[0]
                                count = re.findall(r'\d+', item)
                                if count:
                                    count = float(count[0])
                                else:
                                    count = 1
                                if atom_type in moiety_formula.keys():
                                    moiety_formula[atom_type] += count * multipl_coof
                                else:
                                    moiety_formula[atom_type] = count * multipl_coof
                elif not not_splited:
                    if re.search(r'[A-Za-z]', mol):
                        atom_type = re.findall(r'[A-Za-z]{1,3}', mol)[0]
                        count = re.findall(r'\d+', mol)
                        if count:
                            count = float(count[0])
                        else:
                            count = 1
                        if atom_type in moiety_formula.keys():
                            moiety_formula[atom_type] += count
                        else:
                            moiety_formula[atom_type] = count
                if not_splited:
                    mol = mol.split()
                    for item in mol:
                        if '(' in item:
                            multipl = re.findall(r'^\d+', item)
                            if multipl and item.startswith(multipl[0] + '('):
                                multipl_coof = float(multipl[0])
                            elif item.startswith('('):
                                multipl = re.findall(r'[)]\d+', item)
                                if multipl and item.endswith(multipl[0]):
                                    multipl_coof = float(multipl[0].replace(')', ''))
                                elif item.startswith('(') and item.endswith(')'):
                                    multipl_coof = 1
                                else:
                                    continue
                            else:
                                continue
                            item = item.split('(')[1].split(')')[0].split()
                            for elem in item:
                                if re.search(r'[A-Za-z]', elem):
                                    atom_type = re.findall(r'[A-Za-z]{1,3}', elem)[0]
                                    count = re.findall(r'\d+', elem)
                                    if count:
                                        count = float(count[0])
                                    else:
                                        count = 1
                                    if atom_type in moiety_formula.keys():
                                        moiety_formula[atom_type] += count * multipl_coof
                                    else:
                                        moiety_formula[atom_type] = count * multipl_coof
                        else:
                            if re.search(r'[A-Za-z]', item):
                                atom_type = re.findall(r'[A-Za-z]{1,3}', item)[0]
                                count = re.findall(r'\d+', item)
                                if count:
                                    count = float(count[0])
                                else:
                                    count = 1
                                if atom_type in moiety_formula.keys():
                                    moiety_formula[atom_type] += count
                                else:
                                    moiety_formula[atom_type] = count
            if moiety_formula:
                el_manager, created = ElementsManager.objects.get_or_create(refcode=struct_obj)
                el_manager.save_elements(moiety_formula)
                return 0
    logger_1.warring(
        f'Any chemcal composition was not found\n'
        f'or formula sum and formula moiety have invalid format\n'
        f'and can not be parsed!'
    )
    return 1


def add_universal(model, refcode_obj, cif_block, iucr_dict):
    obj_obj, created = model.objects.get_or_create(refcode=refcode_obj)
    for key, values in iucr_dict.items():
        for data_key in values:
            if data_key in cif_block[1].keys():
                value = cif_block[1][data_key]
                if value == '?':
                    break
                if '(' in value and ')' in value and re.search(r'\d', value) and not re.search(r'[a-zA-Z]', value):
                    value = float(value.split('(')[0])
                if data_key == '_exptl_special_details':
                    for line in value.split('\n'):
                        if key == 'bioactivity' and 'drug' in line:
                            value = line
                            break
                        elif key == 'phase_transitions' and 'phase' in line:
                            value = line
                            break
                        elif key == 'polymorph' and 'polymorph' in line:
                            value = line
                            break
                        elif key == 'sensitivity' and {'decomposition', 'sensitve', 'oxidant', 'stable', 'chromic', 'hygroscopic', 'thermolabile'}.intersection(set(line.split())):
                            value = line
                            break
                        elif key == 'pressure' and ('pressure' in line or re.search(r'[KkGgMmPa]{3}', line)):
                            value = line
                            break
                        elif key == 'disorder' and 'disorder' in line:
                            value = line
                            break
                        elif key == 'melting_point' and {'sublimes', 'above'}.intersection(set(line.split())):
                            value = line
                            break
                        else:
                            value = ''
                if key == 'systematic_name':
                    value = ' '.join(value.split())
                if key == 'wavelength' and type(value) == list:
                    value = float(value[0].split(',')[0])
                if key == 'extinction_coef':
                    value = float(value.split('(')[0])
                if key == 'formula_moiety' or key == 'formula_sum':
                    value = value.replace('\n', '')
                    value = value.replace('\r', '')
                if key in ['r_factor', 'wR_factor', 'gof']:
                    value = float(value) * 100
                if value and value not in ['?', 'none']:
                    setattr(obj_obj, key, value)
                    break
    # if problems with attr values (first of all for COD structures)
    try:
        obj_obj.save()
    except Exception as err:
        key = str(err).split()[1].replace("'", '')
        delattr(obj_obj, key)
        obj_obj.save()
    return obj_obj


def add_all_cif_data(cifs: dict):
    """
    Website https://xstar.sourceforge.net/astar/sf/output/cif_core.htm
    was used to parse the cif file codes
    """
    cif_blocks = cifs.copy()
    for refcode, cif_block in cifs.items():
        logger_1.info(f'{refcode} in progres...')
        struct_obj = add_refcode(cif_block[1], refcode)
        try:
            logger_1.info(f'Add authors')
            add_author(cif_block[1], struct_obj)
            logger_1.info(f'Check space group')
            space_group = get_or_create_space_group(cif_block[1])
            logger_1.info(f'Space group: {space_group}')
            logger_1.info(f'Add unit cell')
            add_cell(cif_block[1], space_group, struct_obj)
            logger_1.info(f'Add reduced cell')
            add_reduced_cell(struct_obj)
            logger_1.info(f'Add compound name')
            add_universal(CompoundName, struct_obj, cif_block, compound_names)
            logger_1.info(f'Add experimental info')
            add_universal(ExperimentalInfo, struct_obj, cif_block, experimental_info)
            logger_1.info(f'Add refinement info')
            add_universal(RefinementInfo, struct_obj, cif_block, refinement_info)
            logger_1.info(f'Add crystal and structure info')
            add_universal(CrystalAndStructureInfo, struct_obj, cif_block, crystal_and_structure_info)
            logger_1.info(f'Add formula')
            add_universal(Formula, struct_obj, cif_block, formula)
            logger_1.info(f'Add journal and publication')
            add_journal_and_publication(cif_block, struct_obj)
            logger_1.info(f'Add element composition')
            add_element_composition(cif_block, struct_obj)
            logger_1.info(f'Information addition from {refcode} is completed successfully!')
        except Exception as err:
            message = f'Caught an exception for structure {refcode}\n{err}'
            logger_1.error(message)
            # delete object if something went wrong
            struct_obj.delete()
            cif_blocks.pop(refcode)
            # if only 1 structure was upload raise error
            if len(cifs) == 1:
                raise Exception(err)
    return cif_blocks
