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

from structure.models import (
    Substructure1, Substructure2, CoordinatesBlock, StructureCode,
    ElementsSet1, ElementsSet2, ElementsSet3, ElementsSet4, ElementsSet5,
    ElementsSet6, ElementsSet7, ElementsSet8
)
from django.db.models import Q
from django.conf import settings
from django_project.loggers import substructure_logger
import cpplib


TEMPLATES = {
    'NO2': ['1 3 2 8 0 1 3 7 0 3 3 8 0 1 14 1 2 2 3', 'Substructure1'],
    'SO2': ['1 3 2 8 0 1 14 16 0 2 14 8 0 1 14 1 2 2 3', 'Substructure1'],
    # 'C_Met': ['', 'Substructure1'],
    # 'N_Met': ['', 'Substructure1'],
    'CS': ['1 2 1 16 0 1 14 6 0 1 14 1 2', 'Substructure1'],
    'ring6': ['1 6 6 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring5': ['1 5 5 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'Ph': ['1 6 6 6 1 3 3 6 1 3 3 6 1 3 3 6 0 3 3 6 1 3 3 6 1 3 3 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring6N1': ['1 6 6 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 7 0 2 14 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring6N2': ['1 6 6 6 0 2 14 6 0 2 14 7 0 2 14 6 0 2 14 6 0 2 14 7 0 2 14 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring5N1': ['1 5 5 6 0 2 14 6 0 2 14 6 0 2 14 6 0 2 14 7 0 2 14 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'ring5N2': ['1 5 5 6 0 2 14 6 0 2 14 7 0 2 14 6 0 2 14 7 0 2 14 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'iPr': ['1 3 2 6 3 4 14 6 1 4 14 6 3 4 14 1 2 2 3', 'Substructure2'],
    'tBu': ['1 4 3 6 3 4 14 6 0 4 4 6 3 4 14 6 3 4 14 1 2 2 3 2 4', 'Substructure2'],
    'C3N': ['1 4 3 6 0 1 14 7 0 3 14 6 0 1 14 6 0 1 14 1 2 2 3 2 4', 'Substructure2'],
    'C2O': ['1 3 2 6 0 1 14 8 0 2 14 6 0 1 14 1 2 2 3', 'Substructure2'],
    'CNO': ['1 3 2 6 0 1 14 7 0 2 14 8 0 1 14 1 2 2 3', 'Substructure2'],
    'C3P': ['1 4 3 6 0 1 14 15 0 3 14 6 0 1 14 6 0 1 14 1 2 2 3 2 4', 'Substructure2'],
    'CO2': ['1 3 2 8 0 1 14 6 0 2 14 8 0 2 14 1 2 2 3', 'Substructure2'],
    # 'C_Hal': ['', 'Substructure2'],
}

SET_ELEMENTS = {
    'hetero': ['B', 'Si', 'P', 'S', 'As', 'Se', 'Ge', 'Sb', 'Te'],
    'other_met': ['Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po'],
    'd_Me4': ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    'd_Me5': ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'],
    'd_Me67': ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Rf', 'Db', 'Sg', 'Bh',
               'Hs', 'Mt'],
    'f_Me': ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
             'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
             'Fm', 'Md', 'No', 'Lr'],
    'AEMet': ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'],
    'AMet': ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'],
    'halogens': ['F', 'Cl', 'Br', 'I', 'At']
}

ELEMENTS_SET_CLASSES = [
    ElementsSet1, ElementsSet2, ElementsSet3, ElementsSet4,
    ElementsSet5, ElementsSet6, ElementsSet7, ElementsSet8
]


def get_fields_list(model):
    return [field.name for field in model._meta.get_fields()]


def reset_old_data_or_create(refcode):
    # reset substructure1 model objects
    fields_substr1 = get_fields_list(Substructure1)
    fields_substr1.remove('id')
    fields_substr1.remove('refcode')
    set_substr1 = dict.fromkeys(fields_substr1, False)
    Substructure1.objects.update_or_create(refcode__refcode=refcode, defaults=set_substr1)
    # reset substructure2 model objects
    fields_substr2 = get_fields_list(Substructure2)
    fields_substr2.remove('id')
    fields_substr2.remove('refcode')
    set_substr2 = dict.fromkeys(fields_substr2, False)
    Substructure2.objects.update_or_create(refcode__refcode=refcode, defaults=set_substr2)


def set_only_CHNO(graphs, substructure_obj=Substructure1, filter_template='refcode__elements__element_set'):
    filtr = {}
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element not in ['id', 'C', 'H', 'N', 'O'] and len(element) <= 3:
                filtr[f'{filter_template}_{i}__isnull'] = True
                filtr[f'{filter_template}_{i}__{element.capitalize()}__isnull'] = True
    structures = graphs.filter(**filtr)
    for structure in structures:
        substr, created = substructure_obj.objects.get_or_create(refcode=structure.refcode)
        substr.only_CHNO = True
        substr.save()


def set_no_C(graphs, substructure_obj=Substructure1, filter_template='refcode__elements__element_set'):
    filtr = {}
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element == 'C':
                filtr[f'{filter_template}_{i}__isnull'] = True
                filtr[f'{filter_template}_{i}__{element.capitalize()}__isnull'] = True
                break
    structures = graphs.filter(**filtr)
    for structure in structures:
        substr, created = substructure_obj.objects.get_or_create(refcode=structure.refcode)
        substr.no_C = True
        substr.save()


def set_elements(
        graphs, attr_name, element_set,
        substructure_obj=Substructure1,
        filter_template='refcode__elements__element_set'
):
    filtr = []
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element in element_set:
                temp = dict()
                temp[f'{filter_template}_{i}__isnull'] = False
                temp[f'{filter_template}_{i}__{element.capitalize()}__isnull'] = False
                filtr.append(Q(**temp))
    custom_filter = filtr[0]
    for i in range(1, len(element_set), 1):
        custom_filter = custom_filter | filtr[i]
    structures = graphs.filter(custom_filter)
    for structure in structures:
        substr, created = substructure_obj.objects.get_or_create(refcode=structure.refcode)
        setattr(substr, attr_name, True)
        substr.save()


def start_dll_and_write(
        template_graph, analyse_data,
        NUM_OF_PROC, attr_name,
        obj_to_save, struct_obj=StructureCode
):
    output = cpplib.SearchMain(template_graph, list(analyse_data), NUM_OF_PROC, False)
    for elem in output:
        refcode = struct_obj.objects.get(id=elem)
        substr, created = obj_to_save.objects.get_or_create(refcode=refcode)
        setattr(substr, attr_name, True)
        substr.save()


def add_substructure_filters(refcodes, NUM_OF_PROC=1):
    substructure_logger.info('Getting structures...')
    structures = CoordinatesBlock.objects.filter(refcode__refcode__in=refcodes)

    substructure_logger.info('Reset old data or create new...')
    for structure in structures:
        reset_old_data_or_create(structure.refcode)

    substructure_logger.info('Add substructure filtration...')
    models = {'Substructure1': Substructure1,
              'Substructure2': Substructure2}
    analyse_data = structures.values_list('graph', flat=True)
    if analyse_data:
        for attr_name, data in TEMPLATES.items():
            template_graph, obj_name = data
            start_dll_and_write(template_graph, analyse_data, NUM_OF_PROC, attr_name, models[obj_name])

    substructure_logger.info('Add element filtration...')
    # set_only_CHNO(structures)
    # set_no_C(structures)
    for attr_name, element_set in SET_ELEMENTS.items():
        set_elements(structures, attr_name, element_set)

    substructure_logger.info('Success!')
