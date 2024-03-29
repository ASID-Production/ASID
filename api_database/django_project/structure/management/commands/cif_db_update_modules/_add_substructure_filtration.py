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
import ctypes
from django.conf import settings
from django_project.loggers import substructure_logger


TEMPLATES = {
    'NO2': ['1 3 2 8 0 7 0 8 0 1 2 2 3', 'Substructure1'],
    'SO2': ['1 3 2 8 0 16 0 8 0 1 2 2 3', 'Substructure1'],
    # 'C_Met': ['', 'Substructure1'],
    # 'N_Met': ['', 'Substructure1'],
    'CS': ['1 2 1 16 0 6 0 1 2', 'Substructure1'],
    'ring6': ['1 6 6 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring5': ['1 5 5 6 0 6 0 6 0 6 0 6 0 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'Ph': ['1 6 6 6 1 6 1 6 1 6 0 6 1 6 1 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring6N1': ['1 6 6 6 0 6 0 6 0 6 0 6 0 7 0 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring6N2': ['1 6 6 6 0 6 0 7 0 6 0 6 0 7 0 1 2 1 6 2 3 3 4 4 5 5 6', 'Substructure2'],
    'ring5N1': ['1 5 5 6 0 6 0 6 0 6 0 7 0 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'ring5N2': ['1 5 5 6 0 6 0 7 0 6 0 7 0 1 2 1 5 2 3 3 4 4 5', 'Substructure2'],
    'iPr': ['1 3 2 6 3 6 0 6 3 1 2 2 3', 'Substructure2'],
    'tBu': ['1 4 3 6 3 6 0 6 3 6 3 1 2 2 3 2 4', 'Substructure2'],
    'C3N': ['1 4 3 6 0 7 0 6 0 6 0 1 2 2 3 2 4', 'Substructure2'],
    'C2O': ['1 3 2 6 0 8 0 6 0 1 2 2 3', 'Substructure2'],
    'CNO': ['1 3 2 6 0 7 0 8 0 1 2 2 3', 'Substructure2'],
    'C3P': ['1 4 3 6 0 15 0 6 0 6 0 1 2 2 3 2 4', 'Substructure2'],
    'CO2': ['1 3 2 8 0 6 0 8 0 1 2 2 3', 'Substructure2'],
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


def dell_old_data(refcode):
    filters1 = Substructure1.objects.filter(refcode__refcode=refcode)
    filters2 = Substructure2.objects.filter(refcode__refcode=refcode)
    filters1.delete()
    filters2.delete()


def set_only_CHNO(graphs):
    filtr = {}
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element not in ['refcode', 'id', 'C', 'H', 'N', 'O']:
                filtr[f'refcode__element_set_{i}__{element.capitalize()}__isnull'] = True
    structures = graphs.filter(**filtr)
    for structure in structures:
        substr, created = Substructure1.objects.get_or_create(refcode=structure.refcode)
        substr.only_CHNO = True
        substr.save()


def set_no_C(graphs):
    filtr = {}
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element not in ['refcode', 'id'] and element == 'C':
                filtr[f'refcode__element_set_{i}__{element.capitalize()}__isnull'] = True
                break
    structures = graphs.filter(**filtr)
    for structure in structures:
        substr, created = Substructure1.objects.get_or_create(refcode=structure.refcode)
        substr.no_C = True
        substr.save()


def set_elements(graphs, attr_name, element_set):
    filtr = []
    for i, model in enumerate(ELEMENTS_SET_CLASSES, start=1):
        for element in get_fields_list(model):
            if element not in ['refcode', 'id'] and element in element_set:
                temp = dict()
                temp[f'refcode__element_set_{i}__{element.capitalize()}__isnull'] = False
                filtr.append(Q(**temp))
    custom_filter = filtr[0]
    for i in range(1, len(element_set), 1):
        custom_filter = custom_filter | filtr[i]
    structures = graphs.filter(custom_filter)
    for structure in structures:
        substr, created = Substructure1.objects.get_or_create(refcode=structure.refcode)
        setattr(substr, attr_name, True)
        substr.save()


def start_dll_and_write(template_graph, analyse_data_c, size, NUM_OF_PROC, attr_name, obj_name):
    template = ctypes.c_char_p(template_graph.encode())
    dll = settings.GET_DLL()
    output = dll.SearchMain(template, analyse_data_c, size, NUM_OF_PROC, False)
    for i in range(1, output[0] + 1):
        refcode = StructureCode.objects.get(id=output[i])
        if obj_name == 'Substructure1':
            substr, created = Substructure1.objects.get_or_create(refcode=refcode)
        elif obj_name == 'Substructure2':
            substr, created = Substructure2.objects.get_or_create(refcode=refcode)
        else:
            raise Exception('Invalid obj_name value!')
        setattr(substr, attr_name, True)
        substr.save()


def add_substructure_filters(refcodes, NUM_OF_PROC=1):
    substructure_logger.info('Getting structures...')
    structures = CoordinatesBlock.objects.filter(refcode__refcode__in=refcodes)

    substructure_logger.info('Delete old data...')
    for structure in structures:
        dell_old_data(structure.refcode)
        Substructure1.objects.get_or_create(refcode=structure.refcode)
        Substructure2.objects.get_or_create(refcode=structure.refcode)

    substructure_logger.info('Add substructure filtration...')
    size = structures.count()
    analyse_data = structures.values_list('graph', flat=True)
    if analyse_data:
        analyse_data_c = (ctypes.c_char_p * size)(*[s.encode() for s in analyse_data])
        for attr_name, data in TEMPLATES.items():
            template_graph, obj_name = data
            start_dll_and_write(template_graph, analyse_data_c, size, NUM_OF_PROC, attr_name, obj_name)

    substructure_logger.info('Add element filtration...')
    set_only_CHNO(structures)
    set_no_C(structures)
    for attr_name, element_set in SET_ELEMENTS.items():
        set_elements(structures, attr_name, element_set)

    substructure_logger.info('Success!')
