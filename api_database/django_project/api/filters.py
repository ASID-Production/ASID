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

from django_filters.rest_framework import FilterSet, filters
from structure.models import (StructureCode, ElementsSet1, ElementsSet2,
                              ElementsSet3, ElementsSet4, ElementsSet5,
                              ElementsSet6, ElementsSet7, ElementsSet8,
                              CENTRINGS)
from gemmi import UnitCell, SpaceGroup, GruberVector


def get_fields_list(model):
    return [field.name for field in model._meta.get_fields()]


def split_element_and_count(element):
    for i in range(len(element)):
        if i < len(element) - 1 and element[i].isalpha() and element[i + 1].isdigit():
            elem = element[:i + 1]
            count = element[i + 1:]
            return elem, count
    return 'error', '404'


def get_reduced_cell(params: list, centring: str) -> list:
    cell = UnitCell(params[0], params[1], params[2], params[3], params[4], params[5])
    if centring.upper() != 'R':
        sg = SpaceGroup(centring.upper() + '1')
    else:
        sg = SpaceGroup(centring.upper() + '3')
    gv = GruberVector(cell, sg)
    gv.niggli_reduce()
    reduced_params = gv.cell_parameters()
    return list(reduced_params)


class StructureFilter(FilterSet):
    # refcode search
    refcode = filters.CharFilter(method='filter_refcode')
    CCDC_number = filters.CharFilter(method='filter_CCDC_number')
    # database search
    CSD = filters.CharFilter(method='filter_CSD')
    ICSD = filters.CharFilter(method='filter_ICSD')
    COD = filters.CharFilter(method='filter_COD')
    user_db = filters.CharFilter(method='filter_user_db')
    # name search
    name = filters.CharFilter(method='filter_name')
    # formula search
    formula = filters.CharFilter(method='filter_formula')
    # element search
    elements = filters.CharFilter(method='filter_elements')
    # element exclude search
    no_elements = filters.CharFilter(method='filter_no_elements')
    # doi search (doi format: 10.7503/cjcu20200475)
    doi = filters.CharFilter(field_name='publication__publication__doi', lookup_expr='iexact')
    # temperature search
    temperature = filters.NumberFilter(
        field_name='experimental_info__structure_determination_temperature',
        lookup_expr='exact'
    )
    # temperature range search
    temperature_range = filters.RangeFilter(field_name='experimental_info__structure_determination_temperature')
    # authors search
    authors = filters.CharFilter(method='filter_authors')
    # cell search (format: a,b,c,al,be,ga,centring,parameter_deviation(value_or_none),angle_deviation(value_or_none); for example: 6,7,8,90,90,90,P,0.015,0.02)
    cell = filters.CharFilter(method='filter_cell')

    class Meta:
        model = StructureCode
        fields = []

    def filter_refcode(self, queryset, name, value):
        if value:
            queryset_temp = queryset
            for element in value.split():
                queryset_temp = queryset.filter(refcode__istartswith=element)
            return queryset_temp
        return queryset

    def filter_CCDC_number(self, queryset, name, value):
        if value:
            for element in value.split():
                queryset = queryset.filter(CCDC_number__icontains=element)
        return queryset

    def filter_ICSD(self, queryset, name, value):
        if value:
            queryset_temp = queryset
            if value == 'False':
                queryset_temp = queryset.filter(ICSD__isnull=True)
            return queryset_temp
        return queryset

    def filter_CSD(self, queryset, name, value):
        if value:
            queryset_temp = queryset
            if value == 'False':
                queryset_temp = queryset.filter(CCDC_number__isnull=True)
            return queryset_temp
        return queryset

    def filter_COD(self, queryset, name, value):
        if value:
            queryset_temp = queryset
            if value == 'False':
                queryset_temp = queryset.filter(COD__isnull=True)
            return queryset_temp
        return queryset

    def filter_user_db(self, queryset, name, value):
        if value:
            queryset_temp = queryset
            if value == 'False':
                queryset_temp = queryset.filter(user__isnull=True)
            return queryset_temp
        return queryset

    def filter_name(self, queryset, name, value):
        if value:
            for element in value.split():
                queryset = queryset.filter(name__systematic_name__icontains=element) | queryset.filter(name__trivial_name__icontains=element)
        return queryset

    def filter_formula(self, queryset, name, value):
        if value:
            filtr = {}
            # for each element of user defined
            for element in value.split():
                # define element and it's quantity
                elem, count = split_element_and_count(element)
                # if there is a problem with spliting, return empty queryset
                if elem == 'error':
                    return queryset.none()
                # find database table with current element
                for i, model in enumerate(
                        [ElementsSet1, ElementsSet2,
                        ElementsSet3, ElementsSet4, ElementsSet5,
                        ElementsSet6, ElementsSet7, ElementsSet8], start=1):
                    # add filter, the quality of the element
                    if elem.capitalize() in get_fields_list(model):
                        filtr[f'element_set_{i}__{elem.capitalize()}__exact'] = float(count)
                        break
                # if element was not found, return empty queryset
                else:
                    return queryset.none()
            return queryset.filter(**filtr)
        return queryset

    def filter_elements(self, queryset, name, value, exclude=False):
        if value:
            filtr = {}
            # for each element of user defined
            for element in value.split():
                # find database table with current element
                for i, model in enumerate(
                        [ElementsSet1, ElementsSet2,
                        ElementsSet3, ElementsSet4, ElementsSet5,
                        ElementsSet6, ElementsSet7, ElementsSet8], start=1):
                    # add filter, that element must be not Null
                    if element.capitalize() in get_fields_list(model):
                        if not exclude:
                            filtr[f'element_set_{i}__{element.capitalize()}__isnull'] = False
                        else:
                            filtr[f'element_set_{i}__{element.capitalize()}__isnull'] = True
                        break
                # if element was not found, return empty queryset
                else:
                    return queryset.none()
            return queryset.filter(**filtr)
        return queryset

    def filter_no_elements(self, queryset, name, value):
        return self.filter_elements(queryset, name, value, True)

    def filter_authors(self, queryset, name, value):
        if value:
            for element in value.split():
                queryset = queryset.filter(authors__family_name__icontains=element)
        return queryset

    def filter_cell(self, queryset, name, value):
        if value:
            raw_params = value.split(',')
            params = list(map(float, raw_params[:6]))

            # get reduced cell
            params = get_reduced_cell(params, raw_params[6].upper())

            params.extend(raw_params[6:])
            abc_diviation = 0.015
            angle_diviation = 0.02
            centrings = dict((v, k) for k, v in CENTRINGS)
            # if diviation is defined by user
            if params[7] != 'none':
                abc_diviation = float(params[7])
            if params[8] != 'none':
                angle_diviation = float(params[8])
            queryset = queryset.filter(
                reduced_cells__a__range=(params[0] - params[0] * abc_diviation, params[0] + params[0] * angle_diviation),
                reduced_cells__b__range=(params[1] - params[1] * abc_diviation, params[1] + params[1] * angle_diviation),
                reduced_cells__c__range=(params[2] - params[2] * abc_diviation, params[2] + params[2] * angle_diviation),
                reduced_cells__al__range=(params[3] - params[3] * abc_diviation, params[3] + params[3] * angle_diviation),
                reduced_cells__be__range=(params[4] - params[4] * abc_diviation, params[4] + params[4] * angle_diviation),
                reduced_cells__ga__range=(params[5] - params[5] * abc_diviation, params[5] + params[5] * angle_diviation),
                cell__centring__exact=centrings[params[6].upper()]
            )
            return queryset
        return queryset
