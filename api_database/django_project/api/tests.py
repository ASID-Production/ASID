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

from django.test import TestCase
from structure.models import StructureCode
from .views import get_search_queryset_with_filtration, start_SearchMain, get_queryset_from_ids
from structure.management.commands.cif_db_update_modules._add_substructure_filtration import TEMPLATES
from itertools import zip_longest
from progress.bar import IncrementalBar


def get_analysed_data(template_graph, queryset, chunk_size, enable_substr_filtr, enable_elem_filtr):
    analyse_data = get_search_queryset_with_filtration(
        template_graph,
        queryset,
        enable_substr_filtr=enable_substr_filtr,
        enable_elem_filtr=enable_elem_filtr
    )
    analyse_data_split = list(zip_longest(*[iter(analyse_data)] * chunk_size, fillvalue=''))
    return analyse_data_split


class Test(TestCase):
    def test_substructure(self):
        only_substructure = False
        only_elements = True
        chunk_size = 10000
        search_sets = {
            'without_filtration': '',
            'substructure_filtration': '',
            'elements_filtration': ''
        }
        queryset = StructureCode.objects.all()
        print(f'Search in {queryset.count()} structures')
        bar = IncrementalBar('Processing substructure tests', max=len(TEMPLATES))
        for group, template_graph in TEMPLATES.items():
            template_graph = template_graph[0]
            # data without filtration
            analyse_data = get_analysed_data(template_graph, queryset, chunk_size, False, False)
            out_refcode_ids = start_SearchMain(
                template_graph, analyse_data, False, 0, False
            )
            refcodes = get_queryset_from_ids(out_refcode_ids, StructureCode)
            search_sets['without_filtration'] = set(refcodes)

            # data with substructure filtration
            if not only_elements:
                analyse_data = get_analysed_data(template_graph, queryset, chunk_size, True, False)
                out_refcode_ids = start_SearchMain(
                    template_graph, analyse_data, False, 0, False
                )
                refcodes = get_queryset_from_ids(out_refcode_ids, StructureCode)
                search_sets['substructure_filtration'] = set(refcodes)
                diff_substr = search_sets['without_filtration'] - search_sets['substructure_filtration']
                # compare
                self.assertIs(
                    len(diff_substr),
                    0,
                    f'''SubstructureFiltrationError: Search of {group} fragment
without filtration has found {len(diff_substr)} structures more,
than search with substructure filtration!'''
                )

            # data with element filtration
            if not only_substructure:
                analyse_data = get_analysed_data(template_graph, queryset, chunk_size, False, True)
                out_refcode_ids = start_SearchMain(
                    template_graph, analyse_data, False, 0, False
                )
                refcodes = get_queryset_from_ids(out_refcode_ids, StructureCode)
                search_sets['elements_filtration'] = set(refcodes)
                diff_elem = search_sets['without_filtration'] - search_sets['elements_filtration']
                # compare results
                self.assertLess(
                    len(diff_elem),
                    len(search_sets['without_filtration']) * 0.0001,
                    f'''ElementsFiltrationError: Search of {group} fragment 
without filtration has found {len(diff_elem)} structures more, 
than search with elements filtration!'''
                )
            bar.next()
        bar.finish()
