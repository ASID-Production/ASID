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

from rest_framework.viewsets import ReadOnlyModelViewSet
from structure.models import StructureCode, CifFile, CoordinatesBlock
from structure.download import create_cif_text
from .serializers import RefcodeShortSerializer, RefcodeFullSerializer, CifUploadSerializer, SearchSerializer
from django.shortcuts import get_object_or_404
from rest_framework import status
from django.http import HttpResponse
from django_filters.rest_framework import DjangoFilterBackend
from .filters import StructureFilter
from .substructure_filtration import set_filter
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
import networkx as nx
from structure.management.commands.cif_db_update_modules._element_numbers import element_numbers
from multiprocessing import cpu_count
import ctypes
from django.conf import settings
from itertools import chain
from structure.management.commands.cif_db_update import main as add_cif_data
import os

NUM_OF_PROC = int(cpu_count() / 2)
MAX_STRS_SIZE = 30000


class StructureViewSet(ReadOnlyModelViewSet):
    queryset = StructureCode.objects.all()
    filter_backends = (DjangoFilterBackend,)
    filterset_class = StructureFilter

    def get_serializer_class(self):
        if self.action == 'list':
            return RefcodeShortSerializer
        return RefcodeFullSerializer

    @action(
        detail=True,
        methods=['GET']
    )
    def download(self, request, pk):
        structure = get_object_or_404(StructureCode, pk=pk)
        filename = f'{structure.refcode}.cif'
        content = create_cif_text(structure)
        response = HttpResponse(content, content_type='text/plain')
        response['Content-Disposition'] = f'attachment; filename={filename}'
        return response

    @action(
        detail=False,
        methods=['POST'],
        permission_classes=[IsAuthenticated],
        serializer_class=[CifUploadSerializer],
    )
    def upload(self, request):
        user = request.user
        count_user_cifs = StructureCode.objects.filter(user=user).count()
        # if user deleted structures before
        while True:
            refcode = 'user-' + str(user.id) + '-' + str(count_user_cifs + 1)
            if StructureCode.objects.filter(refcode=refcode).count() == 0:
                break
            count_user_cifs += 1
        file = request.FILES.get('file')
        # if object is not exists, create new structure object with defined refcode
        if not CifFile.objects.filter(refcode__user=user).filter(old_file_name=file.name).exists():
            refcode_obj = StructureCode.objects.create(user=user, refcode=refcode)
            cif_file_obj = CifFile.objects.create(refcode=refcode_obj, file=file, old_file_name=file.name)
        else:
            cif_file_obj = CifFile.objects.get(refcode__user=user, old_file_name=file.name)
            refcode = cif_file_obj.refcode.refcode
        cif_file_path = os.path.join(settings.BASE_DIR, 'media', str(cif_file_obj.file))
        try:
            add_cif_data(
                args=(cif_file_path,),
                all_data=True,
                user_refcodes={cif_file_path: refcode, }
            )
        except Exception as error_message:
            return Response(
                {'errors': f'Structure information was not added! {error_message}'},
                status=status.HTTP_400_BAD_REQUEST
            )
        return Response(status=status.HTTP_201_CREATED)

    @action(
        detail=False,
        methods=['GET'],
    )
    def search(self, request):
        serializer = SearchSerializer(data=request.data)
        if not serializer.is_valid():
            return Response(
                serializer.errors,
                status=status.HTTP_400_BAD_REQUEST
            )
        graph = nx.Graph()
        graph.add_nodes_from(serializer.data.get('nodes'))
        graph.add_edges_from(serializer.data.get('edges'))
        template_data = get_template_graph(graph)
        analyse_data, size = get_search_queryset()

        array_template = ctypes.c_char_p(template_data.encode())
        array_data = (ctypes.c_char_p * size)(*[s.encode() for s in analyse_data])
        dll = settings.GET_DLL()
        if serializer.data.get('search_type') == 'substructure':
            output = dll.SearchMain(array_template, array_data, size, NUM_OF_PROC, False)
        elif serializer.data.get('search_type') == 'exact':
            output = dll.SearchMain(array_template, array_data, size, NUM_OF_PROC, True)
        else:
            return Response(
                {'errors': 'Unsupported value of "search_type" parameter. Use exact or substructure keywords'},
                status=status.HTTP_400_BAD_REQUEST
            )

        out_refcode_ids = []
        for i in range(1, output[0] + 1):
            out_refcode_ids.append(output[i])
        refcodes = StructureCode.objects.none()
        if len(out_refcode_ids) > MAX_STRS_SIZE:
            for i in range(0, len(out_refcode_ids), MAX_STRS_SIZE):
                qset = StructureCode.objects.filter(id__in=out_refcode_ids[i:i + MAX_STRS_SIZE])
                refcodes = sorted(
                    chain(refcodes, qset),
                    key=lambda structure: structure.refcode, reverse=False)

        else:
            refcodes = StructureCode.objects.filter(id__in=out_refcode_ids).order_by('refcode')
        pages = self.paginate_queryset(refcodes)
        out_serializer = RefcodeShortSerializer(pages, many=True)
        return self.get_paginated_response(out_serializer.data)


def get_search_queryset():
    graphs = CoordinatesBlock.objects.filter(graph__isnull=False)
    analyse_data = graphs.values_list('graph', flat=True)
    return analyse_data, graphs.count()


def get_search_queryset_with_filtration(template):
    graphs = CoordinatesBlock.objects.filter(graph__isnull=False)
    filters = set_filter(template)
    filtrs = dict()
    for filtr, data in filters.items():
        is_true, obj_name = data
        if is_true:
            filtrs[f'refcode__{obj_name.lower()}__{filtr}'] = is_true
    structures = graphs.filter(**filtrs)
    analyse_data = structures.values_list('graph', flat=True)
    return analyse_data, structures.count()


def get_template_graph(graph):
    result = []
    nodes_order = []
    result.append(1)
    result.append(graph.number_of_nodes())
    result.append(graph.number_of_edges())
    # add atoms
    for atom_idx, atom_type in nx.get_node_attributes(graph, "type").items():
        result.append(element_numbers[atom_type])
        result.append(graph.nodes[atom_idx]['Hnum'])
        nodes_order.append(atom_idx)
    # add bonds
    for atom1, atom2 in graph.edges:
        order_id_1 = nodes_order.index(atom1) + 1
        order_id_2 = nodes_order.index(atom2) + 1
        result.append(order_id_1)
        result.append(order_id_2)
    return ' '.join(map(str, result))
