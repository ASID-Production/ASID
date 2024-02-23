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
from itertools import chain, zip_longest
from structure.management.commands.cif_db_update import main as add_cif_data
import os
from asgiref.sync import sync_to_async
from .pagination import LimitPagination
from rest_framework.views import APIView
from rest_framework.renderers import JSONRenderer
from django.core.cache import cache

NUM_OF_PROC = int(cpu_count() / 2)
MAX_STRS_SIZE = 30000
CHUNK_SIZE = 10000  # the number of structures for search in


async def structure_search_view(request):
    view = APIView()
    request = view.initialize_request(request)
    serializer = SearchSerializer(data=request.data)
    if not serializer.is_valid():
        return HttpResponse(
            serializer.errors,
            status=status.HTTP_400_BAD_REQUEST
        )
    response = await structure_search(request, serializer)
    return response


@sync_to_async(thread_sensitive=False)
def structure_search(request, serializer):
    chunk_size = serializer.data.get('chunk_size')
    iter_num = serializer.data.get('iter_num')
    search_type = serializer.data.get('search_type')
    if chunk_size:
        global CHUNK_SIZE
        CHUNK_SIZE = chunk_size

    graph = nx.Graph()
    graph.add_nodes_from(serializer.data.get('nodes'))
    graph.add_edges_from(serializer.data.get('edges'))
    template_data = get_template_graph(graph)
    array_template = ctypes.c_char_p(template_data.encode())
    dll = settings.GET_DLL()

    if search_type == 'substructure':
        exact = False
    elif search_type == 'exact':
        exact = True
    else:
        return Response(
            {'errors': 'Unsupported value of "search_type" parameter. Use exact or substructure keywords'},
            status=status.HTTP_400_BAD_REQUEST
        )
    analyse_data_split = list()
    if chunk_size and request.user.is_authenticated:
        data = cache.get(f'{request.user}-cache')
        if (
                data and
                data['template_data'] == template_data and
                data['search_type'] == search_type and
                data['chunk_size'] == chunk_size
        ):
            analyse_data_split = data['analyse_data_split']
    if not analyse_data_split:
        analyse_data, size = get_search_queryset_with_filtration(template_data)
        # split search for CHUNK_SIZE structures parts and then merge the result
        analyse_data_split = list(zip_longest(*[iter(analyse_data)] * CHUNK_SIZE, fillvalue=''))
        cache.set(f'{request.user}-cache', {
            'analyse_data_split': analyse_data_split,
            'search_type': search_type,
            'template_data': template_data,
            'chunk_size': chunk_size
        })
    out_refcode_ids = []
    partial = False
    for analyse_data in analyse_data_split:
        # if partial return is needed
        if chunk_size and iter_num < len(analyse_data_split):
            partial = True
            analyse_data = analyse_data_split[iter_num]
        array_data = (ctypes.c_char_p * CHUNK_SIZE)(*[s.encode() for s in analyse_data])
        output = dll.SearchMain(array_template, array_data, CHUNK_SIZE, NUM_OF_PROC, exact)
        for j in range(1, output[0] + 1):
            out_refcode_ids.append(output[j])
        # break if partial
        if partial:
            break
    # get queryset on search result
    refcodes = StructureCode.objects.none()
    if len(out_refcode_ids) > MAX_STRS_SIZE:
        for i in range(0, len(out_refcode_ids), MAX_STRS_SIZE):
            qset = StructureCode.objects.filter(id__in=out_refcode_ids[i:i + MAX_STRS_SIZE])
            refcodes = sorted(
                chain(refcodes, qset),
                key=lambda structure: structure.refcode, reverse=False)
    else:
        refcodes = StructureCode.objects.filter(id__in=out_refcode_ids).order_by('refcode')
    paginator = LimitPagination()
    pages = paginator.paginate_queryset(refcodes, request)
    out_serializer = RefcodeShortSerializer(pages, many=True)
    response = paginator.get_paginated_response(out_serializer.data)
    response.accepted_renderer = JSONRenderer()
    response.accepted_media_type = "application/json"
    if partial:
        response.data.update({'max_iter_num': len(analyse_data_split) - 1})
        response.data.move_to_end('max_iter_num', last=False)
    response.renderer_context = {}
    response.render()
    return response


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
        refcode_obj = StructureCode.objects.create(user=user, refcode=refcode)
        cif_file_obj = CifFile.objects.create(refcode=refcode_obj, file=file, old_file_name=file.name)
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
    multitype_atoms = dict()
    multitype = -1
    result.append(1)
    result.append(graph.number_of_nodes())
    result.append(graph.number_of_edges())
    # add atoms
    for atom_idx, atom_type in nx.get_node_attributes(graph, "type").items():
        if len(atom_type.split()) > 1:
            multitype_atoms[multitype] = atom_type.split()
            result.append(multitype)
            multitype -= 1
        else:
            result.append(element_numbers[atom_type])
        result.append(graph.nodes[atom_idx]['Hnum'])
        nodes_order.append(atom_idx)
    # add bonds
    for atom1, atom2 in graph.edges:
        order_id_1 = nodes_order.index(atom1) + 1
        order_id_2 = nodes_order.index(atom2) + 1
        result.append(order_id_1)
        result.append(order_id_2)
    # if multitypes
    if multitype_atoms:
        for multitype_atom, atoms in multitype_atoms.items():
            result.append(multitype_atom)
            for atom in atoms:
                result.append(element_numbers[atom])
        result.append(0)
    return ' '.join(map(str, result))
