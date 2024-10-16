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
from qc_structure.models import QCStructureCode, VaspFile, QCCoordinatesBlock
from qc_structure.vasp import vasp_parser as add_vasp_data
from qc_structure.vasp import get_or_create_space_group as vasp_get_or_create_space_group
from qc_structure.vasp import save_coordinates as vasp_save_coordinates
from qc_structure.export.cif import qc_get_cif_content
from .serializers import (RefcodeShortSerializer, RefcodeFullSerializer, CifUploadSerializer,
                          SearchSerializer, QCRefcodeShortSerializer, QCRefcodeFullSerializer,
                          VaspUploadSerializer, Gen2DImgSerializer)
from django.shortcuts import get_object_or_404
from rest_framework import status
from django.http import HttpResponse
from django_filters.rest_framework import DjangoFilterBackend
from .filters import StructureFilter, QCStructureFilter
from .substructure_filtration import set_filter
from .viewsets import StructureModelViewSet
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
import networkx as nx
from structure.management.commands.cif_db_update_modules._element_numbers import element_numbers
from multiprocessing import cpu_count
from django.conf import settings
from itertools import chain, zip_longest
from structure.management.commands.cif_db_update import main as add_cif_data
import os
from asgiref.sync import sync_to_async
from .pagination import LimitPagination
from rest_framework.views import APIView
from rest_framework.renderers import JSONRenderer
from django.core.cache import cache
import base64
from io import BytesIO
import cpplib
from rest_framework.decorators import api_view
from structure.management.commands.cif_db_update import collect_cif_data
from structure.management.commands.cif_db_update_modules._make_graphs_c import get_data
from structure.management.commands.cif_db_update_modules._add_all_cif_data import get_or_create_space_group
from modules.gen2d.gen2d import main_v2 as gen_smiles_inchi
from modules.gen2d.gen2d import gen2d
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from django.core.files.uploadedfile import TemporaryUploadedFile
from pymatgen.core import Molecule
import itertools

NUM_OF_PROC = int(cpu_count() / 2)
MAX_STRS_SIZE = 30000
CHUNK_SIZE = 10000  # the number of structures for search in
R_H_DIST = 0.95  # distance of R-H bonds in angstroms


async def qc_structure_search_view(request):
    view = APIView()
    request = view.initialize_request(request)
    serializer = SearchSerializer(data=request.data)
    if not serializer.is_valid():
        return Response(
            serializer.errors,
            status=status.HTTP_400_BAD_REQUEST
        )
    response = await structure_search(
        request,
        serializer,
        out_serializer_model=QCRefcodeShortSerializer,
        structure_code_model=QCStructureCode,
        qc=True
    )
    return response


async def structure_search_view(request):
    view = APIView()
    request = view.initialize_request(request)
    serializer = SearchSerializer(data=request.data)
    if not serializer.is_valid():
        return Response(
            serializer.errors,
            status=status.HTTP_400_BAD_REQUEST
        )
    response = await structure_search(request, serializer)
    return response


@sync_to_async(thread_sensitive=False)
def structure_search(
        request,
        serializer,
        qc=False,
        out_serializer_model=RefcodeShortSerializer,
        structure_code_model=StructureCode
):
    chunk_size = serializer.data.get('chunk_size')
    iter_num = serializer.data.get('iter_num')
    search_type = serializer.data.get('search_type')
    partial = True
    if not chunk_size:
        partial = False
        chunk_size = CHUNK_SIZE
    graph = nx.Graph()
    graph.add_nodes_from(serializer.data.get('nodes'))
    graph.add_edges_from(serializer.data.get('edges'))
    template_data = get_template_graph(graph)
    template_data_without_cord = get_template_graph(graph, cord_num=False)

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
    if partial and request.user.is_authenticated:
        data = cache.get(f'{request.user}-cache')
        if (
                data and
                data['template_data'] == template_data and
                data['search_type'] == search_type and
                data['chunk_size'] == chunk_size
        ):
            analyse_data_split = data['analyse_data_split']
    if not analyse_data_split:
        if not qc:
            analyse_data = get_search_queryset_with_filtration(template_data_without_cord)
        else:
            analyse_data = get_search_queryset_with_filtration_qc(template_data_without_cord)
        # split search for chunk_size structures parts and then merge the result
        analyse_data_split = list(zip_longest(*[iter(analyse_data)] * chunk_size, fillvalue=''))
        if partial and request.user.is_authenticated:
            cache.set(f'{request.user}-cache', {
                'analyse_data_split': analyse_data_split,
                'search_type': search_type,
                'template_data': template_data,
                'chunk_size': chunk_size
            })
    out_refcode_ids = []
    for analyse_data in analyse_data_split:
        # if partial return is needed
        if partial and iter_num < len(analyse_data_split):
            analyse_data = analyse_data_split[iter_num]
        output = cpplib.SearchMain(template_data, list(analyse_data), NUM_OF_PROC, exact)
        out_refcode_ids.extend(output)
        # break if partial
        if partial:
            break
    # get queryset on search result
    refcodes = structure_code_model.objects.none()
    if len(out_refcode_ids) > MAX_STRS_SIZE:
        for i in range(0, len(out_refcode_ids), MAX_STRS_SIZE):
            qset = structure_code_model.objects.filter(id__in=out_refcode_ids[i:i + MAX_STRS_SIZE])
            refcodes = sorted(
                chain(refcodes, qset),
                key=lambda structure: structure.refcode, reverse=False)
    else:
        refcodes = structure_code_model.objects.filter(id__in=out_refcode_ids).order_by('refcode')
    paginator = LimitPagination()
    pages = paginator.paginate_queryset(refcodes, request)
    out_serializer = out_serializer_model(pages, many=True)
    response = paginator.get_paginated_response(out_serializer.data)
    response.accepted_renderer = JSONRenderer()
    response.accepted_media_type = "application/json"
    if partial:
        response.data.update({'max_iter_num': len(analyse_data_split) - 1})
        response.data.move_to_end('max_iter_num', last=False)
    response.renderer_context = {}
    response.render()
    return response


def get_img2d(structure, request):
    '''Available formats: img (gif) and cml (ChemDraw).'''
    h = int(request.GET.get('h', 0))
    w = int(request.GET.get('w', 0))
    f = request.GET.get('f', 'img')
    file = int(request.GET.get('file', 0))
    if h and w:
        content = structure.gen_2d_img(size=(w, h), f=f)
    else:
        content = structure.gen_2d_img(f=f)
    if file and content:
        if f == 'img':
            content = BytesIO(base64.b64decode(content.split(',')[-1])).getvalue()
            f = 'gif'
        response = HttpResponse(content, content_type='text/plain', status=200)
        response['Content-Disposition'] = f'attachment; filename={structure.refcode}.{f}'
        return response
    else:
        return HttpResponse(content, content_type='text/plain', status=200)


@api_view(['GET'])
def gen_img2d_view(request):
    serializer = Gen2DImgSerializer(data=request.data)
    if serializer.is_valid():
        try:
            file = request.data.get('file')
            file_format = serializer.data.get('file_format')
            output_format = serializer.data.get('output_format')
            return_type = serializer.data.get('return_type')
            h_size = serializer.data.get('h_size')
            w_size = serializer.data.get('w_size')
            name = serializer.data.get('name')
            # generate smiles and inchi
            inchi = ''
            smiles = ''
            # cif
            if file_format == 'cif':
                data = collect_cif_data(file, dict(), user_refcode='1', use_db=False)
                cif_block = data['1']
                # INFO: "get_or_create_space_group" with "return_only_symops"=True does not change data in database!
                symops_db = get_or_create_space_group(cif_block[1], return_only_symops=True)
                params, coords_types, types, symops = get_data(cif_block, symops_db)
                cpplib_result = cpplib.FindMoleculesInCell(params, symops, coords_types)
            # vasp
            elif file_format == 'vasp':
                if type(file) is TemporaryUploadedFile:
                    vasp_out = Vasprun(file.temporary_file_path())
                else:
                    vasp_out = Vasprun(file.file.getvalue().decode())
                vasp_structure = vasp_out.final_structure
                # INFO: "vasp_get_or_create_space_group" with "return_only_symops"=True does not change data in database!
                symops = vasp_get_or_create_space_group(vasp_structure, return_only_symops=True).split(';')
                spgran = SpacegroupAnalyzer(vasp_structure, symprec=0.02)
                symmed_vasp_struct = spgran.get_refined_structure()
                a, b, c = symmed_vasp_struct.lattice.abc
                al, be, ga = symmed_vasp_struct.lattice.angles
                params = [a, b, c, al, be, ga]
                sites_info = vasp_save_coordinates(None, symmed_vasp_struct, return_only_str_sites=True).split('\n')
                types = []
                atoms_coords_types = []
                for site in sites_info:
                    if site:
                        site_info = site.split()
                        element = site_info[1]
                        types.append(element_numbers[element])
                        coords = [element_numbers[element], float(site_info[2]), float(site_info[3]),
                                  float(site_info[4])]
                        atoms_coords_types.append(tuple(coords))
                cpplib_result = cpplib.FindMoleculesInCell(params, symops, atoms_coords_types)
            # xyz
            elif file_format == 'xyz':
                if type(file) is TemporaryUploadedFile:
                    file_content = file.temporary_file_path()
                    mol = Molecule.from_file(file_content)
                else:
                    file_content = file.file.getvalue().decode()
                    mol = Molecule.from_str(file_content, 'xyz')
                types = list(mol.atomic_numbers)
                coords = tuple(itertools.chain(*mol.cart_coords))
                coords = [coords[i:i + 3] for i in range(0, len(coords), 3)]
                coords_types = []
                for idx, pos in enumerate(coords):
                    coords_types.append(tuple([types[idx], *pos]))
                cpplib_result = cpplib.FindMoleculesWithoutCell(coords_types)

            # change R-H distance
            # sites = []
            # for molc in cpplib_result['xyz_block']:
            #     for atom in molc['atoms']:
            #         sites.append((types[atom['init_idx']], atom['x'], atom['y'], atom['z']))
            #     new_sites = cpplib.himp(sites, 0.95)['xyz_block']
            #     for idx, atom in enumerate(molc['atoms']):
            #         atom['x'], atom['y'], atom['z'] = new_sites[idx]

            # get smiles and inchi
            xyz_mols = cpplib_result['xyz_block']
            data_2d = gen_smiles_inchi(xyz_mols, element_numbers, types)
            if data_2d:
                smiles = data_2d['smiles']
                inchi = data_2d['inchi']

            # get response
            if not smiles and not inchi:
                return Response({'error': 'Failed to create smiles or inchi for this structure!'}, status=status.HTTP_400_BAD_REQUEST)
            # generate 2d representation
            img_type = output_format
            if output_format in ['gif', ]:
                img_type = 'img'
            content = gen2d(smiles=smiles, inchis=[inchi, ], sanitize=False, size=(w_size, h_size), format=img_type)
            if content:
                if img_type == 'img':
                    buffer = BytesIO()
                    content.save(buffer, output_format)
                # return as string
                if return_type == 'string':
                    if img_type == 'img':
                        content = f'data:image/{output_format};base64,' + base64.b64encode(buffer.getvalue()).decode()
                    response = HttpResponse(content, content_type='text/plain', status=200)
                # return as file
                elif return_type == 'file':
                    if img_type == 'img':
                        content = buffer.getvalue()
                    response = HttpResponse(content, content_type='text/plain', status=200)
                    response['Content-Disposition'] = f'attachment; filename={name}.{output_format}'
                return response
            return Response({'error': 'Failed to create 2D graph for this structure!'}, status=status.HTTP_400_BAD_REQUEST)
        except Exception as err:
            return Response({'error': str(err)}, status=status.HTTP_400_BAD_REQUEST)
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class StructureViewSet(StructureModelViewSet):
    filter_backends = (DjangoFilterBackend,)
    filterset_class = StructureFilter

    def get_queryset(self):
        if self.request.user.is_authenticated:
            queryset = (StructureCode.objects.filter(user=self.request.user) |
                        StructureCode.objects.filter(user__isnull=True))
            return queryset
        return StructureCode.objects.filter(user__isnull=True)

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
        detail=True,
        methods=['GET'],
        url_path='export/2d'
    )
    def export_2d(self, request, pk):
        structure = get_object_or_404(StructureCode, pk=pk)
        return get_img2d(structure, request)

    @action(
        detail=False,
        methods=['POST'],
        permission_classes=[IsAuthenticated],
        serializer_class=[CifUploadSerializer],
    )
    def upload(self, request):
        user = request.user
        count_user_cifs = StructureCode.objects.filter(user=user).count()
        # if user removed structures before
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
            if StructureCode.objects.filter(refcode=refcode).exists():
                refcode_obj.delete()
            return Response(
                {'errors': f'Structure information was not added! {error_message}'},
                status=status.HTTP_400_BAD_REQUEST
            )
        out_serializer = RefcodeFullSerializer(refcode_obj)
        return Response(
            out_serializer.data,
            status=status.HTTP_201_CREATED
        )


class QCStructureViewSet(StructureModelViewSet):
    filter_backends = (DjangoFilterBackend,)
    filterset_class = QCStructureFilter

    def get_queryset(self):
        if self.request.user.is_authenticated:
            queryset = (QCStructureCode.objects.filter(user=self.request.user) |
                        QCStructureCode.objects.filter(user__isnull=True))
            return queryset
        return QCStructureCode.objects.filter(user__isnull=True)

    def get_serializer_class(self):
        if self.action == 'list':
            return QCRefcodeShortSerializer
        return QCRefcodeFullSerializer

    @action(
        detail=True,
        methods=['GET'],
        url_path='export/cif'
    )
    def export_cif(self, request, pk):
        qc_structure = get_object_or_404(QCStructureCode, pk=pk)
        filename = f'{qc_structure.refcode}.cif'
        content = qc_get_cif_content(qc_structure)
        response = HttpResponse(content, content_type='text/plain')
        response['Content-Disposition'] = f'attachment; filename={filename}'
        return response

    @action(
        detail=True,
        methods=['GET'],
    )
    def download(self, request, pk):
        qc_structure = get_object_or_404(QCStructureCode, pk=pk)
        filename = f'{qc_structure.refcode}.txt'
        vasp_file = open(os.path.join(settings.BASE_DIR, 'media', str(qc_structure.vasp_file.file)), 'r')
        content = vasp_file.readlines()
        response = HttpResponse(content, content_type='text/xml')
        response['Content-Disposition'] = f'attachment; filename={filename}'
        return response

    @action(
        detail=True,
        methods=['GET'],
        url_path='export/2d'
    )
    def export_2d(self, request, pk):
        structure = get_object_or_404(QCStructureCode, pk=pk)
        return get_img2d(structure, request)

    @action(
        detail=False,
        methods=['POST'],
        permission_classes=[IsAuthenticated],
        serializer_class=[VaspUploadSerializer],
        url_path='upload/vasp'
    )
    def upload_vasp(self, request):
        user = request.user
        count_user_vasp = QCStructureCode.objects.filter(user=user).count()
        while True:
            refcode = 'user-' + str(user.id) + '-vasp-' + str(count_user_vasp + 1)
            if QCStructureCode.objects.filter(refcode=refcode).count() == 0:
                break
            count_user_vasp += 1
        file = request.FILES.get('file')
        refcode_obj = QCStructureCode.objects.create(user=user, refcode=refcode)
        vasp_file_obj = VaspFile.objects.create(refcode=refcode_obj, file=file)
        vasp_file_path = os.path.join(settings.BASE_DIR, 'media', str(vasp_file_obj.file))
        try:
            add_vasp_data(structure_obj=refcode_obj, file=vasp_file_path)
        except Exception as error_message:
            return Response(
                {'errors': f'Structure information was not added! {error_message}'},
                status=status.HTTP_400_BAD_REQUEST
            )
        out_serializer = QCRefcodeFullSerializer(refcode_obj)
        return Response(
            out_serializer.data,
            status=status.HTTP_201_CREATED
        )


def get_search_queryset_with_filtration(template):
    # TODO: переделать формат строки при подструктурной фильтрации!!! cord_min cord_max
    graphs = CoordinatesBlock.objects.filter(graph__isnull=False)
    filters = set_filter(template)
    filtrs = dict()
    for filtr, data in filters.items():
        is_true, obj_name = data
        if is_true:
            filtrs[f'refcode__{obj_name.lower()}__{filtr}'] = is_true
    structures = graphs.filter(**filtrs)
    analyse_data = structures.values_list('graph', flat=True)
    return analyse_data


def get_search_queryset_with_filtration_qc(template):
    graphs = QCCoordinatesBlock.objects.filter(graph__isnull=False)
    filters = set_filter(template)
    filtrs = dict()
    for filtr, data in filters.items():
        is_true, obj_name = data
        if is_true:
            filtrs[f'refcode__qc_{obj_name.lower()}__{filtr}'] = is_true
    structures = graphs.filter(**filtrs)
    analyse_data = structures.values_list('graph', flat=True)
    return analyse_data


def get_template_graph(graph, cord_num=True):
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
        if cord_num:
            result.append(graph.nodes[atom_idx]['cord_min'])
            result.append(graph.nodes[atom_idx]['cord_max'])
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
