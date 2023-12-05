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

import networkx as nx
import matplotlib.pyplot as plt
from django.conf import settings
import ctypes
import json
from ._element_numbers import element_numbers
import os
from django_project.loggers import set_prm_log


def get_coords(cif_block) -> str:
    atomic_sites = ''
    atoms = []
    coords = cif_block[1].GetLoop('_atom_site_label')
    for i, site in enumerate(coords, start=0):
        atoms.append(site[0])
        temp = site.copy()
        for j, element in enumerate(temp, start=0):
            if j == 0:
                temp[j] = element.replace('?', '')
            if '(' in element:
                idx = element.index('(')
                temp[j] = element[:idx]
            if j > 4:
                break
        atomic_sites += ' '.join(temp)
        atomic_sites += ' '
    return atomic_sites, atoms


def get_data(cif_block):
    params = [
        cif_block[1]['_cell_length_a'].split('(')[0], cif_block[1]['_cell_length_b'].split('(')[0],
        cif_block[1]['_cell_length_c'].split('(')[0], cif_block[1]['_cell_angle_alpha'].split('(')[0],
        cif_block[1]['_cell_angle_beta'].split('(')[0], cif_block[1]['_cell_angle_gamma'].split('(')[0]
    ]
    params = list(map(float, params))
    atomic_sites, atoms = get_coords(cif_block)
    atoms_coords = []
    atoms_types = []
    atoms_info = atomic_sites.split()
    for idx, element in enumerate(atoms_info, start=1):
        if idx % 5 == 2:
            atoms_types.append(element_numbers[element])
        elif idx % 5 == 3:
            coords = [atoms_info[idx - 1], atoms_info[idx], atoms_info[idx + 1]]
            atoms_coords.extend(map(float, coords))
    # symops
    symops_list = []
    if '_symmetry_equiv_pos_as_xyz' in cif_block[1].keys():
        symops = cif_block[1].GetLoop('_symmetry_equiv_pos_as_xyz')
        symops_order = symops.GetItemOrder()
        symops_idx = symops_order.index('_symmetry_equiv_pos_as_xyz')
        for symop in symops:
            if symop:
                symops_list.append(symop[symops_idx])
    elif '_symmetry_space_group_name_H-M' in cif_block[1].keys():
        h_m_group = cif_block[1]['_symmetry_space_group_name_H-M']
        file = open(os.path.join(os.path.dirname(__file__), '../symops.json'))
        symops_json = json.load(file)
        for i, item in enumerate(symops_json):
            if i > 0:
                if item['universal_h_m'] == h_m_group:
                    symops_list = item['symops']
                    break
    if not symops_list:
        raise Exception('No symmetry operations were found!')
    return params, atoms_coords, atoms_types, symops_list


def make_graph_c(params, coords, types, refcode, add_graphs_logger, symops):
    dll = settings.GET_DLL()
    c_params = (ctypes.c_float * len(params))(*params)
    c_types = (ctypes.c_int * len(types))(*types)
    c_coords = (ctypes.c_float * len(coords))(*coords)
    c_symops = (ctypes.c_char_p * len(symops))(*[s.encode() for s in symops])
    mols = dll.FindMoleculesInCell(c_params, c_symops, len(symops), c_types, c_coords, len(types))
    # transform to class str
    mols = mols.decode()
    parse_mols = mols.split(';')
    if len(parse_mols) == 1:
        graph_str, = parse_mols
    elif len(parse_mols) == 2:
        graph_str, warnings = parse_mols
        add_graphs_logger.warning(f"FindMoleculesInCellError in {refcode}:\n\t{warnings}")
    else:
        raise Exception(f"Invalid lenth of mols list: {len(parse_mols)}")
    if graph_str.split()[1] == 0:
        raise Exception(f"There are no atoms in graph! May be the structure was unordered")
    return graph_str


def save(params, coords, types, refcode, symops):
    f = open('test_' + refcode + '.txt', 'w')
    f.write(str(params).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(symops).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(len(symops)).replace('[', '').replace(']', '') + ', ')
    f.write(str(types).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(coords).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(len(types)).replace('[', '').replace(']', ''))
    f.close()
    exit()


def print_graph(graphs):
    for item in set(graphs):
        print(item)
        elements = nx.get_node_attributes(item, "element")
        nx.draw(item, with_labels=True, labels=elements, pos=nx.spring_layout(item))
        plt.gca().set_aspect('equal')
        plt.show()


def print_graph_c(string_graph):
    graph = nx.Graph()
    string_graph = string_graph.split()
    num_nodes = int(string_graph[1])
    for idx, el in enumerate(string_graph[3:], start=1):
        # nodes
        if idx < num_nodes * 2 and idx % 2 == 1:
            for key, value in element_numbers.items():
                if value == int(el):
                    element = key
                    break
            graph.add_nodes_from([(int((idx + 1) / 2), {'element': element}), ])
        # edges
        elif idx > num_nodes * 2 and idx % 2 == 1:
            graph.add_edge(int(el), int(string_graph[idx + 3]))
    print_graph([graph, ])


def add_graphs_c(queue, return_dict, proc_num: int):
    # Set up logger
    add_graphs_logger = set_prm_log(proc_num)
    while True:
        try:
            refcode, cif_block = queue.get(block=True, timeout=0.5)
        except Exception:
            # if the queue is empty, terminate the thread
            break
        try:
            add_graphs_logger.info(f"Start processing structure {refcode}")
            params, coords, types, symops = get_data(cif_block)
            add_graphs_logger.info(f"Received atomic coordinates and translation matrix")
            graph_str = make_graph_c(params, coords, types, refcode, add_graphs_logger, symops)
            add_graphs_logger.info(f"Received graph string")
            bonds = []
            add_graphs_logger.info(f"Processing completed {refcode}")
            angles = []
            return_dict[refcode] = [graph_str, bonds, angles]
            add_graphs_logger.info(f"Structure {refcode} successfully added to the resulting list!")
        except Exception as err:
            add_graphs_logger.error(f"Structure {refcode} not added to the resulting list!", exc_info=True)
    add_graphs_logger.info(f"The queue is empty, the thread completed successfully!")
    return 0
