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
import cpplib
from ._element_numbers import element_numbers
from django_project.loggers import set_prm_log
from modules.gen2d.gen2d import main_v2
import re


def get_coords(cif_block) -> str:
    atomic_sites = ''
    atoms = []
    coords = cif_block[1].GetLoop('_atom_site_label')
    order: list = coords.GetItemOrder()
    lable_idx = order.index('_atom_site_label')
    if '_atom_site_type_symbol' in order:
        atom_type_idx = order.index('_atom_site_type_symbol')
    else:
        atom_type_idx = ''
    x_idx = order.index('_atom_site_fract_x')
    y_idx = order.index('_atom_site_fract_y')
    z_idx = order.index('_atom_site_fract_z')
    for site in coords:
        atoms.append(site[lable_idx])
        temp = list()
        temp.append(site[lable_idx])
        if atom_type_idx:
            temp.append(site[atom_type_idx])
        else:
            value = re.findall(r'(^[a-zA-Z]{1,3})', site[lable_idx])
            if value:
                value = value[0]
                temp.append(value)
            else:
                raise Exception('No "_atom_site_type_symbol" key was found in cif file!')
        temp.append(site[x_idx])
        temp.append(site[y_idx])
        temp.append(site[z_idx])
        for j, element in enumerate(temp, start=0):
            if j == 0:
                temp[j] = element.replace('?', '')
            if '(' in element:
                idx = element.index('(')
                temp[j] = element[:idx]
        atomic_sites += ' '.join(temp)
        atomic_sites += ' '
    return atomic_sites, atoms


def get_data(cif_block, symops_db):
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
            element = re.findall(r'[A-Za-z]{1,3}', element)
            if element and element[0] in element_numbers.keys():
                atoms_types.append(element_numbers[element[0]])
            else:
                raise Exception(f'Error: check the element: {element}')
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
    else:
        symops_list = symops_db.split(';')
    return params, atoms_coords, atoms_types, symops_list


def make_graph_c(params, coords, types, refcode, add_graphs_logger, symops):
    cpplib_result = cpplib.FindMoleculesInCell(params, symops, types, coords)
    mols_str_graph = cpplib_result['graph_str']
    xyz_mols = cpplib_result['xyz_block']
    parse_mols = mols_str_graph.split(';')
    if len(parse_mols) == 1:
        graph_str, = parse_mols
    elif len(parse_mols) == 2:
        graph_str, warnings = parse_mols
        add_graphs_logger.warning(f"FindMoleculesInCellError in {refcode}:\n\t{warnings}")
    else:
        raise Exception(f"Invalid length of mols list: {len(parse_mols)}")
    if graph_str.split()[1] == 0:
        raise Exception(f"There are no atoms in graph! May be the structure was unordered")
    # generate data for 2d graph picture
    data_2d = main_v2(xyz_mols, element_numbers, types)
    smiles = ''
    inchi = ''
    if data_2d:
        smiles = data_2d['smiles']
        inchi = data_2d['inchi']
    return graph_str, smiles, inchi


def save(params, coords, types, refcode, symops):
    f = open('test_' + str(refcode) + '.txt', 'w')
    f.write(str(params).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(symops).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(len(symops)).replace('[', '').replace(']', '') + ', ')
    f.write(str(types).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(coords).replace('[', '{').replace(']', '}') + ', ')
    f.write(str(len(types)).replace('[', '').replace(']', ''))
    f.close()
    exit()


def print_graph(graphs):
    import matplotlib.pyplot as plt
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
            refcode, cif_block, symops_db = queue.get(block=True, timeout=0.5)
        except Exception:
            # if the queue is empty, terminate the thread
            break
        if True:
        #try:
            add_graphs_logger.info(f"Start processing structure {refcode}")
            params, coords, types, symops = get_data(cif_block, symops_db)
            add_graphs_logger.info(f"Received atomic coordinates and translation matrix")
            graph_str, smiles, inchi = make_graph_c(params, coords, types, refcode, add_graphs_logger, symops)
            if smiles and inchi:
                add_graphs_logger.info(f"Received graph string and 2D representation")
            else:
                add_graphs_logger.info(f"Build 2D representation failed!")
            add_graphs_logger.info(f"Processing completed {refcode}")
            bonds = []
            angles = []
            return_dict[refcode] = {'graph_str': graph_str, 'bonds': bonds, 'angles': angles, 'smiles': smiles, 'inchi': inchi}
            add_graphs_logger.info(f"Structure {refcode} successfully added to the resulting list!")
        #except Exception as err:
        #    add_graphs_logger.error(f"Structure {refcode} not added to the resulting list!", exc_info=True)
    add_graphs_logger.info(f"The queue is empty, the thread completed successfully!")
    return 0
