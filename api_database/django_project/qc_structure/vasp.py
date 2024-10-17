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

import os
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from structure.models import Spacegroup, SYSTEMS, CENTRINGS, get_elements_list
import json
from django.conf import settings
from api.filters import get_reduced_cell
from django_project.loggers import vasp_logger
from math import cos, sqrt, radians
import networkx as nx
import re
from structure.management.commands.cif_db_update_modules._element_numbers import element_numbers
from structure.management.commands.cif_db_update_modules._make_graphs_c import make_graph_c
from structure.management.commands.cif_db_update_modules._add_substructure_filtration import (TEMPLATES, start_dll_and_write,
                                                                                              set_only_CHNO, set_no_C,
                                                                                              set_elements, SET_ELEMENTS)
from qc_structure.models import (QCStructureCode, QCCell, QCReducedCell, QCFormula,
                                 QCCompoundName, QCElementsManager, QCProperties,
                                 QCCoordinatesBlock, QCSubstructure1, QCSubstructure2,
                                 QCProgram, QCInChI)


def save_program(struct_obj):
    vasp_logger.info('Add program name...')
    prog, created = QCProgram.objects.get_or_create(refcode=struct_obj, vasp=True)


def save_name(struct_obj, syst_name, triv_name):
    vasp_logger.info('Add compound name...')
    name, created = QCCompoundName.objects.get_or_create(refcode=struct_obj)
    if syst_name:
        name.systematic_name = syst_name
    if triv_name:
        name.trivial_name = triv_name
    name.save()


def save_properties(struct_obj, vasp_out):
    vasp_logger.info('Add properties info...')
    energy, units = str(vasp_out.final_energy).split()
    density = vasp_out.final_structure.density
    prop, created = QCProperties.objects.get_or_create(refcode=struct_obj)
    prop.energy = energy
    prop.calculated_density = round(density, 3)
    prop.save()


def get_or_create_space_group(vasp_structure, return_only_symops=False):

    def get_system_id(system):
        system_id = 0
        for syst in SYSTEMS:
            if system in syst:
                system_id = syst[0]
        if not system_id:
            vasp_logger.warning(f'Unknown lattice system: {system}!')
        return system_id

    def get_symops(vasp_hall):
        file = open(os.path.join(settings.BASE_DIR, 'structure', 'management', 'commands', 'symops.json'))
        symops_list = json.load(file)
        for i, symops in enumerate(symops_list):
            if i > 0:
                # lstrip() to delete beginning spaces
                json_hall = symops['hall'].lstrip()
                if json_hall == vasp_hall:
                    operations = ";".join(symops['symops'])
                    return operations
        vasp_logger.error(f'No symmetry was found with hall name {vasp_hall}')
        raise Exception(f'No symmetry was found with hall name {vasp_hall}')

    vasp_logger.info('Determine space group...')
    spgran = SpacegroupAnalyzer(vasp_structure)
    symmetry_info = spgran.get_symmetry_dataset()
    hall = symmetry_info['hall']
    sg = Spacegroup.objects.filter(hall_name=hall)
    if sg.exists():
        if sg.count() == 1:
            space_group = Spacegroup.objects.get(hall_name=hall)
        elif sg.filter(symops__isnull=False).exists():
            space_group = sg.filter(symops__isnull=False).first()
        else:
            space_group = sg.first()
    else:
        space_group_name = symmetry_info['international'].replace('_', '')
        system = spgran.get_crystal_system()
        system_id = get_system_id(system)
        number = int(symmetry_info['number'])
        symops = get_symops(hall)
        if return_only_symops:
            return symops
        space_group, create = Spacegroup.objects.get_or_create(
            name=space_group_name,
            number=int(number),
            system=system_id,
            hall_name=hall,
            symops=symops
        )
    if return_only_symops:
        return space_group.symops
    return space_group


def save_cell(struct_obj, vasp_structure, space_group):

    def get_centring_id(centring):
        for item in CENTRINGS:
            if centring.replace('-', '') in item:
                return item[0]
        raise Exception('Invalid unit cell centring!')

    vasp_logger.info('Add unit cell...')
    spgran = SpacegroupAnalyzer(vasp_structure, symprec=0.02)
    symmed_vasp_struct = spgran.get_refined_structure()
    cif_form_vasp = symmed_vasp_struct.to(fmt='cif', symprec=0.02).split('\n')
    for line in cif_form_vasp:
        line = line.split()
        if line and line[0] == '_cell_formula_units_Z':
            z_val = int(line[1])
    a, b, c = symmed_vasp_struct.lattice.abc
    al, be, ga = symmed_vasp_struct.lattice.angles
    centring = space_group.name[0]
    centring_id = get_centring_id(centring)
    cell, created = QCCell.objects.get_or_create(
        a=float(a), b=float(b), c=float(c),
        al=float(al), be=float(be), ga=float(ga),
        spacegroup=space_group, refcode=struct_obj,
        centring=centring_id, zvalue=z_val
    )
    return symmed_vasp_struct


def save_reduced_cell(struct_obj):
    vasp_logger.info('Add reduced cell...')
    centrings = dict((v, k) for v, k in CENTRINGS)
    centring = centrings[struct_obj.qc_cell.centring]
    params = [struct_obj.qc_cell.a, struct_obj.qc_cell.b, struct_obj.qc_cell.c,
              struct_obj.qc_cell.al, struct_obj.qc_cell.be, struct_obj.qc_cell.ga]
    reduced_params = get_reduced_cell(params, centring)
    a, b, c, al, be, ga = reduced_params
    volume = (
            a * b * c * sqrt(1 + 2 * cos(radians(al)) * cos(radians(be)) *
                             cos(radians(ga)) - cos(radians(al)) ** 2 -
                             cos(radians(be)) ** 2 - cos(radians(ga)) ** 2)
    )
    rc, created = QCReducedCell.objects.get_or_create(
        refcode=struct_obj,
        a=round(a, 3), b=round(b, 3), c=round(c, 3),
        al=round(al, 3), be=round(be, 3), ga=round(ga, 3),
        volume=round(volume, 3),
    )


def save_coordinates(struct_obj, symmed_vasp_struct, return_only_str_sites=False):
    vasp_logger.info('Add coordinates...')
    cif_form_vasp = symmed_vasp_struct.to(fmt='cif', symprec=0.02)
    from pymatgen.io.cif import CifParser
    cif = CifParser.from_str(cif_form_vasp)
    cif_info = list(cif.as_dict().values())[0]
    atom_types = cif_info['_atom_site_type_symbol']
    x_coord = cif_info['_atom_site_fract_x']
    y_coord = cif_info['_atom_site_fract_y']
    z_coord = cif_info['_atom_site_fract_z']
    str_sites = ''
    for idx in range(len(atom_types)):
        str_sites += f'{atom_types[idx]}{idx + 1} {atom_types[idx]} {x_coord[idx]} {y_coord[idx]} {z_coord[idx]}\n'
    if return_only_str_sites:
        return str_sites
    cb_obj, created = QCCoordinatesBlock.objects.get_or_create(refcode=struct_obj)
    cb_obj.coordinates = str_sites
    cb_obj.save()


def save_graph(struct_obj):
    vasp_logger.info('Create and save graph...')
    params = [struct_obj.qc_cell.a, struct_obj.qc_cell.b, struct_obj.qc_cell.c,
              struct_obj.qc_cell.al, struct_obj.qc_cell.be, struct_obj.qc_cell.ga]
    sites_info = struct_obj.qc_coordinates.coordinates.split('\n')
    atoms_types = []
    atoms_coords_types = []
    for site in sites_info:
        if site:
            site_info = site.split()
            element = site_info[1]
            atoms_types.append(element_numbers[element])
            coords = [element_numbers[element], float(site_info[2]), float(site_info[3]), float(site_info[4])]
            atoms_coords_types.append(tuple(coords))
    symops = struct_obj.qc_cell.spacegroup.symops.split(';')
    # create graph
    graph_str, smiles, inchi = make_graph_c(
        params, atoms_coords_types, atoms_types,
        struct_obj, vasp_logger, symops
    )
    graph_str = str(struct_obj.id) + ' ' + graph_str
    # save graph
    if graph_str:
        struct_obj.qc_coordinates.graph = graph_str
        struct_obj.qc_coordinates.save()
    return smiles, inchi


def save_smiles_inchi(structure_obj, qc_smiles, qc_inchi):
    vasp_logger.info('Save smiles and inchi...')
    coord_block = QCCoordinatesBlock.objects.get(refcode=structure_obj)
    if qc_smiles:
        coord_block.smiles = qc_smiles
        coord_block.save()
    if qc_inchi:
        inchi = qc_inchi.split('=')[1].split('/')
        inchi_block = QCInChI.objects.create(refcode=structure_obj, version=inchi[0], formula=inchi[1])
        for item in inchi[2:]:
            if item.startswith('c'):
                inchi_block.connectivity = item
            elif item.startswith('h'):
                inchi_block.hydrogens = item
            elif item.startswith('q'):
                inchi_block.q_charge = item
            elif item.startswith('p'):
                inchi_block.p_charge = item
            elif item.startswith('b'):
                inchi_block.b_stereo = item
            elif item.startswith('t'):
                inchi_block.t_stereo = item
            elif item.startswith('m'):
                inchi_block.m_stereo = item
            elif item.startswith('s'):
                inchi_block.s_stereo = item
            elif item.startswith('i'):
                inchi_block.i_isotopic = item
        inchi_block.save()


def save_formula(structure_obj, symmed_vasp_struct):

    def make_networkx_graph(graph_string):
        graph = nx.Graph()
        string_graph = graph_string.split()
        num_nodes = int(string_graph[1])
        for idx, el in enumerate(string_graph[3:], start=1):
            # nodes
            if idx < num_nodes * 2 and idx % 2 == 1:
                for key, value in element_numbers.items():
                    if value == int(el):
                        element = key
                        break
                graph.add_nodes_from([(int((idx + 1) / 2), {'element': element, 'h_num': int(string_graph[idx + 3])}), ])
            # edges
            elif idx > num_nodes * 2 and idx % 2 == 1:
                graph.add_edge(int(el), int(string_graph[idx + 3]))
        return graph

    vasp_logger.info('Add formula...')
    sum_formula = symmed_vasp_struct.composition

    # moiety_formula
    graph_string = structure_obj.qc_coordinates.graph
    graph = make_networkx_graph(graph_string)
    molecules = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
    moiety_formula = list()
    for mol in molecules:
        atoms = {'H': 0}
        atom_types = list(nx.get_node_attributes(mol, 'element').values())
        h_nums = list(nx.get_node_attributes(mol, 'h_num').values())
        for idx in range(len(atom_types)):
            atom_type = atom_types[idx]
            h_num = h_nums[idx]
            atoms['H'] += h_num
            if atom_type not in atoms.keys():
                atoms[atom_type] = 0
            atoms[atom_type] += 1
        temp = list()
        for atom, count in atoms.items():
            temp.append(f'{atom}{count}')
        moiety_formula.append(' '.join(temp))
    moiety_formula = ','.join(moiety_formula)
    # save formula
    formula_obj, created = QCFormula.objects.get_or_create(refcode=structure_obj)
    formula_obj.formula_moiety = moiety_formula
    formula_obj.formula_sum = str(sum_formula)
    formula_obj.save()


def save_element_sets(structure_obj):
    vasp_logger.info('Add elements...')
    formula = structure_obj.qc_formula.formula_sum
    elements = formula.split()
    elements_from_formula = dict()
    for element in elements:
        atom_type = re.findall(r'[A-Za-z]{1,3}', element)[0]
        count = re.findall(r'\d+', element)
        if count:
            count = float(count[0])
        else:
            count = 1
        elements_from_formula[atom_type] = count
    el_manager, created = QCElementsManager.objects.get_or_create(refcode=structure_obj)
    el_manager.save_elements(elements_from_formula)


def save_substructure(structure_obj):
    vasp_logger.info('Add substructure filtration...')
    models = {'Substructure1': QCSubstructure1,
              'Substructure2': QCSubstructure2}
    graph = structure_obj.qc_coordinates.graph
    if graph:
        for attr_name, data in TEMPLATES.items():
            template_graph, obj_name = data
            start_dll_and_write(template_graph, [graph, ], 1, attr_name, models[obj_name], QCStructureCode)

    graph_query = QCCoordinatesBlock.objects.filter(refcode=structure_obj)
    set_only_CHNO(graph_query, QCSubstructure1, 'refcode__qc_elements__element_set')
    set_no_C(graph_query, QCSubstructure1, 'refcode__qc_elements__element_set')
    for attr_name, element_set in SET_ELEMENTS.items():
        set_elements(graph_query, attr_name, element_set, QCSubstructure1, 'refcode__qc_elements__element_set')


def vasp_parser(structure_obj, file: str, syst_name='', triv_name=''):
    """Read and parse vasprun.xml output file"""
    save_program(structure_obj)
    save_name(structure_obj, syst_name, triv_name)
    vasp_out = Vasprun(file)
    vasp_structure = vasp_out.final_structure
    save_properties(structure_obj, vasp_out)
    space_group = get_or_create_space_group(vasp_structure)
    symmed_vasp_struct = save_cell(structure_obj, vasp_structure, space_group)
    save_reduced_cell(structure_obj)
    save_coordinates(structure_obj, symmed_vasp_struct)
    smiles, inchi = save_graph(structure_obj)
    save_smiles_inchi(structure_obj, smiles, inchi)
    save_formula(structure_obj, symmed_vasp_struct)
    save_element_sets(structure_obj)
    save_substructure(structure_obj)
