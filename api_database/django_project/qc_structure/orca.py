import time

from django_project.loggers import orca_logger
from qc_structure.vasp import save_program, save_name, save_smiles_inchi, save_substructure
from orca_parser import ORCAParse
import re
from structure.management.commands.cif_db_update_modules._element_numbers import element_numbers
from structure.management.commands.cif_db_update_modules._make_graphs_c import make_graph_c
from qc_structure.models import (QCStructureCode, QCCell, QCReducedCell, QCFormula,
                                 QCCompoundName, QCElementsManager, QCProperties,
                                 QCCoordinatesBlock, QCSubstructure1, QCSubstructure2,
                                 QCProgram, QCInChI, QCEnergy, QCInputParameters)


def save_energies(struct_obj, orca_out):
    orca_logger.info('Add energies info...')
    orca_out.parse_energies()
    energy = orca_out.energies
    energy = round(energy[-1], 6)
    energ, created = QCEnergy.objects.get_or_create(refcode=struct_obj)
    energ.energy = energy
    orca_out.parse_free_energy()
    if orca_out.enthalpies:
        energ.enthalpy = round(list(orca_out.enthalpies.values())[-1], 6)
    if orca_out.entropies:
        energ.entropy = round(list(orca_out.entropies.values())[-1], 6)
    if orca_out.AllGibbs:
        energ.gibbs = round(orca_out.Gibbs, 6)
    if "Zero point energy" in orca_out.raw:
        zpe = orca_out.raw.split("Zero point energy")[1].split("Eh")[0].split()[1]
        energ.zpe = round(float(zpe), 6)
    orca_out.parse_HOMO_LUMO()
    energ.homo = orca_out.all_HOMO[-1][-2]
    energ.lumo = orca_out.all_LUMO[-1][-2]
    energ.save()


def save_properties(struct_obj, orca_out):
    orca_logger.info('Add properties info...')
    input_data = orca_out.parse_input()
    prop, created = QCProperties.objects.get_or_create(refcode=struct_obj)
    prop.charge = orca_out.Z
    prop.multiplicity = orca_out.Multiplicity
    dipole = orca_out.parse_dipole()
    if type(dipole) is dict:
        prop.dipole_moment = round(dipole['Magnitude (Debye)'][0], 3)
    prop.save()
    return input_data


def save_input(struct_obj, data):
    orca_logger.info('Add input parameters info...')
    inp, created = QCInputParameters.objects.get_or_create(refcode=struct_obj)
    for field, key in {
        'job': 'Job', 'basis_set': 'BasisSet', 'functional': 'Functional',
        'version': 'version', 'freq': 'Freq', 'dispersion': 'Dispersion',
        'solvation': 'Solvation'
    }.items():
        if key in data.keys() and data[key]:
            setattr(inp, field, data[key])
            print(field)
    inp.save()


def save_coordinates(struct_obj, orca_out, return_only_str_sites=False):
    orca_logger.info('Add coordinates...')
    orca_out.parse_coords()
    atom_types = orca_out.atoms
    coords_final = orca_out.coords[-1]
    str_sites = ''
    for idx in range(len(atom_types)):
        coords = coords_final[idx]
        str_sites += f'{atom_types[idx]}{idx + 1} {coords[0]} {coords[1]} {coords[2]}\n'
    if return_only_str_sites:
        return str_sites
    cb_obj, created = QCCoordinatesBlock.objects.get_or_create(refcode=struct_obj)
    cb_obj.coordinates = str_sites
    cb_obj.is_fractional = False
    cb_obj.save()


def save_formula_and_elements(structure_obj, orca_out):
    atoms = orca_out.atoms
    elements_count = {}
    formula = ''
    for atom in atoms:
        if atom not in elements_count.keys():
            elements_count[atom] = 1
        else:
            elements_count[atom] += 1
    elements_count = dict(sorted(elements_count.items()))
    for element, count in elements_count.items():
        formula += element + str(count)
    formula_obj, created = QCFormula.objects.get_or_create(refcode=structure_obj)
    formula_obj.formula_sum = formula
    formula_obj.save()
    # save elements
    orca_logger.info('Add elements...')
    el_manager, created = QCElementsManager.objects.get_or_create(refcode=structure_obj)
    el_manager.save_elements(elements_count)


def save_graph(struct_obj):
    orca_logger.info('Create and save graph...')
    sites_info = struct_obj.qc_coordinates.coordinates.split('\n')
    atoms_types = []
    atoms_coords_types = []
    for site in sites_info:
        if site:
            site_info = site.split()
            element = re.findall(r'[A-Za-z]{1,3}', site_info[0])[0]
            atoms_types.append(element_numbers[element])
            coords = [element_numbers[element], float(site_info[1]), float(site_info[2]), float(site_info[3])]
            atoms_coords_types.append(tuple(coords))
    # create graph
    graph_str, smiles, inchi = make_graph_c(
        [], atoms_coords_types, atoms_types,
        struct_obj, orca_logger, []
    )
    graph_str = str(struct_obj.id) + ' ' + graph_str
    # save graph
    if graph_str:
        struct_obj.qc_coordinates.graph = graph_str
        struct_obj.qc_coordinates.save()
    return smiles, inchi


def orca_parser(structure_obj, file: str, syst_name='', triv_name=''):
    """Read and parse orca .out output file"""
    orca_logger.info('Add program name...')
    save_program(structure_obj, 'orca')
    orca_logger.info('Add compound name...')
    save_name(structure_obj, syst_name, triv_name)
    # parse orca file
    try:
        orca_out = ORCAParse(file)
    except:
        raise Exception('Error in parsing file')
    if not orca_out.valid:
        raise Exception('Input file is not valid!')
    save_energies(structure_obj, orca_out)
    input_data = save_properties(structure_obj, orca_out)
    save_input(structure_obj, input_data)
    save_coordinates(structure_obj, orca_out)
    smiles, inchi = save_graph(structure_obj)
    orca_logger.info('Save smiles and inchi...')
    save_smiles_inchi(structure_obj, smiles, inchi)
    save_formula_and_elements(structure_obj, orca_out)
    orca_logger.info('Add substructure filtration...')
    save_substructure(structure_obj)
