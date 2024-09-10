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

import copy
from rdkit import Chem
from rdkit.Chem import Draw, rdDetermineBonds, AllChem, rdAbbreviations
from typing import List, Tuple, Dict
from collections import defaultdict
import multiprocessing
import time
import traceback

METAL_IONS = {
    'Li': 1,
    'Na': 1,
    'K': 1,
    'Rb': 1,
    'Cs': 1,
    'Mg': 2,
    'Ca': 2,
    'Sr': 2,
    'Ba': 2
}

ANIONS = {
    'F': -1,
    'Cl': -1,
    'Br': -1,
    'I': -1
}

MILT_ICHARGE_IONS = {
    'Bi': [3, 5]
}

# neutral
MAX_VALENCE = {
    'H': 1,
    'B': 3,
    'C': 4,
    'N': 3,
    'O': 2,
    'F': 1,
    'Al': 4,
    'Si': 4,
    'P': 5,
    'S': 6,
    'Cl': 1,
    'Ge': 4,
    'As': 5,
    'Se': 6,
    'Br': 1,
    'I': 1,
}

IONS = dict(**METAL_IONS, ** ANIONS)


def define_connect_from_graph(mol: Chem.Mol, bonds: List[Tuple]) -> Chem.Mol:
    mol = Chem.RWMol(mol)
    bnds: List[Chem.Bond] = list(mol.GetBonds())
    for bond in bnds:
        mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)
    for bond in bonds:
        mol.AddBond(bond[0], bond[1])
    mol = mol.GetMol()
    return mol


def run_process(mol_copy, mol_charge, return_dict):
    return_dict['error'] = ''
    try:
        rdDetermineBonds.DetermineBonds(mol_copy, charge=mol_charge, allowChargedFragments=True)
    except Exception as error_msg:
        return_dict['error'] = error_msg
    return_dict['mol'] = mol_copy
    return return_dict


def define_bonds_in_molecule_v2(
        rdkit_molecule: Chem.Mol, formula_init: List[str], structure_charge: int = 0, mols_num: int = 1, bonds: List[Tuple] = []
) -> Tuple[Chem.Mol, int, int]:
    charge_order_neg = [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6]
    charge_order_pos = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6]
    charge_order_neg.reverse()
    charge_order_pos.reverse()
    # draw_and_save_molecule(rdkit_molecule)
    init_flag: int = 0 - structure_charge
    while True:
        try:
            mol_copy = copy.deepcopy(rdkit_molecule)
            if rdkit_molecule.HasProp(key='charge') and int(rdkit_molecule.GetProp(key='charge')):
                mol_charge = int(rdkit_molecule.GetProp(key='charge'))
            elif init_flag:
                mol_charge = init_flag
                if init_flag > 0:
                    init_flag -= 1
                else:
                    init_flag += 1
            elif structure_charge >= 0:
                mol_charge = charge_order_neg.pop()
            else:
                mol_charge = charge_order_pos.pop()
            set_C_charge_zero(mol_copy)

            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            # process1 = threading.Thread(target=rdDetermineBonds.DetermineBonds, args=(mol_copy, ), kwargs={'charge': mol_charge, 'allowChargedFragments': True})
            process1 = multiprocessing.Process(target=run_process, args=(mol_copy, mol_charge, return_dict))
            process1.start()
            start = time.time()
            flag = False
            # kill if timeout (10.5 sec)
            while True:
                time.sleep(0.5)
                if not process1.is_alive():
                    flag = True
                    break
                # if the waiting time is exceeded, we interrupt the processes
                if time.time() - start > 10:
                    process1.terminate()
                    raise TimeoutError()

            mol_copy = return_dict.get('mol')
            # if critical error in rdDetermineBonds.DetermineBonds function
            if not flag:
                mol_copy = copy.deepcopy(rdkit_molecule)
                raise Exception('')
            if return_dict['error']:
                raise Exception(return_dict['error'])

            # if composition is not modified
            if formula_init == mol_to_formula(mol_copy):
                structure_charge += mol_charge * mols_num
                break
            elif rdkit_molecule.GetProp(key='charge'):
                raise Exception('ChargeError: Do not find available charge!')
        except Exception as err:
            if (not charge_order_neg) or (not charge_order_pos) or ('Final molecular charge does not match input' not in str(err)):
                if not charge_order_neg or not charge_order_pos:
                    print('Warring: Do not find available charge!')
                mol_copy = define_connect_from_graph(mol_copy, bonds)
                smiles = Chem.MolToSmiles(mol_copy, isomericSmiles=True, allHsExplicit=False)
                smol = Chem.MolFromSmiles(smiles, sanitize=False)
                draw_and_save_molecule(smol)
                return smol, structure_charge, None
    # define_connect_from_graph(mol_copy, bonds)
    return mol_copy, structure_charge, mol_charge


def mol_to_formula(mol: Chem.Mol) -> List[str]:
    formula = []
    comp = defaultdict(lambda: 0)
    for atom in mol.GetAtoms():
        element_symbol = atom.GetSymbol()
        comp[element_symbol] += 1
    sorted_element_keys = sorted(comp.keys())
    for element in sorted_element_keys:
        element += str(comp[element])
        formula.append(element)
    return formula


def set_atom_charge(mol: Chem.Mol):
    atoms: List[Chem.rdchem.Atom] = mol.GetAtoms()
    charge = 0
    for atom in atoms:
        if atom.GetSymbol() in IONS.keys():
            atom.SetFormalCharge(IONS[atom.GetSymbol()])
            charge += IONS[atom.GetSymbol()]
    mol.SetProp('charge', str(charge))
    return mol


def set_C_charge_zero(mol: Chem.Mol):
    atoms: List[Chem.rdchem.Atom] = mol.GetAtoms()
    for atom in atoms:
        if atom.GetSymbol() == 'C':
            atom.SetFormalCharge(0)
    return mol


def draw_and_save_molecule(mol: Chem.Mol, structure_name: str = 'noname'):
    abbrevs = rdAbbreviations.GetDefaultAbbreviations()
    for a in abbrevs:
        if Chem.MolToSmiles(a.mol) == '*C(=O)[OH]':
            a.label = 'COOH'
            a.displayLabel = 'COOH'
            a.displayLabelW = 'HOOC'
    mol_with_abbrev = rdAbbreviations.CondenseMolAbbreviations(mol, abbrevs)
    img = Draw.MolToImage(mol_with_abbrev, size=(1200, 1200))
    img.show()
    img.save(f'{structure_name}_2d.png')


def gen2d(smiles: str = '', inchis: list = [], sanitize: bool = True, size: Tuple[int, int] = (800, 800), format: str = 'img'):
    '''
    Supported formats:
        img - 2d image
        cml - ChemDraw format
    '''
    mol = None
    # first check smiles
    if smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    # if not smiles
    elif inchis:
        for inchi in inchis:
            if mol is None:
                mol = Chem.MolFromInchi(inchi, sanitize=sanitize)
            else:
                temp_mol = Chem.MolFromInchi(inchi)
                mol = Chem.CombineMols(mol, temp_mol)
    if mol is None:
        return 0
    # return formats
    if format == 'cml':
        return Chem.MolToCMLBlock(mol)
    else:
        return Draw.MolToImage(mol, size=size)


def get_symbol_from_element_number(el_number: int, element_numbers: Dict[str, int]) -> str:
    for key, value in element_numbers.items():
        if value == int(el_number):
            element = key
            return element


def construct_xyz_block(data, element_numbers: Dict[str, int], types: List[int], bonds: List[Tuple[int, int]]):
    xyz_block = str(len(data)) + '\n'
    xyz_block += 'cpplib_xyz\n'
    for atom in data:
        element_symbol = get_symbol_from_element_number(types[atom['init_idx']], element_numbers)
        xyz_block += f'{element_symbol} {round(atom["x"], 4)} {round(atom["y"], 4)} {round(atom["z"], 4)}\n'
    # bonds
    #new_bonds = []
    #for bond in bonds:
    #    atom1 = bond[0]
    #    atom2 = bond[1]
    #    for idx, atom in enumerate(data):
    #        if atom['init_idx'] == atom1:
    #            atom3 = idx
    #        elif atom['init_idx'] == atom2:
    #            atom4 = idx
    #    new_bonds.append((atom3, atom4))
    #return xyz_block, new_bonds
    return xyz_block


def reduce_charges(mol: Chem.Mol):
    atoms: List[Chem.rdchem.Atom] = mol.GetAtoms()
    for atom in atoms:
        if atom.GetFormalCharge():
            for neighbor in atom.GetNeighbors():
                a_chg = atom.GetFormalCharge()
                n_chg = neighbor.GetFormalCharge()
                # if charges are different
                if (a_chg > 0 > n_chg) or (a_chg < 0 < n_chg):
                    a_val, n_val = atom.GetTotalValence(), neighbor.GetTotalValence()
                    a_max_val = MAX_VALENCE.get(atom.GetSymbol(), 14)
                    n_max_val = MAX_VALENCE.get(neighbor.GetSymbol(), 14)
                    # if valences allowed
                    if a_val < a_max_val and n_val < n_max_val:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                            # change bond order
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                bond.SetBondType(Chem.BondType.DOUBLE)
                            elif bond.GetBondType() == Chem.BondType.DOUBLE:
                                bond.SetBondType(Chem.BondType.TRIPLE)
                            # change charges
                            if a_chg > 0:
                                atom.SetFormalCharge(a_chg - 1)
                                neighbor.SetFormalCharge(n_chg + 1)
                            else:
                                atom.SetFormalCharge(a_chg + 1)
                                neighbor.SetFormalCharge(n_chg - 1)
                            validate_a = atom.HasValenceViolation()
                            validate_n = neighbor.HasValenceViolation()
                            if not (validate_a and validate_n):
                                print(f'Warring: Bond order validation error ({atom.GetSymbol()}...{neighbor.GetSymbol()})')
    return mol


def reduce_radicals(mol: Chem.Mol):
    excluded_atoms = ['P', ]
    atoms: List[Chem.rdchem.Atom] = mol.GetAtoms()
    for atom in sorted(atoms, key=lambda x: x.GetAtomicNum(), reverse=True):
        if atom.GetNumRadicalElectrons() and atom.GetSymbol() not in excluded_atoms:
            neighbors = atom.GetNeighbors()
            for neighbor in sorted(neighbors, key=lambda x: x.GetAtomicNum(), reverse=True):
                a_rad = atom.GetNumRadicalElectrons()
                n_rad = neighbor.GetNumRadicalElectrons()
                if a_rad and n_rad and neighbor.GetSymbol() not in excluded_atoms:
                    a_val, n_val = atom.GetTotalValence(), neighbor.GetTotalValence()
                    a_max_val = MAX_VALENCE.get(atom.GetSymbol(), 14)
                    n_max_val = MAX_VALENCE.get(neighbor.GetSymbol(), 14)
                    # if valences allowed
                    if a_val < a_max_val and n_val < n_max_val:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                            # change bond order
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                bond.SetBondType(Chem.BondType.DOUBLE)
                            elif bond.GetBondType() == Chem.BondType.DOUBLE:
                                bond.SetBondType(Chem.BondType.TRIPLE)
                            # change radicals
                            atom.SetNumRadicalElectrons(a_rad - 1)
                            neighbor.SetNumRadicalElectrons(n_rad - 1)
                            validate_a = atom.HasValenceViolation()
                            validate_n = neighbor.HasValenceViolation()
                            if not (validate_a and validate_n):
                                print(f'Warring: Bond order validation error ({atom.GetSymbol()}...{neighbor.GetSymbol()})')
    return mol


def reduce_radicals_with_metalls(mol: Chem.Mol):
    metals = ['Li', 'Be', 'Na','Mg','Al','K', 'Ca','Sc', 'Ti', 'V','Cr', 'Mn', 'Fe',
               'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc',
               'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Cs', 'Ba', 'La',
               'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
               'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
               'Pb', 'Bi'
               ]
    atoms: List[Chem.rdchem.Atom] = mol.GetAtoms()
    for atom in atoms:
        if atom.GetNumRadicalElectrons():
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                a_rad = atom.GetNumRadicalElectrons()
                if a_rad and neighbor.GetSymbol() in metals:
                    a_val, n_val = atom.GetTotalValence(), neighbor.GetTotalValence()
                    a_max_val = MAX_VALENCE.get(atom.GetSymbol(), 14)
                    n_max_val = MAX_VALENCE.get(neighbor.GetSymbol(), 14)
                    # if valences allowed
                    if a_val < a_max_val and n_val < n_max_val:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() is Chem.BondType.SINGLE:
                            bond.SetBondType(Chem.BondType.IONIC)
                            # change radicals
                            atom.SetNumRadicalElectrons(a_rad - 1)
    return mol


def validate_mol(mol: Chem.Mol):
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons():
            count += 1
    if count:
        return False
    return True


def main_v2(xyz_mols, element_numbers: Dict[str, int], types: List[int]):
    structure_charge = 0
    xyz_blocks = []
    # make xyz format block
    for xyz_mol in xyz_mols:
        if len(xyz_mol['atoms']):
            xyz_block = construct_xyz_block(xyz_mol['atoms'], element_numbers, types, xyz_mol['bonds'])
            # for test
            xyz_blocks.append([xyz_mol['count'], xyz_block, xyz_mol['bonds']])
    xyz_blocks = sorted(xyz_blocks, key=lambda i: int(i[1][0].split('\n')[0]))
    final_structure = None
    for xyz_block in xyz_blocks:
        mols_num, xyz_block, bonds = xyz_block
        # make rdkit mol
        raw_mol = Chem.MolFromXYZBlock(xyz_block)
        rd_mol: Chem.Mol = Chem.Mol(raw_mol)
        # TODO: добавить заряды на атомах, если были указаны в исходном сифе
        if rd_mol.GetNumAtoms() == 1:
            set_atom_charge(rd_mol)
        formula = mol_to_formula(rd_mol)
        # define bonds and bond orders in each molecule
        rd_mol, structure_charge, mol_charge = define_bonds_in_molecule_v2(rd_mol, formula, structure_charge, mols_num, bonds)
        # merge molecules
        if rd_mol is not None:
            if final_structure is None:
                final_structure = rd_mol
            else:
                final_structure = Chem.CombineMols(final_structure, rd_mol)
    # final generation
    smiles = Chem.MolToSmiles(final_structure, isomericSmiles=True)
    smol = Chem.MolFromSmiles(smiles, sanitize=True)
    set_C_charge_zero(smol)
    Chem.SanitizeMol(smol)
    reduce_radicals(smol)
    reduce_charges(smol)
    reduce_radicals_with_metalls(smol)
    if validate_mol(smol):
        try:
            AllChem.GenerateDepictionMatching3DStructure(smol, reference=final_structure)
        except:
            AllChem.Compute2DCoords(smol)
        # draw the picture
        # to set draw options use rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions (only for Draw.MolToImage)
        # to draw aromatic rings use kekulize=False (only for Draw.MolToImage)
        result_data = {
            'mol_obj_with_coords': final_structure,
            'mol_obj_pretty': smol, 'smiles': Chem.MolToSmiles(smol, isomericSmiles=True),
            'inchi': Chem.MolToInchi(smol),
            'formula': ''.join(mol_to_formula(final_structure))
        }
        return result_data
    return 0
