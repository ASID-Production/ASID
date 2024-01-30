from ... import point_class
import numpy as np
from ..ChemPack import PALETTE, MOLECULE_SYSTEMS
import os
from . import MoleculeClass


class FileParser:

    def __init__(self):
        self.SUPPORTED_FORMATS = {'.xyz': self.parsXyz,
                                  '.pdb': self.parsPdb}

    def parsFile(self, file_path, bond=True, root=None):
        _, ext = os.path.splitext(file_path)
        if ext not in self.SUPPORTED_FORMATS:
            print('File is not supported')
            return None, None
        return self.SUPPORTED_FORMATS[ext](file_path, bond, root)

    def parsXyz(self, file_path, bond, root):
        mol_list, atom_list, bonds_l = None, None, None
        file = open(file_path, 'r')
        mol_sys = MoleculeClass.MoleculeSystem()
        mol = MoleculeClass.Molecule(parent=mol_sys)
        for line in file:
            line_p = [x for x in line[:-1].split(' ') if x]
            if PALETTE.getName(line_p[0]) != 999 and all([True if x.replace('.', '').replace('-', '').isdigit() else False for x in line_p[1:4]]):
                coord = np.array([float(x) for x in line_p[1:]], dtype=np.float32)
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(line_p[0]), parent=mol)
        if bond:
            mol_sys.genBonds()
        if root is not None:
            mol_list = point_class.PointsList(parent=root, name=os.path.basename(file_path).split('.')[0])
            mol.assignPoint(mol_list)
            MOLECULE_SYSTEMS[mol_list] = mol_sys
            atom_list = point_class.PointsList(parent=mol_list, rad=0.25, name='Atoms')
            i = 1
            for atom in mol.children:
                coord = atom.coord.copy()
                point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                          color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                          atom_type=atom.atom_type)
                atom.assignPoint(point)
                point.addProperty('name', f'{point.atom_type}{i}')
                point.addProperty('label', f'{point.atom_type}{i}')
                i += 1
            if bond:
                bonds_l = point_class.PointsList(parent=mol_list, rad=0.1, name='Bonds')
                bonds = []
                for atom in mol:
                    for bond in atom.bonds():
                        if bond not in bonds:
                            bond_l = point_class.PointsList(parent=bonds_l, name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}', rad=bonds_l)
                            b1 = point_class.Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            b2 = point_class.Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            bonds.append(bond)
        return mol_sys, (mol_list, atom_list, bonds_l)

    def parsPdb(self, file_path, bond, root):
        mol_list, atom_list, bonds_l = None, None, None
        file = open(file_path, 'r')
        mol_sys = MoleculeClass.MoleculeSystem()
        mol = MoleculeClass.Molecule(parent=mol_sys)

        for line in file:
            line_ed = [x for x in line[:-1].split(' ') if x]
            line_ed = []
            i = 0
            x = line[i]
            flag = line[0:6].replace(' ', '')

            if flag == 'ATOM' or flag == 'HETATM':
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=np.float32)
                atom_type = line[76:78].replace(' ', '').capitalize()
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(atom_type), parent=mol,
                                          name=f'{line[6:11].replace(" ", "")}-{line[13:17].replace(" ", "")}',
                                          pdb_flag=flag,
                                          pdb_atom_seq=int(line[6:11].replace(" ", "")),
                                          pdb_name=line[12:16].replace(" ", ""),
                                          pdb_alt_loc=line[16],
                                          pdb_res_name=line[17:20],
                                          pdb_chain=line[21],
                                          pdb_res_seq=int(line[22:26]),
                                          pdb_iCode=line[26],
                                          pdb_occupancy=float(line[54:60]),
                                          pdb_tempFactor=float(line[60:66]),
                                          pdb_charge=line[78:80])
                mol.addChild(atom)
        if bond:
            mol_sys.genBonds()
        if root is not None:
            mol_list = point_class.PointsList(parent=root, name=os.path.basename(file_path).split('.')[0])
            mol.assignPoint(mol_list)
            MOLECULE_SYSTEMS[mol_list] = mol_sys
            atom_list = point_class.PointsList(parent=mol_list, rad=0.25, name='Atoms')
            for atom in mol.children:
                coord = atom.coord.copy()
                point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                          color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                          name=atom.name,
                                          atom_type=atom.atom_type,
                                          pdb_flag=atom.pdb_flag,
                                          pdb_atom_seq=atom.pdb_atom_seq,
                                          pdb_name=atom.pdb_name,
                                          pdb_alt_loc=atom.pdb_alt_loc,
                                          pdb_res_name=atom.pdb_res_name,
                                          pdb_chain=atom.pdb_chain,
                                          pdb_res_seq=atom.pdb_res_seq,
                                          pdb_iCode=atom.pdb_Achar,
                                          pdb_occupancy=atom.pdb_occupancy,
                                          pdb_tempFactor=atom.tempFactor,
                                          pdb_charge=atom.pdb_charge)
                atom.assignPoint(point)
                point.addProperty('name', atom.name)
                point.addProperty('label', atom.name)
            if bond:
                bonds_l = point_class.PointsList(parent=mol_list, rad=0.1, name='Bonds')
                bonds = []
                for atom in mol:
                    for bond in atom.bonds():
                        if bond not in bonds:
                            bond_l = point_class.PointsList(parent=bonds_l,
                                                            name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}',
                                                            rad=bonds_l)
                            b1 = point_class.Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            b2 = point_class.Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            bonds.append(bond)
        return mol_sys, (mol_list, atom_list, bonds_l)


PARSER = FileParser()
