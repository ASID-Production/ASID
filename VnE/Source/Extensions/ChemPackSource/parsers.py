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
                                          color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)])
                atom.assignPoint(point)
                point.addProperty('name', f'{PALETTE.getName(atom.atom_type)}{i}')
                point.addProperty('label', f'{PALETTE.getName(atom.atom_type)}{i}')
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
            flag = ''
            while x.isalpha():
                flag += x
                i += 1
                x = line[i]

            if flag == 'ATOM' or flag == 'HETATM':
                coord = np.array([float(line[26:38]), float(line[38:46]), float(line[46:54])], dtype=np.float32)
                atom_type = line[66:78].replace(' ', '').capitalize()
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(atom_type), parent=mol, name=f'{line[i:11].replace(" ", "")}-{line[13:17].replace(" ", "")}')
                point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                          color=PALETTE.point_dict[atom_type])
                atom.assignPoint(point)
                mol.addChild(atom)
                point.addProperty('name', f'{line[i:11].replace(" ", "")}-{line[13:17].replace(" ", "")}')
                point.addProperty('label', point.name)
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
                                          color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)])
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
