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
# ******************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# ******************************************************************************************


from .. import point_class
import numpy as np
from .ChemPackSource import AtomsPalette
import os

PALETTE = AtomsPalette.Palette()

MOLECULE_SYSTEMS = {}

PARSER = None

TREE_MODEL = None


def pars(file, bond=True, root=None):
    '''from .ChemPackSource import MoleculeClass

    def pars_xyz(bond=True):
        i = 1
        for line in file:
            line_ed = [x for x in line[:-1].split(' ') if x]
            if len(line_ed) == 4:
                coord = np.array([float(x) for x in line_ed[1:]], dtype=np.float32)
                if line_ed[0].isdigit():
                    atom = MoleculeClass.Atom(coord.copy(), int(line_ed[0]), parent=mol)
                else:
                    atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(line_ed[0]), parent=mol)
                coord[2] -= 500
                if line_ed[0].isdigit():
                    line_ed[0] = PALETTE.getName(line_ed[0])
                point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                          color=PALETTE.point_dict[line_ed[0]])
                atom.assignPoint(point)
                mol.addChild(atom)
                point.addProperty('name', f'{line_ed[0]}{i}')
                point.addProperty('label', f'{line_ed[0]}{i}')
                i += 1
        if bond:
            MOLECULE_SYSTEMS[molecule_list].genBonds()
        bonds = []
        for atom in mol:
            for bond in atom.bonds():
                if bond not in bonds:
                    bond_l = point_class.PointsList(parent=bond_list)
                    bond_l.addProperty('name', f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}')
                    bond_l.addProperty('rad', bond_list)
                    b1 = point_class.Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(), rad=bond_l,
                                           parent=bond_l)
                    b2 = point_class.Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(), rad=bond_l,
                                           parent=bond_l)
                    bonds.append(bond)
        return

    def pars_pdb(bond=True):
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
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(atom_type), parent=mol)
                coord[2] -= 500
                point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                          color=PALETTE.point_dict[atom_type])
                atom.assignPoint(point)
                mol.addChild(atom)
                point.addProperty('name', f'{line[i:11].replace(" ", "")}-{line[13:17].replace(" ", "")}')
                point.addProperty('label', point.name)
        if bond:
            MOLECULE_SYSTEMS[molecule_list].genBonds()
        bonds = []
        for atom in mol:
            for bond in atom.bonds():
                if bond not in bonds:
                    bond_l = point_class.PointsList(parent=bond_list)
                    bond_l.addProperty('name', f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}')
                    bond_l.addProperty('rad', bond_list)
                    b1 = point_class.Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(), rad=bond_l,
                                           parent=bond_l)
                    b2 = point_class.Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(), rad=bond_l,
                                           parent=bond_l)
                    bonds.append(bond)
        return

    global MOLECULE_SYSTEMS
    molecule_list = point_class.PointsList(parent=TREE_MODEL.getRoot())
    MOLECULE_SYSTEMS[molecule_list] = MoleculeClass.MoleculeSystem()
    molecule_list.addProperty('name', name)
    atom_list = point_class.PointsList(parent=molecule_list)
    atom_list.addProperty('name', 'Atoms')
    atom_list.addProperty('rad', 0.25)
    bond_list = point_class.PointsList(parent=molecule_list)
    bond_list.addProperty('name', 'Bonds')
    bond_list.addProperty('rad', 0.1)

    mol = MoleculeClass.Molecule()
    MOLECULE_SYSTEMS[molecule_list].addChild(mol)
    mol.assignPoint(molecule_list)
    if ext == 'xyz':
        pars_xyz(bond)
    elif ext == 'pdb':
        pars_pdb(bond)
    return molecule_list, (atom_list, bond_list)'''
    global PARSER
    global TREE_MODEL
    if PARSER is None:
        from .ChemPackSource.parsers import PARSER as parser
        PARSER = parser
    if TREE_MODEL is not None:
        if os.path.basename(file) == 'paths.pdb' or os.path.basename(file) == 'CPs.pdb':
            bond = False
        molsys, point_lists = PARSER.parsFile(file, bond=bond, root=root)
        return molsys, point_lists
    else:
        molsys, _ = PARSER.parsFile(file)
        return molsys, None


def createPalette():
    z = False
    for child in TREE_MODEL.getRoot().children:
        if child.name == 'Atom Colors':
            pass
            z = True
            break
    if z:
        pass
    else:
        global PALETTE
        color_pal = point_class.PointsList(parent=TREE_MODEL.getRoot())
        color_pal.addProperty('name', 'Atom Colors')
        for atom in PALETTE.palette:
            color = np.array(PALETTE.getColor(atom) + [255], dtype=np.float32) / 255
            point = point_class.Point(parent=color_pal, color=color)
            point.addProperty('name', atom)
            PALETTE.addPoint(atom, point)


def test_cont():
    from .ChemPackSource import contacts

    mol_sys = None
    for mol_l in MOLECULE_SYSTEMS:
        if mol_l.pick == 1.0:
            mol_sys = MOLECULE_SYSTEMS[mol_l]
    if mol_sys is not None:
        contacts.test(mol_sys)


def find_sub():
    from .ChemPackSource import DrawerWidget
    DrawerWidget.show()
    return


def DbSearch():
    from .ChemPackSource import Db_viewer
    Db_viewer.show()
    return


def execute(model):
    from PySide6.QtWidgets import QFileDialog
    filter = '*.xyz *.pdb'
    filename = QFileDialog.getOpenFileName(filter=filter)
    if type(filename) == tuple:
        filename = filename[0]

    if filename != '':

        _, points_lists = pars(filename, root=TREE_MODEL.getRoot())
        if _ is not None:
            model.insertRow(model.rowCount())
            if points_lists is not None:
                if points_lists[1] is not None:
                    atoms_index = model.index(0, 0, by_point=points_lists[1])
                    model.attachObserver(atoms_index, 'Sphere')
                if points_lists[2] is not None:
                    bonds_index = model.index(0, 0, by_point=points_lists[2])
                    model.attachObserver(bonds_index, 'Bond')
    return 0


def save():
    from .ChemPackSource.save_file import SAVE_FILE
    from PySide6.QtWidgets import QFileDialog
    import os

    formats = SAVE_FILE.getFormats()
    filters = f'*.{" *.".join([x for x in formats])}'
    filename = QFileDialog.getSaveFileName(filter=filters)

    if type(filename) == tuple:
        filename = filename[0]

    if filename != '':
        format = os.path.basename(filename).split('.')[1]
        mol_sys = None
        for mol_l in MOLECULE_SYSTEMS:
            if mol_l.pick == 1.0 and mol_l.isValid():
                mol_sys = MOLECULE_SYSTEMS[mol_l]
                break
        if mol_sys is not None:
            SAVE_FILE.save(mol_sys, filename, format)
    return


def setup(menu, model):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    createPalette()

    cmenu = menu.addMenu('ChecmPack')

    action_test = QAction('Find Sub')
    action_test.triggered.connect(find_sub)
    cmenu.addAction(action_test)

    action_DB = QAction('DB Search')
    action_DB.triggered.connect(DbSearch)
    cmenu.addAction(action_DB)

    open_action = QAction('Open')
    open_action.setShortcut('Ctrl+O')
    cmenu.addAction(open_action)
    open_action.triggered.connect(lambda: execute(model))
    save_action = QAction('Save')
    save_action.setShortcut('Ctrl+S')
    cmenu.addAction(save_action)
    save_action.triggered.connect(save)
    actions = [open_action, action_test, action_DB, save_action]
    return actions
