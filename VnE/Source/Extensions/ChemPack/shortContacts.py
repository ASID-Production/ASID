#  Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
#  This file is part of ASID - Atomistic Simulation Instruments and Database
#  For more information see <https://github.com/ASID-Production/ASID>
#  #
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  #
#      http://www.apache.org/licenses/LICENSE-2.0
#  #
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#  #
#  ******************************************************************************************
#   Author:      Alexander A. Korlyukov (head)
#   ORCID:       0000-0002-5600-9886
#   Author:      Alexander D. Volodin (author of cpplib)
#   ORCID:       0000-0002-3522-9193
#   Author:      Petr A. Buikin (author of api_database)
#   ORCID:       0000-0001-9243-9915
#   Author:      Alexander R. Romanenko (author of VnE)
#   ORCID:       0009-0003-5298-6836
#  #
#  ******************************************************************************************

from ..ChemPack.MoleculeClass import MoleculeSystem, Molecule, Atom
from ..ChemPack.ui.select_mol_dialog import SelectMolDialog
from ...point_class import PointsList, Point
import os.path as opath
import numpy as np

class ContMol:
    def __init__(self, mol: Molecule, symm, ind):
        self.mol = mol
        self.symm = symm
        self.atoms = []
        self.point_list = None
        self.ind = ind

        self.drain_contacts = []
        self.source_contacts = []

        for a in self.mol.children:
            self.atoms.append(ContAtom(a, self.symm, self))
        self.visible = False

    def vis(self, parent: PointsList, scene):
        from ..ChemPack.parsers import FileParser
        from ..ChemPack import TREE_MODEL

        scene.visible_mols.append(self)
        mol_list = PointsList(parent=parent, name=str(self.symm))
        self.mol.genBonds()
        self.mol.assignPoint(mol_list)
        atom_list = PointsList(parent=mol_list, rad=0.15, name='Atoms')
        i = 1
        for atom in self.mol.children:
            point = FileParser._pointCreation(atom_list, atom)
            atom.assignPoint(point)
            i += 1
        bonds_l = PointsList(parent=mol_list, rad=0.05, name='Bonds')
        bond_atm = []
        for atom in self.mol:
            for bond in atom.bonds():
                atm = bond.get(atom)
                if atm not in bond_atm:
                    bond_l = PointsList(parent=bonds_l,
                                        name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}',
                                        rad=bonds_l)
                    b1 = Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                               rad=bond_l,
                               parent=bond_l)
                    b2 = Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                               rad=bond_l,
                               parent=bond_l)
                    bond.assignPoint((b1, b2))
            bond_atm.append(atom)
        atoms_index = TREE_MODEL.index(0, 0, by_point=atom_list)
        TREE_MODEL.attachObserver(atoms_index, 'Sphere')
        atoms_index = TREE_MODEL.index(0, 0, by_point=bonds_l)
        TREE_MODEL.attachObserver(atoms_index, 'Bond')
        self.visible = True


class ContAtom:
    def __init__(self, atm: Atom, symm, parent: ContMol=None):
        self.atom = atm
        self.symm = symm
        self.parent = parent
        self.point = None

class Contact:
    def __init__(self, source: ContAtom, dest: ContAtom):
        self.source = source
        self.dest = dest
        self.dist = np.linalg.norm(self.dest.atom.coord - self.source.atom.coord)
        self.point_list = None
        self.dp, self.cps, self.cpd = None, None, None
        self.visible = False

    def vis(self, parent: PointsList=None):
        from ..ChemPack import PALETTE, TREE_MODEL
        if self.point_list:
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.dp)
            TREE_MODEL.attachObserver(atoms_index, 'Sphere')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cps)
            TREE_MODEL.attachObserver(atoms_index, 'Bond')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cpd)
            TREE_MODEL.attachObserver(atoms_index, 'Bond')
            self.visible = True
            return
        else:
            self.point_list = PointsList(parent=parent, name=f'{self.source.atom.name}:{self.source.symm}--{self.dest.atom.name}:{self.dest.symm} = {self.dist}')
            self.dp = Point(parent=self.point_list, color=np.array([*PALETTE.getColor(self.dest.atom.atom_type), 255], dtype=np.float32)/255, name=f'{self.dest.atom.name}:{self.dest.symm}', rad=0.10, coord=self.dest.atom.coord.copy())
            self.cps = Point(parent=self.point_list, color=np.array([*PALETTE.getColor(self.source.atom.atom_type), 255], dtype=np.float32)/255, name=f'{self.source.atom.name}:{self.source.symm}--{self.dest.atom.name}:{self.dest.symm}', rad=0.025, freq=5, hfreq=1, coord=self.source.atom.coord.copy())
            self.cpd = Point(parent=self.point_list, color=np.array([*PALETTE.getColor(self.dest.atom.atom_type), 255], dtype=np.float32)/255, name=f'{self.source.atom.name}:{self.source.symm}--{self.dest.atom.name}:{self.dest.symm}', rad=0.025, freq=5, hfreq=1, coord=self.dest.atom.coord.copy())
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.dp)
            TREE_MODEL.attachObserver(atoms_index, 'Sphere')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cps)
            TREE_MODEL.attachObserver(atoms_index, 'Bond')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cpd)
            TREE_MODEL.attachObserver(atoms_index, 'Bond')
            self.visible = True


    def hide(self):
        from ..ChemPack import TREE_MODEL
        if self.visible:
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.dp)
            TREE_MODEL.detachObserver(atoms_index, 'Sphere')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cps)
            TREE_MODEL.detachObserver(atoms_index, 'Bond')
            atoms_index = TREE_MODEL.index(0, 0, by_point=self.cpd)
            TREE_MODEL.detachObserver(atoms_index, 'Bond')
            self.visible = False

    def delete(self, scene):
        if self.point_list:
            from ..ChemPack import TREE_MODEL
            self.point_list.__del__()
            self.point_list = None
            TREE_MODEL.update()
        try:
            scene.contacts.remove(self)
        except ValueError:
            return

class Scene:

    def __init__(self, center: list[ContMol], parent_point_list, mol_sys: MoleculeSystem):
        self.cent = center
        self.surr_mols = {}
        self.contacts = []
        self.visible_mols = []
        self.points_list = PointsList(parent_point_list, name='Short contacts')
        self.mols_list = PointsList(self.points_list, name='Mols')
        self.contacts_list = PointsList(self.points_list, name='Contacts')
        self.mol_sys = mol_sys

        _ = ['cif_cell_a', 'cif_cell_b', 'cif_cell_c', 'cif_cell_al', 'cif_cell_be', 'cif_cell_ga']
        self.cell = [float(getattr(self.cent[0].mol.children[0], x)) for x in _]
        self.symms = [x[1] for x in self.cent[0].mol.children[0].cif_sym_codes]
        self.ver_radius = 5.0

    def genContacts(self, mol: ContMol, r=5.0):
        import cpplib
        from ..ChemPack.parsers import PARSER
        tuples = [[(at.atom_type, *at.cif_frac_coords) for at in y.mol.children] for y in self.cent]
        coords = [(*at.cif_frac_coords, r) for at in mol.mol.children]
        self.mol_sys.addChild(mol.mol)
        mol.vis(self.mols_list, self)
        for ind, tm in enumerate(tuples):
            ret = cpplib.Cluster(self.cell, self.symms, tm, coords, self.ver_radius)
            exclude = self.visible_mols

            dec = PARSER.fracToDec(*self.cell, [x['point_frac'][:3] for x in ret['points']])
            # dec = [x['point_cart'][:3] for x in ret['points']]
            new_mol = Molecule()
            cif = ['cif_space_group', 'cif_sym_codes', 'cif_cell_a', 'cif_cell_b', 'cif_cell_c', 'cif_cell_al',
                   'cif_cell_be', 'cif_cell_ga']
            cif = {x: getattr(self.cent[ind].mol.children[0], x) for x in cif}
            for i, nat in enumerate(ret['points']):
                if nat['shift'] == (0, 0, 0) and nat['symmref'] == 0:
                    at = self.cent[ind].mol.children[nat['index']]
                    cif['cif_frac_coords'] = at.cif_frac_coords.copy()
                    new_mol.addChild(Atom(at.coord.copy(), at.atom_type, name=at.name, creation_code=(0, 0, 0, 0), **cif))
                else:
                    cif['cif_frac_coords'] = np.array(nat['point_frac'])
                    new_mol.addChild(Atom(dec[i], self.cent[ind].mol.children[nat['index']].atom_type,
                                          name=f'{self.cent[ind].mol.children[nat["index"]].name}_{nat["symmref"]},{nat["shift"][0]},{nat["shift"][1]},{nat["shift"][2]}',
                                          creation_code=(*[nat['symmref'], *nat["shift"]],), **cif, cif_uniq=False))
            new_mol.genBonds()
            new_mols = new_mol.splitMol()
            mol_surr = []
            surr_codes = []
            surr = []
            for m in self.surr_mols:
                surr_codes += [(x.ind, *x.symm) for x in self.surr_mols[m]]
                surr += self.surr_mols[m]
            surr_codes += [(x,0,0,0,0) for x in range(len(self.cent))]
            surr += [x for x in self.cent]
            for m in new_mols:
                code = (ind, *m.children[0].creation_code)
                if code == (mol.ind, *mol.symm):
                    continue
                elif code in surr_codes:
                    i = surr_codes.index(code)
                    mol_surr.append(surr[i])
                else:
                    mol_surr.append(ContMol(m, m.children[0].creation_code, ind))
            self.surr_mols[mol] = mol_surr
            for c in mol.source_contacts:
                c.source.parent.drain_contacts.remove(c)
                c.delete(self)
            for m in mol_surr:
                if not m in exclude:
                    self.createContacts(mol, m)
    def createContacts(self, source: ContMol, dest: ContMol):
        for a in source.atoms:
            for ad in dest.atoms:
                c = Contact(a, ad)
                if c.dist > 5:
                    continue
                self.contacts.append(c)
                source.drain_contacts.append(c)
                dest.source_contacts.append(c)

def execute():
    from ..ChemPack import MOLECULE_SYSTEMS, TREE_MODEL
    from ..ChemPack import loadMolSys
    from PySide6.QtWidgets import QFrame, QLabel, QLineEdit, QHBoxLayout, QPushButton
    scene = None
    def process(molsys=None):
        nonlocal scene
        molsys = DIALOG.curr_sys
        list_obj, molsys = molsys
        try:
            r = float(DIALOG.rad.text())
        except ValueError:
            r = 2.5
        if not scene:
            mols = [ContMol(x, (0,0,0,0), i) for i, x in enumerate(molsys.children[0].splitMol())]
            new_molsys = MoleculeSystem()
            point_list = PointsList(parent=TREE_MODEL.getRoot(), name='test')
            TREE_MODEL.update()
            MOLECULE_SYSTEMS[point_list] = new_molsys
            scene = Scene(mols, point_list, new_molsys)
            scene.genContacts(mols[0], r)
            for c in scene.contacts:
                if c.dist < float(text_edit2.text()):
                    c.vis(scene.contacts_list)
            ...
        else:
            for c in scene.contacts:
                if c.dp and c.dp.pick:
                    if not c.dest.parent.visible:
                        c.dest.parent.mol.genBonds()
                        scene.genContacts(c.dest.parent, r)
            for c in scene.contacts:
                if c.dist < float(text_edit2.text()):
                    c.vis(scene.contacts_list)
        return
    def setRad():
        for c in scene.contacts:
            if c.dist < float(text_edit2.text()):
                c.vis(scene.contacts_list)
            else:
                c.hide()
    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)

    frame = QFrame()
    frame2 = QFrame()

    layout = QHBoxLayout()
    layout2 = QHBoxLayout()
    label = QLabel("Search radius:")
    label2 = QLabel("Visible radius:")
    text_edit = QLineEdit()
    text_edit2 = QLineEdit()
    button = QPushButton(text='Generate contact')
    button2 = QPushButton(text='Set visible radius')
    layout.addWidget(label)
    layout.addWidget(text_edit)
    layout2.addWidget(label2)
    layout2.addWidget(text_edit2)
    frame.setLayout(layout)
    frame2.setLayout(layout2)
    DIALOG.rad = text_edit

    button.pressed.connect(process)
    button2.pressed.connect(setRad)

    DIALOG.main_layout.insertWidget(DIALOG.main_layout.count() - 1, frame)
    DIALOG.main_layout.insertWidget(DIALOG.main_layout.count() - 1, frame2)
    DIALOG.main_layout.insertWidget(DIALOG.main_layout.count() - 1, button2)
    DIALOG.main_layout.insertWidget(DIALOG.main_layout.count() - 1, button)
    DIALOG.show()
    pass

def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Short contacts')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions