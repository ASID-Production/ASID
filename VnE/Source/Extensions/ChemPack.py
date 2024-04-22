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

UNIFORM_MODEL = None


def loadMolSys(molsys, points_lists):
    if TREE_MODEL is not None:
        uniform_model = UNIFORM_MODEL
        model = TREE_MODEL

        coords = np.array([a.coord for a in points_lists[1].children])
        cent = np.sum(coords, axis=0) / len(coords)
        if uniform_model is not None:
            root = uniform_model.root()
            root.rotation_point = cent.copy()
        if molsys is not None:
            model.insertRow(model.rowCount())
            if points_lists is not None:
                if points_lists[1] is not None:
                    atoms_index = model.index(0, 0, by_point=points_lists[1])
                    model.attachObserver(atoms_index, 'Sphere')
                if points_lists[2] is not None:
                    bonds_index = model.index(0, 0, by_point=points_lists[2])
                    model.attachObserver(bonds_index, 'Bond')


def pars(file, bond=True, root=None):
    global PARSER
    global TREE_MODEL
    if PARSER is None:
        from .ChemPackSource.parsers import PARSER as parser
        PARSER = parser
    if TREE_MODEL is not None:
        if os.path.basename(file) == 'paths.pdb' or os.path.basename(file) == 'CPs.pdb':
            bond = False
        molsys, points_lists = PARSER.parsFile(file, bond=bond, root=root)
        loadMolSys(molsys, points_lists)
        return molsys, points_lists
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


def SymOp():
    from .ChemPackSource import Sym_op
    Sym_op.show()
    return


def exportData():
    from .ChemPackSource import exportData
    exportData.execute()
    return

def attachCPprop():
    from .ChemPackSource import attach_cpprop
    attach_cpprop.execute()
    return

def execute(model, uniform_model=None):
    from PySide6.QtWidgets import QFileDialog
    global PARSER
    if PARSER is None:
        from .ChemPackSource.parsers import PARSER as parser
        PARSER = parser
    filter = PARSER.SUPPORTED_FORMATS.keys()
    filter = ' '.join([f'*{x}' for x in filter])
    filename = QFileDialog.getOpenFileName(filter=filter)
    if type(filename) == tuple:
        filename = filename[0]

    if filename != '':
        _, points_lists = pars(filename, root=TREE_MODEL.getRoot())
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


def setup(menu, model, uniform_model=None, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    global UNIFORM_MODEL
    TREE_MODEL = model
    UNIFORM_MODEL = uniform_model

    createPalette()

    cmenu = menu.addMenu('ChemPack')

    action_test = QAction('Find Sub')
    action_test.triggered.connect(find_sub)
    cmenu.addAction(action_test)

    action_DB = QAction('DB Search')
    action_DB.triggered.connect(DbSearch)
    cmenu.addAction(action_DB)

    open_action = QAction('Open')
    open_action.setShortcut('Ctrl+O')
    cmenu.addAction(open_action)
    open_action.triggered.connect(lambda: execute(model, uniform_model))
    save_action = QAction('Save')
    save_action.setShortcut('Ctrl+S')
    cmenu.addAction(save_action)
    save_action.triggered.connect(save)

    action_sym_op = QAction('Symmetry operations')
    action_sym_op.triggered.connect(SymOp)
    cmenu.addAction(action_sym_op)

    action_export = QAction('Export Data')
    action_export.triggered.connect(exportData)
    cmenu.addAction(action_export)

    action_attach = QAction('Attach CPprop')
    action_attach.triggered.connect(attachCPprop)
    cmenu.addAction(action_attach)

    actions = [open_action, action_test, action_DB, save_action, action_sym_op, action_export, action_attach]
    return actions
