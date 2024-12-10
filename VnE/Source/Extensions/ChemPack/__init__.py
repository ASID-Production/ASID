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


from ... import point_class
import numpy as np
from . import AtomsPalette
import os
import types

import debug

PALETTE = AtomsPalette.Palette()

MOLECULE_SYSTEMS = {}

PARSER = None

TREE_MODEL = None

UNIFORM_MODEL = None

MAIN_WIDGET = None

MAIN_MENU = None


def loadMolSys(molsys, points_lists):
    if TREE_MODEL is not None:
        uniform_model = UNIFORM_MODEL
        model = TREE_MODEL
        if len(points_lists[1].children) > 0:
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
        from .parsers import PARSER as parser
        PARSER = parser
    if TREE_MODEL is not None:
        if root is None:
            root = TREE_MODEL.getRoot()
        if os.path.basename(file) == 'paths.pdb' or os.path.basename(file) == 'CPs.pdb':
            bond = False
        ret = PARSER.parsFile(file, bond=bond, root=root)
        if isinstance(ret, types.GeneratorType):
            for molsys, lists in ret:
                if molsys and lists:
                    loadMolSys(molsys, lists)
            return molsys, lists
        else:
            if ret[0]:
                loadMolSys(*ret)
            return ret
    else:
        ret = PARSER.parsFile(file)
        if isinstance(ret, types.GeneratorType):
            return list(ret)
        else:
            return ret


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
    from . import contacts

    mol_sys = None
    for mol_l in MOLECULE_SYSTEMS:
        if mol_l.pick == 1.0:
            mol_sys = MOLECULE_SYSTEMS[mol_l]
    if mol_sys is not None:
        contacts.test(mol_sys)


def find_sub():
    from . import DrawerWidget
    DrawerWidget.show()
    return


def DbSearch():
    from . import Db_viewer
    Db_viewer.show()
    return


def SymOp():
    from . import Sym_op
    Sym_op.show()
    return


def exportDataF():
    from . import exportData as exportData
    exportData.execute()
    return


def attachCPprop():
    from . import attach_cpprop
    attach_cpprop.execute()
    return


def parseMultiWfn():
    from PySide6.QtWidgets import QFileDialog

    global PARSER
    global TREE_MODEL
    if PARSER is None:
        from .parsers import PARSER as parser
        PARSER = parser
    filename, _ = QFileDialog.getOpenFileName(caption='Molecule file', filter='*.pdb')
    if not filename:
        pass
    else:
        molsys, points_lists = None, None
        if TREE_MODEL is not None:
            molsys, points_lists = PARSER.parsPdb(filename, bond=True, root=TREE_MODEL.getRoot())
            loadMolSys(molsys, points_lists)
        else:
            molsys, _ = PARSER.parsPdb(filename)

    filename, _ = QFileDialog.getOpenFileName(caption='Critical points file', filter='*.pdb')
    if not filename:
        pass
    else:
        cp_molsys, points_lists = None, None
        if TREE_MODEL is not None:
            cp_molsys, cp_points_lists = PARSER.parsPdb(filename, bond=False, root=TREE_MODEL.getRoot())
            cp_points_lists[1].rad = 0.15
            loadMolSys(cp_molsys, cp_points_lists)
        else:
            cp_molsys, _ = PARSER.parsPdb(filename, bond=False)

    filename, _ = QFileDialog.getOpenFileName(caption='Paths file', filter='*.pdb')
    if not filename:
        pass
    else:
        molsys, points_lists = None, None
        if TREE_MODEL is not None:
            molsys, points_lists = PARSER.parsPdb(filename, bond=False, root=TREE_MODEL.getRoot())
            points_lists[1].rad = 0.05
            loadMolSys(molsys, points_lists)
        else:
            molsys, _ = PARSER.parsPdb(filename, bond=False)

    from . import attach_cpprop

    attach_cpprop.execute(cp_points_lists)

    return molsys, points_lists


def parseWinxpro():
    from PySide6.QtWidgets import QFileDialog

    global PARSER
    global TREE_MODEL
    if PARSER is None:
        from .parsers import PARSER as parser
        PARSER = parser
    filename, _ = QFileDialog.getOpenFileName(filter='*.out')
    if not filename:
        return None, None
    molsys, points_lists = None, None
    if TREE_MODEL is not None:
        gen = PARSER.parsWinxproOut(filename, bond=True, root=TREE_MODEL.getRoot())
        for molsys, points_lists in gen:
            loadMolSys(molsys, points_lists)
    else:
        molsys, _ = PARSER.parsWinxproOut(filename)

    filename, _ = QFileDialog.getOpenFileName(filter='*.bpc')
    if not filename:
        return None, None

    if TREE_MODEL is not None:
        gen = PARSER.parsWinxproPaths(filename, bond=False, root=TREE_MODEL.getRoot())
        for molsys, points_lists in gen:
            loadMolSys(molsys, points_lists)
    else:
        molsys, _ = PARSER.parsWinxproPaths(filename)
    return molsys, points_lists

def applyClsf():
    from . import applyCls
    applyCls.execute()

def parsAIMALLsumviz():
    from PySide6.QtWidgets import QFileDialog

    global PARSER
    global TREE_MODEL
    if PARSER is None:
        from .parsers import PARSER as parser
        PARSER = parser
    filename, _ = QFileDialog.getOpenFileName(filter='*.sumviz')
    if not filename:
        return None, None
    molsys, points_lists = None, None
    if TREE_MODEL is not None:
        gen = PARSER.parsAIMALLsumviz(filename, bond=True, root=TREE_MODEL.getRoot())
        for molsys, points_lists in gen:
            if type(molsys) is dict:
                TREE_MODEL.insertRow(TREE_MODEL.rowCount())
            else:
                loadMolSys(molsys, points_lists)
    else:
        molsys, _ = PARSER.parsWinxproOut(filename)
    return molsys, points_lists


def export2dDiagram():
    from . import Db_viewer

    Db_viewer.exportGif()


def symmPOSCAR():
    from . import symm_poscar
    symm_poscar.execute(pars)


def execute(model, uniform_model=None):
    from PySide6.QtWidgets import QFileDialog
    global PARSER
    if PARSER is None:
        from .parsers import PARSER as parser
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
    from .save_file import SAVE_FILE
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

def assemble():
    from . import assemble_cif
    assemble_cif.execute()

def setup(menu, model, uniform_model=None, *args, main_widget=None, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    global UNIFORM_MODEL
    global MAIN_WIDGET
    global MAIN_MENU
    TREE_MODEL = model
    UNIFORM_MODEL = uniform_model
    MAIN_WIDGET = main_widget
    MAIN_MENU = kwargs.get('main_menu', None)

    createPalette()

    from . import calcWidget

    cmenu = menu.addMenu('ChemPack')

    action_test = QAction('Find Sub')
    action_test.triggered.connect(find_sub)
    cmenu.addAction(action_test)

    from . import Db_viewer
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
    action_export.triggered.connect(exportDataF)
    cmenu.addAction(action_export)

    qtaim = cmenu.addMenu('QTAIM tools')

    action_winx = QAction('Parse WinXPRO')
    action_winx.triggered.connect(parseWinxpro)
    qtaim.addAction(action_winx)

    action_cls = QAction('Apply WinXPRO .cls')
    action_cls.triggered.connect(applyClsf)
    qtaim.addAction(action_cls)

    action_aimall = QAction('Parse AIMALL sumviz')
    action_aimall.triggered.connect(parsAIMALLsumviz)
    qtaim.addAction(action_aimall)

    action_multiwfn = QAction('Parse Multiwfn')
    action_multiwfn.triggered.connect(parseMultiWfn)
    qtaim.addAction(action_multiwfn)

    action_2d_export = QAction('Export 2d diagram')
    action_2d_export.triggered.connect(export2dDiagram)
    cmenu.addAction(action_2d_export)

    action_symm_poscar = QAction('Symm POSCAR')
    action_symm_poscar.triggered.connect(symmPOSCAR)
    cmenu.addAction(action_symm_poscar)

    action_assemble = QAction('Assemble molecule')
    action_assemble.triggered.connect(assemble)
    cmenu.addAction(action_assemble)

    actions = [open_action, action_test, action_DB, save_action, action_sym_op, action_export, action_winx, action_cls, action_aimall, action_2d_export, action_multiwfn, action_symm_poscar, action_assemble]
    return actions
