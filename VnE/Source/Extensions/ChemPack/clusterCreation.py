from VnE.Source.Extensions.ChemPack.MoleculeClass import MoleculeSystem, Molecule, Atom
from VnE.Source.Extensions.ChemPack.ui.select_mol_dialog import SelectMolDialog
import os.path as opath

def execute():
    from ..ChemPack import MOLECULE_SYSTEMS
    from ..ChemPack import loadMolSys
    from PySide6.QtWidgets import QFrame, QLabel, QLineEdit, QHBoxLayout

    def process(molsys):

        from ..ChemPack.parsers import PARSER
        import numpy as np
        from ..ChemPack import TREE_MODEL
        import cpplib
        list_obj, molsys = molsys
        try:
            r = float(DIALOG.rad.text())
        except ValueError:
            r = 0.0
        _ = ['cif_cell_a', 'cif_cell_b', 'cif_cell_c', 'cif_cell_al', 'cif_cell_be', 'cif_cell_ga']
        cell = [float(getattr(molsys.children[0].children[0], x)) for x in _]
        symm = [x[1] for x in molsys.children[0].children[0].cif_sym_codes]
        tuples = [(at.atom_type, *at.cif_frac_coords) for at in molsys.children[0]]
        coords = [(*at.cif_frac_coords, r) for at in molsys.children[0]]
        ver_radius = float(10)
        ret = cpplib.Cluster(cell, symm, tuples, coords, ver_radius)
        dec = PARSER.fracToDec(*cell, [x[:3] for x in ret['return']])
        new_mol_sys = MoleculeSystem()
        new_mol_sys.name = molsys.name + ' cluster'
        new_mol = Molecule(new_mol_sys)
        cif = ['cif_space_group', 'cif_sym_codes', 'cif_cell_a', 'cif_cell_b', 'cif_cell_c', 'cif_cell_al', 'cif_cell_be', 'cif_cell_ga']
        cif = {x: getattr(molsys.children[0].children[0], x) for x in cif}
        for i, nat in enumerate(ret['return']):
            if nat[4:] == [0,0,0,0]:
                at = molsys.children[0].children[nat[3]]
                new_mol.addChild(Atom(at.coord, at.atom_type, name=at.name, creation_code=(0, 0, 0, 0), cif_frac_coords=at.cif_frac_coords, **cif))
            new_mol.addChild(Atom(dec[i], molsys.children[0].children[nat[3]].atom_type, name=f'{molsys.children[0].children[nat[3]].name}_{nat[4]},{nat[5]},{nat[6]},{nat[7]}', creation_code=(*nat[4:],), cif_frac_coords=np.array(nat[:3], dtype=np.float32), **cif))
        lists = PARSER.parsMolSys(new_mol_sys, True, TREE_MODEL.getRoot())

        cell_coords = [[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1],
                       [0, 1, 1],
                       [1, 1, 0],
                       [1, 0, 1],
                       [1, 1, 1]]
        cell_dec_coords = PARSER.fracToDec(*cell, cell_coords)

        PARSER.createCellList(lists[1][0], cell_dec_coords)
        loadMolSys(new_mol_sys, lists[1])

        return

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)

    frame = QFrame()

    layout = QHBoxLayout()
    label = QLabel("Cluster radius:")
    text_edit = QLineEdit()
    layout.addWidget(label)
    layout.addWidget(text_edit)
    frame.setLayout(layout)
    DIALOG.rad = text_edit

    DIALOG.main_layout.insertWidget(DIALOG.main_layout.count() - 1, frame)
    DIALOG.show()
    pass

def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Create cluster')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions