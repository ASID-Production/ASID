from VnE.Source.Extensions.ChemPack.searchProcess import parser

DIALOG = None

def execute():
    from ..ChemPack import MOLECULE_SYSTEMS
    from ..ChemPack import TREE_MODEL
    from .ui import select_mol_dialog
    from PySide6.QtWidgets import QFileDialog
    from .MoleculeClass import Atom, MoleculeSystem, Molecule
    import os
    import numpy as np
    import cpplib
    from . import parsers
    from . import loadMolSys

    mol_sys = None

    def payload(mol_tuple):

        def fracToDec(a, b, c, al, be, ga, coords):
            al = (al/180)*np.pi
            be = (be/180)*np.pi
            ga = (ga/180)*np.pi

            sin = np.sin
            cos = np.cos
            cot = lambda x: np.tan(x)**-1
            csc = lambda x: np.sin(x)**-1

            mat = np.array([[a*sin(be)*np.sqrt(1-(cot(al)*cot(be) - csc(al)*csc(be)*cos(ga))**2), 0, 0],
                            [a*csc(al)*cos(ga) - a*cot(al)*cos(be), b*sin(al), 0],
                            [a*cos(be), b*cos(al), c]])
            mat = mat.transpose()
            for i in range(len(coords)):
                coords[i] = (coords[i] @ mat).astype(dtype=np.float32)
            return coords

        point_list, mol_sys = mol_tuple
        mol = mol_sys.children[0]
        atom_list = None
        for l in point_list.children:
            if l.name == 'Atoms':
                atom_list = l
                break

        if atom_list is None:
            return

        file, _ = QFileDialog.getOpenFileName(caption='WinXPRO *.cls file', filter='*.cls')
        new_mol_sys = MoleculeSystem()
        new_mol_sys.file_name = file
        new_mol_sys.open_func = lambda: payload(mol_tuple)
        new_mol_sys.name = os.path.basename(file).split('.')[0]
        new_mol = Molecule(new_mol_sys)
        new_atoms = []

        for atom in mol:
            props = {x: y for x, y in atom.__dict__.items() if x[0] != '_'}
            new_atom = Atom(props.pop('coord'), props.pop('atom_type'), new_mol, **props)
            new_atoms.append(new_atom)


        file = open(file, 'r')
        for line in file:
            if 'symm.       translations        appl.for atoms  No.in list' in line:
                for line in file:
                    sym = []
                    data = [x for x in line.split(' ') if x]
                    if len(data) > 3:
                        atom = mol.children[int(data[4]) - 1]
                        props = {x: y for x, y in atom.__dict__.items() if x[0] != '_'}
                        new_atom = Atom(props.pop('coord'), props.pop('atom_type'), new_mol, **props)
                        new_atoms.append(new_atom)
                        new_atom.name = data[6]
                        shift = [lambda n: f'x{n:+},y,z', lambda n: f'x,y{n:+},z', lambda n: f'x,y,z{n:+}']
                        sym.append(new_atom.cif_sym_codes[int(data[0])-1][1])
                        for tr in data[1:4]:
                            num = int(float(tr))
                            if num == 0:
                                shift.pop(0)
                                continue
                            else:
                                sym.append(shift.pop(0)(num))
                        frac = (new_atom.atom_type, *new_atom.cif_frac_coords)
                        for s in sym:
                            frac = cpplib.GenSymm([frac], 0b0, [s])
                            if len(frac) != 1:
                                frac = frac[1]
                            else:
                                frac = frac[0]
                        new_atom.cif_frac_coords = frac[1:]

        cell = [new_atoms[0].cif_cell_a,
                new_atoms[0].cif_cell_b,
                new_atoms[0].cif_cell_c,
                new_atoms[0].cif_cell_al,
                new_atoms[0].cif_cell_be,
                new_atoms[0].cif_cell_ga]

        cell.append([x.cif_frac_coords for x in new_atoms])
        new_coords = fracToDec(*cell)
        for i, atm in enumerate(new_atoms):
            atm.coord = new_coords[i]

        ret = parsers.PARSER.parsMolSys(new_mol_sys, True, root=TREE_MODEL.getRoot())
        loadMolSys(*ret)

    global DIALOG
    DIALOG = select_mol_dialog.SelectMolDialog(MOLECULE_SYSTEMS, payload)
    DIALOG.show()
