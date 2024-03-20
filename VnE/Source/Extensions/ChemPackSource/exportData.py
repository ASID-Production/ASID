def execute():
    from ..ChemPack import MOLECULE_SYSTEMS
    from PySide6.QtWidgets import QFileDialog
    import cpplib
    import csv

    mol_sys = None

    for mol_l in MOLECULE_SYSTEMS:
        if mol_l.isValid() and mol_l.pick:
            mol_sys = MOLECULE_SYSTEMS[mol_l]
            break

    if mol_sys is not None:
        filepath = QFileDialog.getSaveFileName(filter='*.csv')
        if not filepath[0]:
            return
        atoms = mol_sys.children[0].children
        arg = [[x.atom_type, *list(x.coord)] for x in atoms]
        ret_dist = [[int(x) if i != 2 else float(x) for i, x in enumerate(pair.split(':'))] for pair in cpplib.GenBondsEx(arg).split('\n')[:-1]]
        args = [
            [x.atom_type for x in atoms],
            [list(x.coord) for x in atoms],
            [0, 0, 0, 0.0, 0.0, 0.0, -180.0, 180]
        ]
        #angles = cpplib.FindAngleWC(*args)
        with open(filepath[0], 'w', newline='') as out:
            writer = csv.writer(out, delimiter=';')
            lines = []
            for pair in ret_dist:
                line = [f'{atoms[pair[0]].name}--{atoms[pair[1]].name}', str(round(pair[2], 4))]
                lines.append(line)
            writer.writerow(['Pair', 'Dist'])
            writer.writerows(lines)