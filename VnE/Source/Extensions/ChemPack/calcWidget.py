from ..ChemPack import MAIN_WIDGET, MOLECULE_SYSTEMS, TREE_MODEL
from . import contacts
from PySide6.QtWidgets import QLabel, QFrame, QSizePolicy

class CalcWidget(QLabel):

    def __init__(self, *args, **kwargs):
        QLabel.__init__(self, *args, **kwargs)
        self.sel = {}
        TREE_MODEL.item_selected.connect(self.add)
        TREE_MODEL.item_deselected.connect(self.remove)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)

    def calc(self):
        atoms = []
        for atms in self.sel.values():
            atoms += atms

        if len(atoms) == 4:
            angle = contacts.angle(atoms, False)
            line = f'{atoms[0].name}--{atoms[1].name}--{atoms[2].name}--{atoms[3].name}: {round(angle, 1)}'
        elif len(atoms) == 3:
            angle = contacts.angle(atoms, False)
            line = f'{atoms[0].name}--{atoms[1].name}--{atoms[2].name}: {round(angle, 1)}'
        elif len(atoms) == 2:
            dist = contacts.dist(*atoms)
            line = f'{atoms[0].name}--{atoms[1].name}: {round(dist, 3)}'
        else:
            line = [str(x.name) for x in atoms]
            line = '--'.join(line)
        self.setText(line)

    def add(self, index):
        atoms = []
        point = index.internalPointer()
        for mol_sys in MOLECULE_SYSTEMS.values():
            atoms += mol_sys.findProp('_point', point)
        atoms_with_coords = []
        for atom in atoms:
            if atom.coord is not None:
                atoms_with_coords.append(atom)
        if atoms_with_coords:
            self.sel[point] = atoms_with_coords
            self.calc()

    def remove(self, index):
        point = index.internalPointer()
        res = self.sel.pop(point, None)
        if res is not None:
            self.calc()


opengl_frame = MAIN_WIDGET.findChild(QFrame, 'opengl_frame')
label = CalcWidget()
opengl_frame.layout().addWidget(label)
a = 0