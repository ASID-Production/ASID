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

from .ui.select_mol_dialog import SelectMolDialog


def execute():
    from . import MOLECULE_SYSTEMS
    from . import parsers
    from . import point_class
    from . import TREE_MODEL
    import numpy as np

    def process(molsys):

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

        list_obj, molsys = molsys
        import cpplib

        mol = molsys.children[0]
        a = mol.children[0]
        cell = [a.cif_cell_a, a.cif_cell_b, a.cif_cell_c, a.cif_cell_al, a.cif_cell_be, a.cif_cell_ga]
        symm = [x[1] for x in a.cif_sym_codes]
        coords = [tuple([x.atom_type, *x.cif_frac_coords]) for x in mol.children]
        b = []
        c = []
        for x in coords:
            b += x[1:]
            c.append(x[0])
        ret = cpplib.compaq(cell, symm, coords)['xyz_block']
        ncoords = fracToDec(*cell, ret)
        for i, atm in enumerate(mol.children):
            atm.cif_frac_coords = ret[i]
            atm.coord = ncoords[i]
        l = point_class.PointsList()
        lists = parsers.PARSER.parsMolSys(molsys, root=l)[1]
        moll, atoms, bonds = lists
        MOLECULE_SYSTEMS.pop(moll)
        for old_l in list_obj.children:
            if old_l.name == 'Atoms':
                oatm = old_l
            if old_l.name == 'Bonds':
                obonds = old_l

        ind1 = TREE_MODEL.index(0,0, by_point=oatm)
        orow = ind1.row()
        par = ind1.parent()
        TREE_MODEL.removeRow(ind1.row(), parent=par)
        ind2 = TREE_MODEL.index(0,0, by_point=obonds)
        TREE_MODEL.removeRow(ind2.row(), parent=par)
        list_obj.addChild(atoms)
        list_obj.addChild(bonds)
        TREE_MODEL.insertRow(orow, parent=par)
        TREE_MODEL.insertRow(orow, parent=par)
        ind1 = TREE_MODEL.index(0, 0, by_point=atoms)
        ind2 = TREE_MODEL.index(0, 0, by_point=bonds)
        TREE_MODEL.attachObserver(ind1, 'Sphere')
        TREE_MODEL.attachObserver(ind2, 'Bond')
        a = 0



        return

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()
    pass