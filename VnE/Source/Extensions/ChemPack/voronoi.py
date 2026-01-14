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


def execute():
    from . import MAIN_WIDGET, MOLECULE_SYSTEMS, TREE_MODEL
    from ... import point_class
    from PySide6.QtOpenGLWidgets import QOpenGLWidget
    from ..ChemPack.ui.select_mol_dialog import SelectMolDialog
    import numpy as np
    import cpplib

    def fracToDec(a, b, c, al, be, ga, coords):
        al = (al / 180) * np.pi
        be = (be / 180) * np.pi
        ga = (ga / 180) * np.pi

        sin = np.sin
        cos = np.cos
        cot = lambda x: np.tan(x) ** -1
        csc = lambda x: np.sin(x) ** -1

        mat = np.array([[a * sin(be) * np.sqrt(1 - (cot(al) * cot(be) - csc(al) * csc(be) * cos(ga)) ** 2), 0, 0],
                        [a * csc(al) * cos(ga) - a * cot(al) * cos(be), b * sin(al), 0],
                        [a * cos(be), b * cos(al), c]])
        mat = mat.transpose()
        for i in range(len(coords)):
            coords[i] = (coords[i] @ mat).astype(dtype=np.float32)
        return coords

    def process(mol_sys):
        opengl_widget = MAIN_WIDGET.findChild(QOpenGLWidget, "OpenGLWidget")
        sel = opengl_widget.selection_model.selection()
        sel = [x.indexes()[0].internalPointer() for x in sel]
        sel = [x._atom for x in sel if x._atom]
        args = {'cell_params': None,
                'symms': None,
                'atoms': None,
                'bools': None,
                'cutoff': None}

        atoms = mol_sys[1].children[0].children
        sym_codes = atoms[0].cif_sym_codes


        if sym_codes:
            cell = [atoms[0].cif_cell_a,
                    atoms[0].cif_cell_b,
                    atoms[0].cif_cell_c,
                    atoms[0].cif_cell_al,
                    atoms[0].cif_cell_be,
                    atoms[0].cif_cell_ga]
        else:
            ...

        symms = [x[1] for x in sym_codes]
        data = [tuple([x.atom_type, *list(x.cif_frac_coords)]) for x in atoms]
        bools = [True if x in sel else False for x in atoms]

        res = cpplib.VoronoiCalculation(cell, symms, data, bools, 6.0)

        points_list = point_class.PointsList(TREE_MODEL.getRoot(), name=f'Voronoi polyhedron')
        edges_l = point_class.PointsList(points_list, rad=0.01  , color=[0,0,0,1], name='Edges')
        vert_l = point_class.PointsList(points_list, rad=0.01  , color=[0,0,0,1], name='Vertex')
        poly_l = point_class.PointsList(points_list, color=[0.0,1.0,0.0,0.5], name='Polygons')
        for v in res['vertexes']:
            point_class.Point(vert_l, coord=v, color=vert_l, rad=vert_l)
        for i, v in enumerate(res['vertexes']):
            edge_list = point_class.PointsList(edges_l, color=edges_l, rad=edges_l)
            point_class.Point(edge_list, coord=vert_l.children[i], color=edge_list, rad=edge_list)
            point_class.Point(edge_list, coord=vert_l.children[i][i - 1], color=edge_list, rad=edge_list)
        for poly in res['polygons']:
            poly_list = point_class.PointsList(poly_l, color=poly_l, rad=poly_l)
            for vert in poly:
                point_class.Point(poly_list, coord=vert_l.children[vert], color=poly_list, rad=poly_list)

        TREE_MODEL.insertRow(TREE_MODEL.rowCount())
        atoms_index = TREE_MODEL.index(0, 0, by_point=edges_l)
        TREE_MODEL.attachObserver(atoms_index, 'Line')
        atoms_index = TREE_MODEL.index(0, 0, by_point=poly_l)
        TREE_MODEL.attachObserver(atoms_index, 'Plane')

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()

def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Voronoi')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions