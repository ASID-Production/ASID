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

def triangulation(points_coords, points):
    from scipy.spatial import ConvexHull
    from itertools import combinations
    triang = ConvexHull(points_coords)
    hull = [[points[y] for y in x] for x in triang.simplices]
    edges = [list(combinations(x, 2)) for x in triang.simplices]
    edges = [y for x in edges for y in x]
    edges = map(sorted, edges)
    edges = list(set(tuple(x) for x in edges))
    edges = [[points[y] for y in x] for x in edges]
    return hull, edges


def execute():
    from . import MAIN_WIDGET, MOLECULE_SYSTEMS, TREE_MODEL
    from ... import point_class
    from PySide6.QtOpenGLWidgets import QOpenGLWidget
    import numpy as np

    opengl_widget = MAIN_WIDGET.findChild(QOpenGLWidget, "OpenGLWidget")
    sel = opengl_widget.selection_model.selection()
    sel = [x.indexes()[0].internalPointer() for x in sel]
    for p in sel:
        if p._atom and len(p._atom.bonds()) >= 4:
            points_coords = np.array([x.get(p._atom).coord for x in p._atom.bonds()])
            points = [x.get(p._atom).point() for x in p._atom.bonds()]
            hull, edges = triangulation(points_coords, points)

            points_list = point_class.PointsList(TREE_MODEL.getRoot(), name=f'Polyhedron {p.name}')
            edges_l = point_class.PointsList(points_list, rad=0.01  , color=[0,0,0,1], name='Edges')
            poly_l = point_class.PointsList(points_list, color=[0.0,1.0,0.0,0.5], name='Polygons')
            for edge in edges:
                edge_list = point_class.PointsList(edges_l, color=edges_l, rad=edges_l)
                for vert in edge:
                    p = point_class.Point(edge_list, coord=vert, color=edge_list, name=vert, rad=edge_list)
            for poly in hull:
                poly_list = point_class.PointsList(poly_l, color=poly_l, rad=poly_l)
                for vert in poly:
                    p = point_class.Point(poly_list, coord=vert, color=poly_list, name=vert, rad=poly_list)

            TREE_MODEL.insertRow(TREE_MODEL.rowCount())
            atoms_index = TREE_MODEL.index(0, 0, by_point=edges_l)
            TREE_MODEL.attachObserver(atoms_index, 'Line')
            atoms_index = TREE_MODEL.index(0, 0, by_point=poly_l)
            TREE_MODEL.attachObserver(atoms_index, 'Plane')
