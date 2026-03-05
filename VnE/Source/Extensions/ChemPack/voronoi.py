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

import numpy as np


class Polyhedron:
    def __init__(self, polygons, edges, points, center, atom, symop, volume, area):
        for p in points:
            p.polyhedra.append(self)
        self.polygons = polygons
        self.edges = edges
        self.points = points
        self.center = center
        self.atom = atom
        self.symop = symop
        self.volume = volume
        self.area = area
        self.llist = None

    @classmethod
    def mergePolyhedra(cls, polyh):

        polygons = []
        edges = []
        points = []
        center = []
        atom = None
        symop = (0,0,0)
        volume = 0
        area = 0

        poly = []
        poly_dict = {}
        for i, p1 in enumerate(polyh):
            edges += p1.edges
            points += p1.points
            center.append(p1.center)
            atom = p1.atom
            volume += p1.volume
            for p in p1.polygons:
                area += p.area
                symops = [tuple([tuple([round(v, 5) for v in x]) for x in s]) for s in p.symops]
                poly1 = frozenset([(p.atoms[0], symops[0]), (p.atoms[1], symops[1])])
                for poly2 in poly:
                    if poly1 == poly2:
                        area -= 2*p.area
                        polygons.remove(poly_dict[poly1])
                        break
                else:
                    polygons.append(p)
                    poly_dict[poly1] = p
                    poly.append(poly1)
        '''for p in points:
            col = p.elliminateCollision(points)
            for rp in col:
                points.remove(rp)'''
        center = sum(center)/len(center)
        return cls(polygons, edges, points, center, atom, symop, volume, area)

    @classmethod
    def fromCppLib(cls, d, atoms, cell):
        def validate(p):
            if not p['edges'] and not p['polygons'] and not p['vertices']:
                return False
            return True
        ret = []
        polyhedra = {i: x for i, x in enumerate(d['voronoi_cells']['polyhedra']) if validate(x)}
        vert_list = d['voronoi_cells']['vertices']
        points = {}
        edges = {}
        polygons = {}
        for polyh in polyhedra:
            parea = 0
            polyhi = polyh
            polyh = polyhedra[polyh]
            polyh_atom = d['voronoi_cells']['polygons'][polyh['polygons'][0]]['atoms'][0] if d['voronoi_cells']['polygons'][polyh['polygons'][0]]['atoms'][0] in d['voronoi_cells']['polygons'][polyh['polygons'][1]]['atoms'] else d['voronoi_cells']['polygons'][polyh['polygons'][0]]['atoms'][1]
            p_list = []
            e_list = []
            pol_list = []
            for pi in polyh['vertices']:
                p = points.get(pi, None)
                if not p:
                    p = Point(np.array(vert_list[pi]), cell)
                    points[pi] = p
                p_list.append(p)
            for ei in polyh['edges']:
                e = edges.get(ei, None)
                if not e:
                    pl = d['voronoi_cells']['edges'][ei]
                    e = Edge([points[pl[0]], points[pl[1]]])
                    edges[ei] = e
                e_list.append(e)
            for pi in polyh['polygons']:
                p = polygons.get(pi, None)
                if not p:
                    area = d['voronoi_cells']['polygons'][pi]['area']
                    vert = [points[pi] for pi in d['voronoi_cells']['polygons'][pi]['vertices']]
                    edg = [edges[ei] for ei in d['voronoi_cells']['polygons'][pi]['edges']]
                    at = d['voronoi_cells']['polygons'][pi]['atoms']
                    pa = at.index(polyh_atom) - 1
                    tr = [np.array(d['unit_cell'][at[0]]['shift']), np.array(d['unit_cell'][at[1]]['shift'])]
                    tr_a = np.array(d['voronoi_cells']['polygons'][pi]['shift'])
                    symmref = [d['unit_cell'][at[0]]['symmref'], d['unit_cell'][at[1]]['symmref']]
                    at = [atoms[d['unit_cell'][at[0]]['index']], atoms[d['unit_cell'][at[1]]['index']]]
                    symops = [at[0].symCodeMat(symmref[0]), at[0].symCodeMat(symmref[1])]
                    symops[0][:-1, 3] += tr[0]
                    symops[1][:-1, 3] += tr[1]
                    symops[pa][:-1, 3] += tr_a
                    solid_angle = d['voronoi_cells']['polygons'][pi]['solid_angle']
                    p = Polygon(edg, at, symops, vert, area, solid_angle)
                    parea += area
                    polygons[pi] = p
                pol_list.append(p)
            tr = np.array(polyh['center']) - atoms[d['unit_cell'][polyhi]['index']].cif_frac_coords
            symmref = d['unit_cell'][polyhi]['symmref']
            symop = atoms[d['unit_cell'][polyhi]['index']].symCodeMat(symmref)
            symop[:-1,3] += tr
            p = cls(pol_list, e_list, p_list, np.array(polyh['center']), atoms[d['unit_cell'][polyhi]['index']], symop, polyh['volume'], parea)
            ret.append(p)
        return ret

    def createList(self, parent=None, colors=(None, None, None), new=False):
        if not self.llist or new:
            from ... import point_class
            v_color = colors[0] if colors[0] else [1,0,0,1]
            e_color = colors[1] if colors[1] else [0,0,0,1]
            p_color = colors[2] if colors[2] else [0,1,0,0.5]
            self.llist = point_class.PointsList(parent, name=f'{self.atom.name}')
            self.llist.addProperty('volume', self.volume)
            self.llist.addProperty('area', self.area)
            vl = point_class.PointsList(self.llist, rad=0.1  , color=v_color, name='Vertices')
            for i, v in enumerate(self.points):
                v.createPoint(str(i), vl, new=new)
            el = point_class.PointsList(self.llist, rad=0.01, color=e_color, freq=1, hfreq=1, name='Edges')
            for e in self.edges:
                e.createList(el, new=new)
            pl = point_class.PointsList(self.llist, rad=0.01, color=p_color, name='Polygons')
            for p in self.polygons:
                p.createList(pl, new=new)
        return self.llist

    def replacePoint(self, old, new):
        if old in self.points:
            self.points.remove(old)
            self.points.append(new)

    def copyTo(self, symop):
        points = {p: p.copyTo(symop) for p in self.points}
        edges = {e: Edge([points[e.points[0]], points[e.points[1]]]) for e in self.edges}
        polygons = [Polygon([edges[e] for e in p.edges], p.atoms, [symop @ p.symops[0], symop @ p.symops[1]], [points[point] for point in p.points], p.area, p.solid_angle) for p in self.polygons]
        new_center = np.array((np.append(self.center, 1.0) @ symop.T)[:-1])
        new_symop = symop @ self.symop

        copy = Polyhedron(polygons, list(edges.values()), list(points.values()), new_center, self.atom, new_symop, self.volume, self.area)
        return copy

class Polygon:
    def __init__(self, edges, atoms, symops, points, area=0, solid_angle=0):
        for p in points:
            p.polygons.append(self)
        self.edges = edges
        self.atoms = atoms
        self.symops = symops
        self.points = points
        self.area = area
        self.solid_angle = solid_angle
        self.llist = None

    def createList(self, parent=None, new=False):
        if not self.llist or new:
            from ... import point_class
            points = [x.createPoint() for x in self.points]
            if parent:
                self.llist = point_class.PointsList(parent, rad=parent, color=parent, name=f'{self.atoms[0].name}-{self.atoms[1].name}', area=self.area)
            else:
                self.llist = point_class.PointsList(parent, rad=0.01, color=[0, 1, 0, 0.5], name=f'{self.atoms[0].name}-{self.atoms[1].name}', area=self.area, solid_angle=self.solid_angle)
            self.llist.addProperty('area', self.area)
            self.llist.addProperty('solid_angle', self.solid_angle)
            a = point_class.Point(self.llist, coord=points[0].coord, rad=self.llist, color=self.llist, name=f'{points[0].name}')
            point_class.Point(self.llist, coord=points[1].coord, rad=self.llist, color=self.llist, name=f'{points[1].name}')
            b = point_class.Point(self.llist, coord=points[2].coord, rad=self.llist, color=self.llist, name=f'{points[2].name}')
            if len(points) > 3:
                for p in points[3:]:
                    c = p
                    point_class.Point(self.llist, coord=a.coord, rad=self.llist, color=self.llist, name=f'{a.name}')
                    point_class.Point(self.llist, coord=b.coord, rad=self.llist, color=self.llist, name=f'{b.name}')
                    b = point_class.Point(self.llist, coord=c.coord, rad=self.llist, color=self.llist, name=f'{c.name}')
        return self.llist

    def replacePoint(self, old, new):
        if old in self.points:
            self.points.remove(old)
            self.points.append(new)

class Edge:
    def __init__(self, points):
        for p in points:
            p.edges.append(self)
        self.points = points
        self.llist = None

    def createList(self, parent=None, new=False):
        if not self.llist or new:
            from ... import point_class
            p1 = self.points[0].createPoint()
            p2 = self.points[1].createPoint()
            if parent:
                self.llist = point_class.PointsList(parent, rad=parent, color=parent, freq=parent, hfreq=parent, name=f'{p1.name}-{p2.name}')
            else:
                self.llist = point_class.PointsList(parent, rad=0.01, color=[0,0,0,1], freq=1, hfreq=1, name=f'{p1.name}-{p2.name}')
            point_class.Point(self.llist, coord=p1.coord, rad=self.llist, color=self.llist, freq=self.llist, hfreq=self.llist, name=f'{p1.name}')
            point_class.Point(self.llist, coord=p2.coord, rad=self.llist, color=self.llist, freq=self.llist, hfreq=self.llist, name=f'{p2.name}')
        return self.llist


    def replacePoint(self, old, new):
        if old in self.points:
            self.points.remove(old)
            self.points.append(new)

class Point:
    def __init__(self, coord, cell):
        self.coord = coord
        self.edges = []
        self.polygons = []
        self.polyhedra = []
        self.cell = cell
        self.lpoint = None

    def createPoint(self, name='0', parent=None, new=False):
        if not self.lpoint or new:
            from ... import point_class
            from .MoleculeClass import fracToDec
            coord = fracToDec(*self.cell, [self.coord])[0]
            if parent:
                self.lpoint = point_class.Point(parent, coord=coord, rad=parent, color=parent, name=f'{name}')
            else:
                self.lpoint = point_class.Point(coord=coord, rad=0.1  , color=[1,0,0,1], name=f'{name}')

        return self.lpoint


    def checkForCollision(self, points):
        col = []
        for p in points:
            if not p is self:
                if np.all(np.isclose(p.coord, self.coord, atol=0.0001)):
                    col.append(p)
        return col

    def elliminateCollision(self, points):
        col = self.checkForCollision(points)
        for p in col:
            for o in p.edges + p.polygons + p.polyhedra:
                o.replacePoint(p, self)
        return col

    def copyTo(self, symop):
        new_coord = np.array((np.append(self.coord, 1.0) @ symop.T)[:-1])
        return Point(new_coord, self.cell)

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
        n = (cos(al)-(cos(ga) * cos(be)))/sin(ga)
        p = (1-cos(al)**2-cos(be)**2-cos(ga)**2+2*cos(al)*cos(be)*cos(ga))**0.5
        mat = np.array([[a, b*cos(ga), c*cos(be)],
                        [0, b*sin(ga), c*(cos(al)-cos(be)*cos(ga))/sin(ga)],
                        [0, 0, c*p/sin(ga)]])
        for i in range(len(coords)):
            coord = np.array(coords[i])[...,np.newaxis]
            coord = mat @ coord
            coord = coord.transpose().squeeze()
            coords[i] = coord
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

        polyh = Polyhedron.fromCppLib(res, atoms, cell)

        points_list = point_class.PointsList(TREE_MODEL.getRoot(), name=f'Voronoi polyhedron')
        colors_list = point_class.PointsList(points_list, name='Colors')
        atoms_poly_list = point_class.PointsList(points_list, name='Atoms')
        sum_poly_list = point_class.PointsList(points_list, name='Polyhedra')
        v_color = point_class.Point(colors_list, name='Vertex color', color=[1, 0, 0, 1])
        e_color = point_class.Point(colors_list, name='Edge color', color=[0, 0, 0, 1])
        p_color = point_class.Point(colors_list, name='Polygon color', color=[0, 1, 0, 0.5])
        colors = [v_color, e_color, p_color]
        lists = []
        for i, p in enumerate(polyh):
            if not np.all(np.isclose(np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]), p.symop, atol=0.0001)):
                p = p.copyTo(np.linalg.inv(p.symop))
                polyh[i] = p
            lists.append(p.createList(atoms_poly_list, colors))
        sum_poly = Polyhedron.mergePolyhedra(polyh)
        ls = sum_poly.createList(sum_poly_list, colors, new=True)

        TREE_MODEL.insertRow(TREE_MODEL.rowCount())
        for l in lists:
            index = TREE_MODEL.index(0, 0, by_point=ls.children[1])
            TREE_MODEL.attachObserver(index, 'Line')
            index = TREE_MODEL.index(0, 0, by_point=ls.children[2])
            TREE_MODEL.attachObserver(index, 'Plane')
        '''edges_l = point_class.PointsList(points_list, rad=0.01  , color=[0,0,0,1], name='Edges')
        vert_l = point_class.PointsList(points_list, rad=0.1  , color=[1,0,0,1], name='Vertex')
        poly_l = point_class.PointsList(points_list, color=[0.0,1.0,0.0,0.5], name='Polygons')

        atoms_list = point_class.PointsList(TREE_MODEL.getRoot(), name=f'Atoms voronoi polyhedrons')
        a_lists = {atoms[i]: {'ind': i} for i, b in enumerate(bools) if b}
        for a in [x for x in atoms if x in sel]:
            a_lists[a]['edges'] = point_class.PointsList(atoms_list, rad=0.01, color=[0, 0, 0, 1], name='Edges')
            a_lists[a]['vertex'] = point_class.PointsList(atoms_list, rad=0.1, color=[1, 0, 0, 1], name='Vertex')
            a_lists[a]['polygon'] = point_class.PointsList(atoms_list, color=[0.0, 1.0, 0.0, 0.5], name='Polygons')

        vert_tr = [False for x in res['vertices']]
        for i, p in enumerate(res['polyhedra']):
            a = atoms[i]
            cent = np.array(p['center'])
            tr = np.round(a.cif_frac_coords - cent)
            for vi in p['vertices']:
                if not vert_tr[vi]:
                    vert_tr[vi] = True
                    #res['vertices'][vi] = tuple(tr + np.array(res['vertices'][vi]))

        res['vertexes_dec'] = fracToDec(*cell, [list(x) for x in res['vertices']])

        sum_area = 0

        for i, v in enumerate(res['vertexes_dec']):
            point_class.Point(vert_l, coord=v, color=vert_l, rad=vert_l, name=str(i))
        for poly in res['polygons']:
            #sum_area += poly['area']
            #poly_list = point_class.PointsList(poly_l, color=poly_l, rad=poly_l, name=f"{atoms[poly['atoms'][0]].name}-{atoms[poly['atoms'][1]].name}")
            poly_list = point_class.PointsList(poly_l, color=poly_l, rad=poly_l)
            poly = poly['vertices']
            for i, vert in enumerate(poly):
                edge_list = point_class.PointsList(edges_l, color=edges_l, rad=edges_l)
                point_class.Point(edge_list, coord=vert_l.children[vert], color=edge_list, rad=edge_list)
                point_class.Point(edge_list, coord=vert_l.children[poly[i-1]], color=edge_list, rad=edge_list)
            point_class.Point(poly_list, coord=vert_l.children[poly[0]], color=poly_list, rad=poly_list, name=poly[0])
            point_class.Point(poly_list, coord=vert_l.children[poly[1]], color=poly_list, rad=poly_list, name=poly[1])
            point_class.Point(poly_list, coord=vert_l.children[poly[2]], color=poly_list, rad=poly_list, name=poly[2])
            a, b, c = poly[0], poly[1], poly[2]
            for i in range(3, len(poly)):
                point_class.Point(poly_list, coord=vert_l.children[a], color=poly_list, rad=poly_list, name=str(a))
                point_class.Point(poly_list, coord=vert_l.children[c], color=poly_list, rad=poly_list, name=str(c))
                point_class.Point(poly_list, coord=vert_l.children[poly[i]], color=poly_list, rad=poly_list, name=str(poly[i]))
                a, b, c = a, c, poly[i]'''

    global DIALOG
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.show()

def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('PVD')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions