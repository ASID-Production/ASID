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
                    polygons[pi] = p
                parea += p.area
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
    from PySide6.QtWidgets import QCheckBox, QFileDialog
    from ..ChemPack.ui.select_mol_dialog import SelectMolDialog
    from ..ChemPack.MoleculeClass import Atom
    import numpy as np
    import cpplib

    def saveCsv(sum_p, p_list):

        def symFromMat(mat):
            end = ['x', 'y', 'z', '']
            line = []
            for i in range(3):
                l = ''
                for j in range(4):
                    if mat[i][j] != 0:
                        if j == 3:
                            l += f'{mat[i][j]:+.15g}' if l else f'{mat[i][j]:.15g}'
                        else:
                            if mat[i][j] == 1:
                                l += '+' + end[j] if l else end[j]
                            elif mat[i][j] == -1:
                                l += '-' + end[j]
                            else:
                                l += f'{mat[i][j]:+.15g}' + end[j] if l else f'{mat[i][j]:.15g}' + end[j]
                line.append(l)
            return ','.join(line)

        out, _ = QFileDialog.getSaveFileName(caption='Table', filter='*.csv')
        if not out:
            return
        out = open(out, 'w')
        out.write(f'Sum polyhedra\nArea;{sum_p.area}\nVolume;{sum_p.volume}\n')
        for p in p_list:
            out.write(f'{p.atom.name};{symFromMat(p.symop)}\nArea;{p.area}\nVolume;{p.volume}\nAtom 1;Atom 2;Area;Solid angle;Solid angle fraction\n')
            s = 0
            for pol in p.polygons:
                s += pol.solid_angle
            for pol in p.polygons:
                out.write(f'{pol.atoms[0].name}({symFromMat(pol.symops[0])});{pol.atoms[1].name}({symFromMat(pol.symops[1])});{pol.area};{pol.solid_angle};{pol.solid_angle/s}\n')
        out.close()

    def saveObj(sum_p, p_list):
        from .MoleculeClass import fracToDec
        out, _ = QFileDialog.getSaveFileName(caption='Obj file', filter='*.obj')
        if not out:
            return
        ...
        out = open(out, 'w')
        v = []
        vn = []
        f = []
        iv = 1
        atoms = [p.atom for p in p_list]
        for p in sum_p.polygons:
            z = 0
            verts = []
            for vert in p.points:
                z += 1
                v.append(fracToDec(*vert.cell, [vert.coord])[0])
                verts.append(v[-1])
            ac = p.atoms[0].coord if p.atoms[0] in atoms else p.atoms[1].coords
            n = np.cross(verts[1]-verts[0], verts[-1]-verts[0])
            n = n/np.linalg.norm(n)
            s = np.sign(np.dot(n, verts[0] - ac))
            n *= s
            vn.append(n)
            f += [[(0+iv, None, len(vn)), (1+iv+x, None, len(vn)), (2+iv+x, None, len(vn))] for x in range(z-2)]
            iv += z
        for vert in v:
            out.write(f'v  {vert[0]:>12.4f}{vert[1]:>12.4f}{vert[2]:>12.4f}\n')
        out.write('\n')
        for n in vn:
            out.write(f'vn {n[0]:>12.4f}{n[1]:>12.4f}{n[2]:>12.4f}\n')
        out.write('\n')
        for face in f:
            out.write(f'f  {face[0][0]:>6d}//{face[0][2]:<6d}{face[1][0]:>6d}//{face[1][2]:<6d}{face[2][0]:>6d}//{face[2][2]:<6d}\n')
        out.write('\n')
        out.close()




    def process(mol_sys):
        opengl_widget = MAIN_WIDGET.findChild(QOpenGLWidget, "OpenGLWidget")
        sel = opengl_widget.selection_model.selection()
        sel_ind = []
        [sel_ind := sel_ind + [y.internalPointer() for y in x.indexes()] for x in sel]
        atoms = mol_sys[1].children[0].children
        sym_codes = atoms[0].cif_sym_codes
        sel = [x._atom for x in sel_ind if x._atom in atoms]
        if not sel:
            return



        if atoms[0].cif_cell_a == 1 and atoms[0].cif_cell_b == 1 and atoms[0].cif_cell_c == 1:
            coords = np.array([x.coord for x in atoms])
            a, b, c = abs(max(coords[:,0]) - min(coords[:,0])), abs(max(coords[:,1]) - min(coords[:,1])), abs(max(coords[:,2]) - min(coords[:,2]))
            a, b, c = a + 0.25, b + 0.25, c + 0.25,
            coords[:,0] /= a
            coords[:,1] /= b
            coords[:,2] /= c
            for i, at in enumerate(atoms):
                at.cif_cell_a = a
                at.cif_cell_b = b
                at.cif_cell_c = c
                at.cif_frac_coords = coords[i]
        cell = [atoms[0].cif_cell_a,
                atoms[0].cif_cell_b,
                atoms[0].cif_cell_c,
                atoms[0].cif_cell_al,
                atoms[0].cif_cell_be,
                atoms[0].cif_cell_ga]

        symms = [x[1] for x in sym_codes]
        data = [tuple([x.atom_type, *list(x.cif_frac_coords)]) for x in atoms]
        bools = [True if x in sel else False for x in atoms]

        res = cpplib.VoronoiCalculation(cell, symms, data, bools, 15.0)

        polyh = Polyhedron.fromCppLib(res, atoms, cell)

        points_list = point_class.PointsList(TREE_MODEL.getRoot(), name='Voronoi polyhedron')
        colors_list = point_class.PointsList(points_list, name='Colors')
        atoms_poly_list = point_class.PointsList(points_list, name='Atoms')
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
        if CHECKBOX.isChecked():
            saveCsv(sum_poly, polyh)
        if CHECKBOX_OBJ.isChecked():
            saveObj(sum_poly, polyh)
        ls = sum_poly.createList(points_list, colors, new=True)
        ls.name = 'Polyhedra'

        TREE_MODEL.insertRow(TREE_MODEL.rowCount())
        for l in lists:
            index = TREE_MODEL.index(0, 0, by_point=ls.children[1])
            TREE_MODEL.attachObserver(index, 'Line')
            index = TREE_MODEL.index(0, 0, by_point=ls.children[2])
            TREE_MODEL.attachObserver(index, 'Plane')

    global DIALOG
    global CHECKBOX
    global CHECKBOX_OBJ
    CHECKBOX = QCheckBox('Save .csv')
    CHECKBOX_OBJ = QCheckBox('Save .obj')
    DIALOG = SelectMolDialog(MOLECULE_SYSTEMS, process)
    DIALOG.layout().insertWidget(1, CHECKBOX)
    DIALOG.layout().insertWidget(2, CHECKBOX_OBJ)
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