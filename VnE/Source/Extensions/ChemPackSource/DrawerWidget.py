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

from abc import ABC, abstractmethod
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import *
from PySide6.QtUiTools import loadUiType
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from ..ChemPack import PALETTE
from OpenGL.GL import *
import numpy as np
import time
from .contacts import Node, Pack, Condition, dist, angle, avgDiff, maxMeanPlaneDiff
from .MoleculeClass import Atom
from ... import point_class
from ...Facade import RenderFacade
from ...Scenes import Scene
from ...UniformBuffers import SceneUniformBuffer
from ... import Observers

SINGLE_OBSERVER = None
DRAW_WIDGET = None


class Drawing:

    def __init__(self, points: point_class.PointsList):
        self.points = points
        self.connections = {}
        self.contacts = {}
        self.point_node = {}
        self.node_point = {}
        self.pack = Pack()
        self.conditions_d = {'contacts': [],
                             'angle': [],
                             'avgDiff': [],
                             'maxMeanPlaneDiff': []}
        self.conditions_value_d = {'contacts': [],
                                   'angle': [],
                                   'avgDiff': [],
                                   'maxMeanPlaneDiff': []}

    def add_point(self, point):
        if point not in self.connections:
            self.connections[point] = {}
            self.contacts[point] = {}
            self.point_node[point] = Node(self.pack, str_atom=Atom([0,0,0], atom_type=point.atom_type, name=point.label))
            self.node_point[self.point_node[point]] = point

    def changeAtomType(self, point):
        self.point_node[point].setStructAtom(Atom([0,0,0], atom_type=point.atom_type, name=point.label))

    def labelContact(self, ind):
        try:
            dist = self.conditions_value_d['contacts'][ind]
            node1, node2 = self.conditions_d["contacts"][ind].nodes
            point1, point2 = self.node_point[node1], self.node_point[node2]
        except IndexError:
            return [None, None]
        return [f'{point1.label}:{point2.label}', str(dist)]

    def labelAngle(self, ind):
        try:
            angle = self.conditions_value_d['angle'][ind]
            cond = self.conditions_d["angle"][ind]
        except IndexError:
            return [None, None]
        if cond.split == -1:
            node1, node2, node3 = cond.nodes
            point1, point2, point3 = self.node_point[node1], self.node_point[node2], self.node_point[node3]
            return [f'{point1.label}-{point2.label}-{point3.label}', str(angle)]
        else:
            points = [self.node_point[x] for x in cond.nodes]
            return ['-'.join([x.label for x in points[:cond.split]]) + ':' + '-'.join([x.label for x in points[cond.split:]]), str(angle)]

    def labelAvgDiff(self, ind):
        try:
            dist = self.conditions_value_d['avgDiff'][ind]
            cond = self.conditions_d["avgDiff"][ind]
            nodes = cond.nodes
            split = cond.split
            points = [self.node_point[x] for x in nodes]
        except IndexError:
            return [None, None]
        return ['-'.join([x.label for x in points[:split]]) + ':' + '-'.join([x.label for x in points[split:]]), str(dist)]

    def labelMaxMeanPlaneDiff(self, ind):
        try:
            dist = self.conditions_value_d['maxMeanPlaneDiff'][ind]
            cond = self.conditions_d["maxMeanPlaneDiff"][ind]
            nodes = cond.nodes
            points = [self.node_point[x] for x in nodes]
        except IndexError:
            return [None, None]
        return ['-'.join([x.label for x in points]), str(dist)]

    def setDataContact(self, ind, dist):
        try:
            self.conditions_d['contacts'][ind].setCompFunc(lambda x: x <= dist)
            self.conditions_value_d['contacts'][ind] = dist
        except IndexError:
            return

    def setDataAngle(self, ind, angle):
        try:
            self.conditions_value_d['angle'][ind] = angle
            self.conditions_d['angle'][ind].setCompFunc(lambda x: x >= angle or x <= 180-angle)
        except IndexError:
            return

    def setDataAvgDiff(self, ind, dist):
        try:
            self.conditions_d['avgDiff'][ind].setCompFunc(lambda x: x <= dist)
            self.conditions_value_d['avgDiff'][ind] = dist
        except IndexError:
            return

    def setDataMaxMeanPlaneDiff(self, ind, dist):
        try:
            self.conditions_d['maxMeanPlaneDiff'][ind].setCompFunc(lambda x: x <= dist)
            self.conditions_value_d['maxMeanPlaneDiff'][ind] = dist
        except IndexError:
            return

    def add_connection(self, point1, point2):
        point1, linep1 = point1
        point2, linep2 = point2

        if point1 in self.connections:
            if point2 not in self.connections[point1]:
                self.connections[point1][point2] = linep2
        else:
            self.add_point(point1)
        if point2 in self.connections:
            if point1 not in self.connections[point2]:
                self.connections[point2][point1] = linep1
        else:
            self.add_point(point2)
        self.point_node[point1].addConnect(self.point_node[point2])

    def add_contact(self, point1, point2):
        point1, linep1 = point1
        point2, linep2 = point2

        self.conditions_d['contacts'].append(Condition(lambda x: x <= 0.0, dist, [self.point_node[point1], self.point_node[point2]]))
        self.conditions_value_d['contacts'].append(0.0)

        if point1 in self.contacts:
            if point2 not in self.contacts[point1]:
                self.contacts[point1][point2] = linep2
        else:
            self.add_point(point1)
        if point2 in self.contacts:
            if point1 not in self.contacts[point2]:
                self.contacts[point2][point1] = linep1
        else:
            self.add_point(point2)

    def addAngle(self, *points, split=-1):
        self.conditions_value_d['angle'].append(0.0)
        self.conditions_d['angle'].append(Condition(lambda x: x >= 0.0 or x <= 0.0, lambda *x: angle(x, rad=False, split=split), [self.point_node[x] for x in points]))
        self.conditions_d['angle'][-1].split = split

    def addAvgDiff(self, *points, split=1):
        self.conditions_value_d['avgDiff'].append(0.0)
        self.conditions_d['avgDiff'].append(
            Condition(lambda x: x <= 0.0, lambda *x: avgDiff(x, split=split), [self.point_node[x] for x in points]))
        self.conditions_d['avgDiff'][-1].split = split

    def addMaxMeanPlaneDiff(self, *points):
        self.conditions_value_d['maxMeanPlaneDiff'].append(0.0)
        self.conditions_d['maxMeanPlaneDiff'].append(
            Condition(lambda x: x <= 0.0, lambda *x: maxMeanPlaneDiff(*x), [self.point_node[x] for x in points]))

    def removeAngle(self, ind):
        if ind == -1:
            return
        cond = self.conditions_d['angle'].pop(ind)

    def removeContact(self, ind):
        if ind == -1:
            return
        cond = self.conditions_d['contacts'].pop(ind)
        point1, point2 = cond.nodes
        self.conditions_value_d['contacts'].pop(ind)
        self.contacts[point1][point2].parent.destroy()
        self.contacts[point1].pop(point2)
        self.contacts[point2].pop(point1)

    def removeAvgDiff(self, ind):
        if ind == -1:
            return
        cond = self.conditions_d['avgDiff'].pop(ind)

    def removeMaxMeanPlaneDiff(self, ind):
        if ind == -1:
            return
        cond = self.conditions_d['maxMeanPlaneDiff'].pop(ind)
        self.conditions_value_d['maxMeanPlaneDiff'].pop(ind)

    def remove_point(self, point):
        for point2 in self.connections[point].copy():
            line = self.connections[point].pop(point2)
            self.connections[point2].pop(point)
            line.parent.destroy()
        for point2 in self.contacts[point].copy():
            line = self.contacts[point].pop(point2)
            self.contacts[point2].pop(point)
            line.parent.destroy()
        self.connections.pop(point)
        self.contacts.pop(point)
        node = self.point_node.pop(point)
        self.node_point.pop(node)
        node.delete()
        point.destroy()


class SolutionOrganizer:

    def __init__(self, *conditions):
        self.conditions = list(conditions)
        self.separate = {}
        self.separateConditions()
        self.headers = {key: '' for key in self.separate}
        self.createHeader()
        self.tables = []

    def createHeader(self):
        if not self.conditions:
            return
        for packs in self.separate:
            header_packs = ';'.join(['-'.join([y.struct_atom.name if 'name' in y.struct_atom.__dict__ else str(None) for y in x.nodes]) for x in packs])
            conditions = self.separate[packs]
            header_conditions = ';'.join(['--'.join([y.struct_atom.name if 'name' in y.struct_atom.__dict__ else str(None) for y in x.nodes]).__add__(';Value') for x in conditions])
            header = f'{header_packs};{header_conditions}'
            self.headers[packs] = header

    def fillTable(self):
        if not self.conditions:
            return
        for packs in self.separate:
            solutions = self.conditions[0].intersection(*self.separate[packs])
            lines = []
            for i in range(len(solutions[0][packs[0]])):
                solution = {}
                for pack in packs:
                    for cond in solutions:
                        if pack in cond:
                            solution[pack] = cond[pack]
                            break
                packs_line = ';'.join(['--'.join([y.real_atoms[solution[y.pack][i]].point().name for y in x.nodes]) for x in packs])
                conditions_value = [str(round(x.value_func(*[y.real_atoms[solution[y.pack][i]] for y in x.nodes]), 3)) for x in self.separate[packs]]
                conditions_line = ['--'.join([y.real_atoms[solution[y.pack][i]].point().name for y in x.nodes]) for x in self.separate[packs]]
                zip_line = [f'{conditions_line[i]};{conditions_value[i]}' for i in range(len(conditions_line))]
                zip_line = ';'.join(zip_line)
                lines.append(f'{packs_line};{zip_line}')
            lines.insert(0, self.headers[packs])
            self.tables.append('\n'.join(lines))

    def getTable(self):
        return '\n\n'.join(self.tables)

    def separateConditions(self):
        connections = {}
        if not self.conditions:
            return
        if len(self.conditions) == 1:
            self.separate[tuple(self.conditions[0].packs.keys())] = self.conditions
        for condition in self.conditions:
            for pack in condition.packs:
                if pack in connections:
                    if condition not in connections[pack]:
                        connections[pack].append(condition)
                else:
                    connections[pack] = [condition]

        rec_connection = {key: x.copy() for key, x in connections.items()}
        rec_conds = self.conditions.copy()

        def rec(cond, res=None, packs=None):
            if packs is None:
                packs = []
            if res is None:
                res = []
            if cond in res:
                return res, packs
            else:
                res.append(cond)
                rec_conds.remove(cond)
            for pack in cond.packs:
                if pack not in packs:
                    packs.append(pack)
                for cond2 in rec_connection[pack]:
                    res, packs = rec(cond2, res, packs)
            return res, packs

        while rec_conds:
            value, key = rec(rec_conds[0])
            self.separate[tuple(key)] = value


class SimpleDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self, parent=None):
        super().__init__(parent=parent)


class TableModel(QAbstractTableModel):
    def __init__(self, parent=None, label=None, data=None, setter=None):
        self.label = label
        self.data_label = data
        self.setter = setter
        self.selected = []
        self.selected_len = 0
        self.rows = 0
        self.columns = 2
        self.props = {}
        super().__init__(parent)

    def rowCount(self, parent=QModelIndex(), *args, **kwargs) -> int:
        return self.rows

    def columnCount(self, parent=QModelIndex(), *args, **kwargs):
        return self.columns

    def removeRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginRemoveRows(parent, row, row + count - 1)
        self.rows -= count
        self.endRemoveRows()
        return True

    def insertRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginInsertRows(parent, row, row + count - 1)
        self.rows += count
        self.endInsertRows()
        return True

    def insertColumns(self, column, count, parent=QModelIndex(), *args, **kwargs):
        self.beginInsertColumns(parent, column, column + count - 1)
        self.columns += count
        self.endInsertColumns()
        return True

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = ...):
        if role != Qt.ItemDataRole.DisplayRole and role != Qt.ItemDataRole.ToolTipRole:
            return None
        if section == 0 and orientation == Qt.Horizontal:
            return 'Atoms'
        elif section == 1 and orientation == Qt.Horizontal:
            return 'Distance'

    def data(self, index: QModelIndex, role: int = ...):
        if role == Qt.ItemDataRole.DisplayRole:
            if index.column() == 0:
                if self.label is not None:
                    label = self.label(index.row())
                else:
                    return None
            if index.column() != 0:
                if self.data_label is not None:
                    label = self.data_label(index.row())
                else:
                    return None
            return label

    def setData(self, index, value, role=None) -> bool:
        if role == Qt.ItemDataRole.EditRole:
            try:
                value = float(value)
            except ValueError:
                return False
            self.setter(index.row(), value)
            return True
        elif role == 99:
            property = self.data(index, role=99)[0]
            for i, selected in enumerate(self.selected):
                if property in selected.internalPointer().getProperties().keys():
                    pass
                else:
                    selected.internalPointer().addProperty(property, None)
                selected.internalPointer().__setattr__(property, value)

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        if index.column() != 1:
            return Qt.ItemFlag.ItemIsEnabled
        if index.column() == 1:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable


class TableView(QtWidgets.QTableView):

    def __init__(self, remove, parent=None, model=None):
        super().__init__(parent)
        self._remove = remove
        self.setItemDelegate(SimpleDelegate(self))
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)
        self.horizontalHeader().hide()
        self.verticalHeader().hide()
        if model is not None:
            self.setModel(model)
            self.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            self.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)

    def mousePressEvent(self, e: QtGui.QMouseEvent) -> None:
        super().mousePressEvent(e)

    def remove(self):
        self._remove(self.currentIndex().row())
        self.model().removeRows(self.currentIndex().row(), 1)

    def showContextMenu(self, pos: QPoint):
        menu = QtWidgets.QMenu('Context Menu', self)
        action2 = QtGui.QAction('Delete', self)

        action2.triggered.connect(self.remove)
        menu.addAction(action2)
        menu.exec(self.mapToGlobal(pos))
        return


class UniqueId:

    def __init__(self):
        self.nums = []
        self.max = 0

    def get(self):
        if self.nums:
            return self.nums.pop(0)
        else:
            self.max += 1
            return self.max

    def remove(self, id):
        self.nums.append(id)
        self.nums.sort()


class aEvent(ABC):

    class Timer:

        def __init__(self):
            self.start_time = None
            self.check_l = []
            self.last_time = None
            self.run_b = False

        def start(self):
            self.check_l = []
            self.start_time = time.perf_counter()
            self.check_l.append(self.start_time)
            self.run_b = True
            return 0

        def stop(self):
            if self.start_time is None:
                self.start()
                self.check_l.append(self.start_time)
            else:
                self.check_l.append(time.perf_counter())
            self.start_time = None
            self.last_time = self.check_l[-1] - self.check_l[-2]
            self.run_b = False
            return self.last_time

        def check(self):
            ret_time = self.stop()
            self.start()
            return ret_time

        def time(self):
            if self.last_time is None:
                return 0
            else:
                return self.last_time

        def run(self):
            return self.run_b

    def __init__(self, *args, exclusive=False, exclusive_group=0, **kwargs):
        self.exclusive = exclusive
        self.exclusive_group = exclusive_group
        self.widget = None

    @abstractmethod
    def attach(self, widget):
        pass

    @abstractmethod
    def detach(self):
        pass

    @abstractmethod
    def assertEvent(self, event: QtCore.QEvent, widget):
        pass


class aDrawWidgetEvent(aEvent, ABC):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.timer = self.Timer()
        self.tol = 0.025

    def attach(self, widget):
        if widget is self.widget:
            return
        elif self.widget is not None:
            self.detach()
        if self.exclusive:
            for event in widget.events:
                if event.exclusive and event.exclusive_group == self.exclusive_group:
                    event.detach()
        widget.events.append(self)
        self.widget = widget

    def detach(self):
        if self.widget is not None:
            try:
                self.widget.events.remove(self)
            except ValueError:
                pass
            finally:
                self.widget = None


class Drag(aDrawWidgetEvent):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.drag_point = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.drag_point = widget.select(event.localPos())
            self.timer.start()

        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.LeftButton:
            if self.timer.run():
                self.timer.stop()
            if event.modifiers() == QtCore.Qt.NoModifier and self.timer.time() > self.tol:
                if self.drag_point is not None:
                    self.drag_point.coord = np.array([((event.localPos().x() / self.widget.width()) * 2 - 1),
                                                      ((event.localPos().y() / self.widget.height()) * (-2) + 1) / (
                                                                  self.widget.width() / self.widget.height()),
                                                      0], dtype=np.float32)

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            self.drag_point = None
            self.timer.stop()


class HighLight(aDrawWidgetEvent):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.point = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseMove:
            point = self.widget.select(event.localPos())
            if point is not self.point and point is not None:
                point.pick = 1.0
            if self.point is not None and self.point is not point:
                self.point.pick = 0.0
            self.point = point


class Draw(aDrawWidgetEvent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.atom_type = 'C'
        self.atom_type_label = 'C'
        self.atom_color = PALETTE.getColor('C')

        self.point1 = None
        self.point2 = None
        self.line1 = None
        self.line2 = None
        self.line_list = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.timer.start()
            self.point1 = self.widget.select(event.localPos())

        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.LeftButton:
            if self.timer.run():
                self.timer.stop()
            if self.timer.time() > self.tol:
                pos = event.localPos()
                if self.point1 is None:
                    id = self.widget.ids.get()
                    self.point1 = point_class.Point(parent=self.widget.p_list,
                                                    coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()),
                                                                    0], dtype=np.float32),
                                                    rad=self.widget.p_list,
                                                    label=self.atom_type_label + str(id),
                                                    id=id,
                                                    atom_type=self.atom_type,
                                                    color=self.atom_color)
                    self.widget.drawing.add_point(self.point1)
                if self.point2 is None:
                    id = self.widget.ids.get()
                    self.point2 = point_class.Point(parent=self.widget.p_list,
                                                    coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()),
                                                                    0], dtype=np.float32),
                                                    rad=self.widget.p_list,
                                                    color=self.atom_color,
                                                    atom_type=self.atom_type,
                                                    label=self.atom_type_label + str(id),
                                                    id=id)
                    self.widget.drawing.add_point(self.point2)

                if self.line_list is None:
                    self.line_list = point_class.PointsList(parent=self.widget.l_list)
                    self.line1 = point_class.Point(parent=self.line_list, coord=self.point1, color=self.point1,
                                                   rad=self.point1)
                    self.line2 = point_class.Point(parent=self.line_list, coord=self.point2, color=self.point2,
                                                   rad=self.point2)
                    self.widget.drawing.add_connection((self.point1, self.line1), (self.point2, self.line2))
                else:
                    point2 = self.widget.select(pos)
                    if point2 is None or point2 is self.point2:
                        self.point2.coord = np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                      ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                  self.widget.width() / self.widget.height()),
                                                      0], dtype=np.float32)
                    else:
                        self.point2.coord = point2

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            if self.timer.run():
                self.timer.stop()
            if self.timer.time() > self.tol:
                point2 = self.widget.select(event.localPos())
                if point2 is not None and (point2 is not self.point2 and point2 is not self.point1):
                    self.widget.ids.remove(self.point2.id)
                    self.widget.drawing.remove_point(self.point2)
                    self.point2 = point2
                    self.line_list = point_class.PointsList(parent=self.widget.l_list)
                    self.line1 = point_class.Point(parent=self.line_list, coord=self.point1, color=self.point1,
                                                   rad=self.point1)
                    self.line2 = point_class.Point(parent=self.line_list, coord=self.point2, color=self.point2,
                                                   rad=self.point2)
                    self.widget.drawing.add_connection((self.point1, self.line1), (self.point2, self.line2))

                elif self.point2 is None:
                    pos = event.localPos()
                    point_a = self.widget.select(pos)
                    if point_a is None:
                        id = self.widget.ids.get()
                        point_a = point_class.Point(parent=self.widget.p_list,
                                                    coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()), 0],
                                                                   dtype=np.float32),
                                                    color=self.atom_color,
                                                    rad=self.widget.p_list,
                                                    id=id,
                                                    atom_type=self.atom_type,
                                                    label=self.atom_type_label + str(id))
                        self.widget.drawing.add_point(point_a)
                    else:
                        point_a.label = self.atom_type_label + str(point_a.id)
                        point_a.color = self.atom_color
                        point_a.atom_type = self.atom_type
                        self.widget.drawing.changeAtomType(point_a)

                self.point1 = None
                self.point2 = None
                self.line1 = None
                self.line2 = None
                self.line_list = None

    def setType(self, atom_type):
        if type(atom_type) is list:
            self.atom_color = np.array([190, 85, 255, 255], dtype=np.float32)/255
            self.atom_type = atom_type
            self.atom_type_label = 'X'
        elif type(atom_type) is str:
            self.atom_type = atom_type
            self.atom_type_label = atom_type
            self.atom_color = np.array(PALETTE.getColor(atom_type) + [255], dtype=np.float32)/255


class Clear(aDrawWidgetEvent):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.point_d = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.timer.start()
            self.point_d = self.widget.select(event.localPos())

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            if self.timer.run():
                self.timer.stop()
            if self.timer.time() > self.tol:
                if self.point_d is not None:
                    self.widget.ids.remove(self.point_d.id)
                    self.widget.drawing.remove_point(self.point_d)
            else:
                self.point_d = None
            self.point_d = None


class Contact(aDrawWidgetEvent):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.old_color = None
        self.pc1 = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.LeftButton:
                point = self.widget.select(event.localPos())
                if point is not None:
                    if self.pc1 is None:
                        self.pc1 = point
                        self.old_color = point.color
                        point.color = np.array([0, 1, 0, 1], dtype=np.float32)
                    elif self.pc1 is not None and point is not self.pc1:
                        dl_list = point_class.PointsList(parent=self.widget.dl_list)
                        dlp1 = point_class.Point(parent=dl_list, color=self.pc1, coord=self.pc1, rad=self.pc1)
                        dlp2 = point_class.Point(parent=dl_list, color=point, coord=point, rad=point)
                        self.widget.drawing.add_contact((self.pc1, dlp1), (point, dlp2))
                        self.pc1.color = self.old_color
                        self.widget.newContact.emit()
                        self.pc1 = None

            if event.buttons() == QtCore.Qt.RightButton:
                if self.pc1 is not None:
                    self.pc1.color = self.old_color
                    self.pc1 = None
                    self.old_color = None


class Angle(aDrawWidgetEvent):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        '''self.pa1 = None
        self.pa2 = None
        self.old_colors = []'''

        self.p1 = []
        self.p2 = []
        self.old_colors = [[], []]
        self.add = self.p1, self.old_colors[0]

    def assertEvent(self, event: QtCore.QEvent, widget):
        '''if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.LeftButton:
                point = self.widget.select(event.localPos())
                if point is not None:
                    if self.pa1 is None:
                        self.pa1 = point
                        self.old_colors.append(self.pa1.color)
                        point.color = np.array([0, 1, 0, 1], dtype=np.float32)
                    elif self.pa2 is None and point is not self.pa1:
                        self.pa2 = point
                        self.old_colors.append(self.pa2.color)
                        point.color = np.array([0, 1, 0, 1], dtype=np.float32)
                    elif point is not self.pa1 and point is not self.pa2:
                        self.widget.drawing.addAngle(self.pa1, self.pa2, point)
                        self.pa1.color = self.old_colors[0]
                        self.pa2.color = self.old_colors[1]
                        self.widget.newAngle.emit()
                        self.pa1 = None
                        self.pa2 = None

            if event.buttons() == QtCore.Qt.RightButton:
                if self.pa1 is not None:
                    self.pa1.color = self.old_colors[0]
                    self.pa1 = None
                if self.pa2 is not None:
                    self.pa2.color = self.old_colors[1]
                    self.pa2 = None
                self.old_colors = []'''

        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.LeftButton:
                point = self.widget.select(event.localPos())
                if point is not None and point not in self.add[0]:
                    self.add[0].append(point)
                    self.add[1].append(point.color)
                    point.color = np.array([0, 1, 0, 1], dtype=np.float32)

            if event.buttons() == QtCore.Qt.RightButton:
                if self.add[0]:
                    point = self.add[0].pop()
                    color = self.add[1].pop()
                    point.color = color
                elif self.add[0] is self.p2:
                    self.add = self.p1, self.old_colors[0]
                    for point in self.add[0]:
                        point.color = np.array([0, 1, 0, 1], dtype=np.float32)

        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
                if self.add[0] is self.p2 and len(self.add[0]) == 0:
                    for point, color in zip(self.p1, self.old_colors[0]):
                        point.color = color
                    points = self.p1
                    self.widget.drawing.addAngle(*points, split=-1)
                    self.widget.newAngle.emit()
                    self.p1 = []
                    self.p2 = []
                    self.old_colors = [[], []]
                    self.add = self.p1, self.old_colors[0]
                    return
                if self.add[0] and len(self.add[0]) >= 3:
                    for point in self.add[0]:
                        point.color = np.array([1, 1, 0, 1], dtype=np.float32)
                if self.add[0] is self.p1 and len(self.add[0]) >= 3:
                    self.add = self.p2, self.old_colors[1]
                elif self.add[0] is not self.p1 and len(self.add[0]) >= 3:
                    for point, color in zip(self.p2, self.old_colors[1]):
                        point.color = color
                    for point, color in zip(self.p1, self.old_colors[0]):
                        point.color = color
                    points = self.p1 + self.p2
                    split = len(self.p1)
                    self.widget.drawing.addAngle(*points, split=split)
                    self.widget.newAngle.emit()
                    self.p1 = []
                    self.p2 = []
                    self.old_colors = [[], []]
                    self.add = self.p1, self.old_colors[0]


class AvgDiff(aDrawWidgetEvent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.p1 = []
        self.p2 = []
        self.old_colors = [[], []]
        self.add = self.p1, self.old_colors[0]

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.LeftButton:
                point = self.widget.select(event.localPos())
                if point is not None and point not in self.add[0]:
                    self.add[0].append(point)
                    self.add[1].append(point.color)
                    point.color = np.array([0, 1, 0, 1], dtype=np.float32)

            if event.buttons() == QtCore.Qt.RightButton:
                if self.add[0]:
                    point = self.add[0].pop()
                    color = self.add[1].pop()
                    point.color = color
                elif self.add[0] is self.p2:
                    self.add = self.p1, self.old_colors[0]
                    for point in self.add[0]:
                        point.color = np.array([0, 1, 0, 1], dtype=np.float32)
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
                if self.add[0]:
                    for point in self.add[0]:
                        point.color = np.array([1, 1, 0, 1], dtype=np.float32)
                if self.add[0] is self.p1:
                    self.add = self.p2, self.old_colors[1]
                else:
                    for point, color in zip(self.p2, self.old_colors[1]):
                        point.color = color
                    for point, color in zip(self.p1, self.old_colors[0]):
                        point.color = color
                    points = self.p1 + self.p2
                    split = len(self.p1)
                    self.widget.drawing.addAvgDiff(*points, split=split)
                    self.widget.newAvgDiff.emit()
                    self.p1 = []
                    self.p2 = []
                    self.old_colors = [[], []]
                    self.add = self.p1, self.old_colors[0]


class MaxMeanPlaneDiff(aDrawWidgetEvent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.selected_points = []
        self.old_colors = []

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.LeftButton:
                point = self.widget.select(event.localPos())
                if point is not None and point not in self.selected_points:
                    self.selected_points.append(point)
                    self.old_colors.append(point.color)
                    point.color = np.array([0, 1, 0, 1], dtype=np.float32)

            if event.buttons() == QtCore.Qt.RightButton:
                if self.selected_points:
                    point = self.selected_points.pop()
                    point.color = self.old_colors.pop()

        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
                for point, color in zip(self.selected_points, self.old_colors):
                    point.color = color

                points = self.selected_points
                self.widget.drawing.addMaxMeanPlaneDiff(*points)
                self.widget.newMaxMeanPlaneDiff.emit()

                self.selected_points = []
                self.old_colors = []


class DrawerGL(QOpenGLWidget):

    newContact = Signal(name='newContact')
    newAngle = Signal(name='newAngle')
    newAvgDiff = Signal(name='avgDiff')
    newMaxMeanPlaneDiff = Signal(name='newMaxMeanPlaneDiff')

    class AtomType:
        def __init__(self, atom_type):
            self.atom_type = atom_type

        def setAtomType(self, atom_type):
            self.atom_type = atom_type

        def __str__(self):
            if type(self.atom_type) is list:
                return 'X'
            else:
                return str(self.atom_type)

    def __init__(self, *args, parent=None):
        super().__init__(parent)
        self.surface_format = QtGui.QSurfaceFormat()
        self.surface_format.setSamples(8)
        self.setFormat(self.surface_format)
        self.error = []
        self._events = []
        self.installEventFilter(self)
        self.ready = False
        self.pressed = False
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

        self.checked_func = 'drag'
        self._checked_mod = self.AtomType('C')

        self.points = point_class.PointsList()
        self.p_list = point_class.PointsList(parent=self.points, rad=0.02, color=np.array([0.5, 0.5, 0.5, 1.0], dtype=np.float32))
        self.l_list = point_class.PointsList(parent=self.points)
        self.dl_list = point_class.PointsList(parent=self.points)
        self.drawing = Drawing(self.points)

        self.facade = None
        self.initializeGL()

        self.ids = UniqueId()

    def paintGL(self):
        if self.ready:
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            self.facade.drawScene(self.scene)
        else:
            self.ready = True
            self.initializeGL()
        return

    def resizeGL(self, w, h):
        if self.ready:
            self.makeCurrent()
            self.uniform_buffer.wh = [self.width(), self.height()]

    def select(self, pos):
        mod = self.uniform_buffer.aspect_ratio
        tol = np.array([1, 1])
        pos = [(pos.x() / self.width()) * 2 - 1, (pos.y() / self.height()) * (-2) + 1]
        selected = self.p_list.select(pos, tol=tol, mod=mod)
        if selected is not None:
            selected = selected[0]
            #print(selected.seq)
        else:
            #print(None)
            pass
        return selected

    def initializeGL(self):
        if self.ready:
            self.setMouseTracking(True)

            self.facade = RenderFacade(self)
            self.scene = self.facade.addScene(Scene)

            global SINGLE_OBSERVER
            SINGLE_OBSERVER = Observers.SingleObserver(self.facade, self.scene)

            self.label_observer = SINGLE_OBSERVER.getObserver(Observers.LabelObserver)
            self.facade.changePipelineUniforms(self.label_observer.getPipeline(), 'const_scale', 0.25)

            self.p_list.attach(SINGLE_OBSERVER.getObserver(Observers.SphereObserver))
            self.p_list.attach(SINGLE_OBSERVER.getObserver(Observers.LabelObserver))
            self.l_list.attach(SINGLE_OBSERVER.getObserver(Observers.LineObserver))
            self.dl_list.attach(SINGLE_OBSERVER.getObserver(Observers.DashedLineObserver))
            #self.l_list.attach(SINGLE_OBSERVER.getObserver(Observers.BondsObserver))
            #self.points.attach(SINGLE_OBSERVER.getObserver(Observers.LineObserver))

            self.uniform_buffer_id = self.facade.addUniformBufferToScene(self.scene, uniform_buffer_cls=SceneUniformBuffer)
            self.uniform_buffer = self.facade.getInst(self.uniform_buffer_id)
            self.uniform_buffer: SceneUniformBuffer
            self.uniform_buffer.perspective = np.array([[1,0,0,0],
                                                        [0,1,0,0],
                                                        [0,0,1,0],
                                                        [0,0,0,1]], dtype=np.float32)
            self.uniform_buffer.wh = [self.width(), self.height()]
            glClearColor(1.0, 1.0, 1.0, 1.0)
            self.update()

        else:
            super().initializeGL()

    @property
    def events(self):
        return self._events

    def eventFilter(self, obj, event):
        """
        Event handling (catching) function (overloaded)
        :param obj: DrawerGL
        :param event: QEvent object with event info
        :type event: QtCore.QEvent
        :return:
        """

        for wid_event in self._events:
            wid_event.assertEvent(event, self)

        self.update()
        return super().eventFilter(obj, event)


DrawWidgetTypes = loadUiType(r'.\Source\Extensions\ChemPackSource\ui\Drawer_model_ui.ui')


class DrawWidget(DrawWidgetTypes[1], DrawWidgetTypes[0]):

    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle('Subgraph search')
        self.parent = parent
        self.openGl_drawer = DrawerGL(parent=self.frame_7)
        self.verticalLayout_6.addWidget(self.openGl_drawer)

        self.organizer = None

        self.contacts_model = TableModel(parent=self,
                                         label=lambda x: self.openGl_drawer.drawing.labelContact(x).__getitem__(0),
                                         data=lambda x: self.openGl_drawer.drawing.labelContact(x).__getitem__(1),
                                         setter=lambda ind, dist: self.openGl_drawer.drawing.setDataContact(ind, dist))

        self.angle_model = TableModel(parent=self,
                                      label=lambda x: self.openGl_drawer.drawing.labelAngle(x).__getitem__(0),
                                      data=lambda x: self.openGl_drawer.drawing.labelAngle(x).__getitem__(1),
                                      setter=lambda ind, angle: self.openGl_drawer.drawing.setDataAngle(ind, angle))

        self.avgDiff_model = TableModel(parent=self,
                                         label=lambda x: self.openGl_drawer.drawing.labelAvgDiff(x).__getitem__(0),
                                         data=lambda x: self.openGl_drawer.drawing.labelAvgDiff(x).__getitem__(1),
                                         setter=lambda ind, dist: self.openGl_drawer.drawing.setDataAvgDiff(ind, dist))

        self.maxMeanPlaneDiff_model = TableModel(parent=self,
                                        label=lambda x: self.openGl_drawer.drawing.labelMaxMeanPlaneDiff(x).__getitem__(0),
                                        data=lambda x: self.openGl_drawer.drawing.labelMaxMeanPlaneDiff(x).__getitem__(1),
                                        setter=lambda ind, dist: self.openGl_drawer.drawing.setDataMaxMeanPlaneDiff(ind, dist))

        self.contacts_view = TableView(remove=lambda x: self.openGl_drawer.drawing.removeContact(x), parent=self, model=self.contacts_model)

        self.angle_view = TableView(remove=lambda x: self.openGl_drawer.drawing.removeAngle(x), parent=self, model=self.angle_model)

        self.avgDiff_view = TableView(remove=lambda x: self.openGl_drawer.drawing.removeAvgDiff(x), parent=self, model=self.avgDiff_model)

        self.maxMeanPlaneDiff_view = TableView(remove=lambda x: self.openGl_drawer.drawing.removeMaxMeanPlaneDiff(x), parent=self, model=self.maxMeanPlaneDiff_model)

        self.openGl_drawer.newContact.connect(lambda: self.contacts_model.insertRows(self.contacts_model.rowCount(), 1))
        self.openGl_drawer.newAngle.connect(lambda: self.angle_model.insertRows(self.angle_model.rowCount(), 1))
        self.openGl_drawer.newAvgDiff.connect(lambda: self.avgDiff_model.insertRows(self.avgDiff_model.rowCount(), 1))
        self.openGl_drawer.newMaxMeanPlaneDiff.connect(lambda: self.maxMeanPlaneDiff_model.insertRows(self.maxMeanPlaneDiff_model.rowCount(), 1))

        self.verticalLayout_8.addWidget(self.contacts_view)
        self.verticalLayout_9.addWidget(self.angle_view)
        self.verticalLayout_9.addWidget(self.angle_view)
        self.verticalLayout_11.addWidget(self.avgDiff_view)
        self.verticalLayout_10.addWidget(self.maxMeanPlaneDiff_view)

        HighLight().attach(self.openGl_drawer)
        drag_event = Drag(exclusive=True, exclusive_group=1)
        drag_event.attach(self.openGl_drawer)
        self.draw_event = Draw(exclusive=True, exclusive_group=1)
        clear_event = Clear(exclusive=True, exclusive_group=1)
        contact_event = Contact(exclusive=True, exclusive_group=1)
        angle_event = Angle(exclusive=True, exclusive_group=1)
        avgDiff_event = AvgDiff(exclusive=True, exclusive_group=1)
        maxMeanPlaneDiff_event = MaxMeanPlaneDiff(exclusive=True, exclusive_group=1)

        self.pushButton_7.pressed.connect(self.exportTable)
        self.pushButton_10.pressed.connect(lambda: drag_event.attach(self.openGl_drawer))
        self.pushButton.pressed.connect(lambda: clear_event.attach(self.openGl_drawer))
        self.pushButton_9.pressed.connect(lambda: contact_event.attach(self.openGl_drawer))
        self.pushButton_8.pressed.connect(lambda: angle_event.attach(self.openGl_drawer))
        self.pushButton_12.pressed.connect(lambda: avgDiff_event.attach(self.openGl_drawer))
        self.pushButton_13.pressed.connect(lambda: maxMeanPlaneDiff_event.attach(self.openGl_drawer))
        self.pushButton_6.pressed.connect(lambda: self.findSub(self.openGl_drawer.drawing))
        self.dbSearchButton.pressed.connect(lambda: self.dbSearch(self.openGl_drawer.drawing))

        self.pushButton_5.pressed.connect(lambda: self.changeDraw('C'))
        self.pushButton_4.pressed.connect(lambda: self.changeDraw('H'))
        self.pushButton_3.pressed.connect(lambda: self.changeDraw('N'))
        self.pushButton_2.pressed.connect(lambda: self.changeDraw('O'))
        self.pushButton_11.pressed.connect(self.setCustomAtomType)
        #self.pushButton_9.pressed.connect(lambda: self.checked_func.change_attr(func=1))
        #self.pushButton_8.pressed.connect(lambda: self.checked_func.change_attr(func=2))
        #self.dbSearchButton.pressed.connect(lambda: self.db_search(self.openGl_drawer.drawing.struct))

        #self.pushButton_6.pressed.connect(lambda: parent.substr_search(self.openGl_drawer.drawing.struct,
        #                                                               self.openGl_drawer.drawing.conditions))
        #self.pushButton_7.pressed.connect(self.openGl_drawer.save_template)
        self.gridLayout.addWidget(self.pushButton_6, 0, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_10, 1, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton, 2, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_11, 6, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_9, 7, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_8, 8, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_12, 9, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_13, 10, 0, 1, 2)

        #self.tableWidget_2.itemChanged.connect(lambda x: self.asd(x, tab=1))
        #self.tableWidget_3.itemChanged.connect(lambda x: self.asd(x, tab=2))

    def changeDraw(self, atom_type):
        self.draw_event.setType(atom_type)
        self.draw_event.attach(self.openGl_drawer)

    def setCustomAtomType(self):
        dialog = QtWidgets.QDialog(parent=self)
        dialog.setWindowTitle('Custom atom type')
        lineEd = QtWidgets.QLineEdit(parent=dialog)
        dialog.setLayout(QtWidgets.QVBoxLayout())
        dialog.layout().addWidget(lineEd)
        lineEd.returnPressed.connect(dialog.accept)
        func = lambda: self.getCustomAtomType(lineEd)
        dialog.accepted.connect(func)
        dialog.show()

    def getCustomAtomType(self, lineEd: QtWidgets.QLineEdit):
        text = lineEd.text()
        if text == '':
            self.changeDraw('H')
        types = text.split(',')
        self.changeDraw(types)

    def findSub(self, draw: Drawing):
        from ..ChemPack import MOLECULE_SYSTEMS, TREE_MODEL
        from .contacts import Pack, Node, findSubGraph, Condition
        from .MoleculeClass import Atom
        mol_sys = None

        for mol_l in MOLECULE_SYSTEMS:
            if mol_l.pick == 1.0 and mol_l.isValid():
                mol_sys = MOLECULE_SYSTEMS[mol_l]
                break

        if mol_sys is not None:

            sub_packs = draw.pack.splitPack()

            points = draw.points.children[0].children.copy()
            atoms = {}
            for point in points:
                atoms[point] = Atom(coord=np.array([0,0,0]), atom_type=point.atom_type, name=point.label)
            conditions = []
            a = [conditions.__iadd__(draw.conditions_d[x]) for x in draw.conditions_d]

            atoms = {}
            bonds = {}
            for atom in mol_sys.children[0].children:
                atoms[atom] = atom
                bonds[atom] = []
                for bond in atom.bonds():
                    atom2 = bond.parents()[bond.parents().index(atom)-1]
                    bonds[atom].append(atom2)
                    a2_bonds = bonds.get(atom2, None)
                    if a2_bonds is None:
                        bonds[atom2] = [atom]
                    else:
                        if atom not in a2_bonds:
                            a2_bonds.append(atom)

            points = mol_sys.children[0].children.copy()
            pack = Pack()
            mem = {}

            def rec(point, pack, mem, bonds):
                if point not in mem:
                    points.remove(point)
                    node = Node(pack, str_atom=atoms[point])
                    mem[point] = node
                else:
                    return mem[point]
                for point_2 in bonds[point]:
                    node_2 = rec(point_2, pack, mem, bonds)
                    node.addConnect(node_2)
                return node

            while points:
                rec(points[0], pack, mem, bonds)
            a = []
            for sub_pack in sub_packs:
                a.append(findSubGraph(pack, sub_pack))
            result = [[x[2] for x in y if x[2]] for y in a]
            for i, res in enumerate(result):
                sub_pack = sub_packs[i]
                for rnode_l in res:
                    for rnode in rnode_l:
                        sub_pack.addSolution(rnode.pack)

            if conditions:
                ret = Condition.intersection(*conditions, eval=True)
                for cond in ret:
                    for sub_pack in cond:
                        for node in sub_pack.nodes:
                            for sol in cond[sub_pack]:
                                point = node.real_atoms[sol].point()
                                if point.pick is None:
                                    point.addProperty('pick', 1.0)
                                    continue
                                node.real_atoms[sol].point().pick = 1.0
            else:
                for sp in a:
                    for sol in sp:
                        if sol[2]:
                            for node in sol[2][0].pack.nodes:
                                point = node.assigned_node.struct_atom.point()
                                if point.pick is None:
                                    point.addProperty('pick', 1.0)
                                    continue
                                point.pick = 1.0

            def visCond(condition_type, conditions, result, root):
                if TREE_MODEL is None:
                    return
                if condition_type == 'contacts':
                    contacts_list = point_class.PointsList(parent=root, name='Contacts', color=np.array([0, 0, 0, 1], dtype=np.float32))
                    cond_l = point_class.PointsList(parent=contacts_list, color=contacts_list, name='Pairs')
                    labels_l = point_class.PointsList(parent=contacts_list, color=contacts_list, name='Labels')
                    for i in range(len(conditions)):
                        cond = conditions[i]
                        res = result[i]
                        length = len(res[list(res.keys())[0]])
                        ind1 = cond.packs_list.index(cond.nodes[0].pack)
                        ind2 = cond.packs_list.index(cond.nodes[1].pack)
                        for sol in range(length):
                            sol1 = res[cond.packs_list[ind1]][sol]
                            sol2 = res[cond.packs_list[ind2]][sol]
                            value = round(cond.value_func(cond.nodes[0].real_atoms[sol1], cond.nodes[1].real_atoms[sol2]), 2)
                            point1 = cond.nodes[0].real_atoms[sol1].point()
                            point2 = cond.nodes[1].real_atoms[sol2].point()

                            cl = point_class.PointsList(parent=cond_l, color=cond_l,
                                                        name=f'{point1.name}_{point2.name}')
                            lp = point_class.Point(parent=labels_l, coord=(point1.coord+point2.coord)/2, color=contacts_list, name=cl, label=str(value))
                            pc1 = point_class.Point(parent=cl, coord=point1, color=cl, name=point1, label=point1)
                            pc2 = point_class.Point(parent=cl, coord=point2, color=cl, name=point2, label=point2)
                    mol_ind = TREE_MODEL.index(-1, 0, by_point=root)
                    TREE_MODEL.insertRow(TREE_MODEL.rowCount(parent=mol_ind), parent=mol_ind)
                    cont_ind = TREE_MODEL.index(TREE_MODEL.rowCount(parent=mol_ind)-1, 0, parent=mol_ind)
                    ind = TREE_MODEL.index(0, 0, parent=cont_ind)
                    TREE_MODEL.attachObserver(ind, 'Dashed line')
                if condition_type == 'angle':
                    angles_cond_l = point_class.PointsList(parent=root, name='Angles', color=np.array([0, 0, 0, 1], dtype=np.float32))
                    angles_l = point_class.PointsList(parent=angles_cond_l, color=np.array([0, 1, 0, 0.5], dtype=np.float32), name='Angles')
                    planes_l = point_class.PointsList(parent=angles_cond_l, color=np.array([1, 0, 0, 0.5], dtype=np.float32), name='Planes')
                    labels_l = point_class.PointsList(parent=angles_cond_l, color=np.array([0, 0, 0, 1], dtype=np.float32), name='Labels')
                    for i in range(len(conditions)):
                        cond = conditions[i]
                        res = result[i]
                        length = len(res[list(res.keys())[0]])
                        if cond.split == -1:
                            for sol in range(length):
                                a1, a2, a3 = [cond.nodes[x].real_atoms[res[cond.packs_list[cond.packs_list.index(cond.nodes[x].pack)]][sol]] for x in range(3)]
                                p1, p2, p3 = a1.point(), a2.point(), a3.point()
                                vec1, vec2 = p1.coord - p2.coord, p3.coord - p2.coord
                                vec1, vec2 = 1.5 * vec1 / np.linalg.norm(vec1), 1.25 * vec2 / np.linalg.norm(vec2)
                                coord1, coord2 = p2.coord + vec1, p2.coord + vec2
                                plane = point_class.PointsList(parent=angles_l, name=f'{p1.name}_{p2.name}_{p3.name}', color=angles_l)
                                pp1 = point_class.Point(parent=plane, color=plane, coord=p2.coord+vec1)
                                pp2 = point_class.Point(parent=plane, color=plane, coord=p2)
                                pp3 = point_class.Point(parent=plane, color=plane, coord=p2.coord+vec2)
                                value = round(cond.value_func(a1,a2,a3), 1)
                                lp = point_class.Point(parent=labels_l, color=labels_l, name=f'{p1.name}_{p2.name}_{p3.name}', coord=(pp1.coord+pp3.coord)/2, label=str(value))
                        else:
                            from .contacts import meanPlaneEq
                            for sol in range(length):
                                atoms = [cond.nodes[x].real_atoms[
                                                  res[cond.packs_list[cond.packs_list.index(cond.nodes[x].pack)]][sol]]
                                              for x in range(len(cond.nodes))]
                                mean_plane1 = meanPlaneEq(*[x.point() for x in atoms[:cond.split]])
                                mean_plane2 = meanPlaneEq(*[x.point() for x in atoms[cond.split:]])

                                max_coord = np.array([max([y.point().coord[x] for y in atoms]) for x in range(3)], dtype=np.float32)
                                min_coord = np.array([min([y.point().coord[x] for y in atoms]) for x in range(3)], dtype=np.float32)
                                for mean_plane, split in zip((mean_plane1, mean_plane2), ((0, cond.split), (cond.split, len(atoms)))):
                                    a, b, c, d = mean_plane
                                    norm = np.array([a,b,c], dtype=np.float32)
                                    coord_max = max_coord - norm * ((a * max_coord[0] + b * max_coord[1] + c * max_coord[2] + d) / np.linalg.norm(norm)**2)
                                    coord_min = min_coord - norm * ((a * min_coord[0] + b * min_coord[1] + c * min_coord[2] + d) / np.linalg.norm(norm)**2)
                                    vec1 = coord_max - coord_min
                                    vec2 = np.cross(norm, vec1)
                                    c1 = coord_max + np.linalg.norm(vec1)*(vec2/np.linalg.norm(vec2))
                                    c2 = coord_max - np.linalg.norm(vec1)*(vec2/np.linalg.norm(vec2))
                                    c3 = coord_min - np.linalg.norm(vec1)*(vec2/np.linalg.norm(vec2))
                                    c4 = coord_min + np.linalg.norm(vec1)*(vec2/np.linalg.norm(vec2))

                                    p1_l = point_class.PointsList(parent=planes_l, color=planes_l, name='_'.join([x.point().name for x in atoms[split[0]:split[1]]]))
                                    for c in (c1,c2,c3,c1,c4,c3):
                                        point = point_class.Point(parent=p1_l, coord=c, color=planes_l)

                    mol_ind = TREE_MODEL.index(-1, 0, by_point=root)
                    TREE_MODEL.insertRow(TREE_MODEL.rowCount(parent=mol_ind), parent=mol_ind)
                    cont_ind = TREE_MODEL.index(TREE_MODEL.rowCount(parent=mol_ind) - 1, 0, parent=mol_ind)
                    ind = TREE_MODEL.index(0, 0, parent=cont_ind)
                    TREE_MODEL.attachObserver(ind, 'Plane')
                    ind = TREE_MODEL.index(1, 0, parent=cont_ind)
                    TREE_MODEL.attachObserver(ind, 'Plane')
                if condition_type == 'avgDiff':
                    avgDiff_list = point_class.PointsList(parent=root, name='avgDiff',
                                                          color=np.array([0, 0, 0, 1], dtype=np.float32))
                    cent_list = point_class.PointsList(parent=avgDiff_list, name='Cent', color=np.array([1, 0, 1, 1], dtype=np.float32), rad=0.15)
                    cond_l = point_class.PointsList(parent=avgDiff_list, color=avgDiff_list, name='Pairs')
                    labels_l = point_class.PointsList(parent=avgDiff_list, color=avgDiff_list, name='Labels')
                    for i in range(len(conditions)):
                        cond = conditions[i]
                        res = result[i]
                        length = len(res[list(res.keys())[0]])
                        for sol in range(length):
                            atoms = [cond.nodes[x].real_atoms[res[cond.packs_list[cond.packs_list.index(cond.nodes[x].pack)]][sol]] for x in range(len(cond.nodes))]
                            points = [x.point() for x in atoms]
                            value = round(cond.value_func(*atoms), 2)
                            coord1 = np.sum([x.coord for x in atoms[:cond.split]], axis=0)/cond.split
                            cent_point1 = point_class.Point(parent=cent_list, name='_'.join([str(x.name) for x in points]), color=cent_list, rad=cent_list, coord=coord1)
                            coord2 = np.sum([x.coord for x in atoms[cond.split:]], axis=0)/(len(atoms)-cond.split)
                            cent_point2 = point_class.Point(parent=cent_list, name='_'.join([str(x.name) for x in points]), color=cent_list, rad=cent_list, coord=coord2)

                            cl = point_class.PointsList(parent=cond_l, color=cond_l,
                                                        name=f'{cent_point1.name}_{cent_point2.name}')
                            lp = point_class.Point(parent=labels_l, coord=(cent_point1.coord + cent_point2.coord) / 2,
                                                   color=avgDiff_list, name=cl, label=str(value))
                            pc1 = point_class.Point(parent=cl, coord=cent_point1, color=cl, name=cent_point1, label=cent_point1)
                            pc2 = point_class.Point(parent=cl, coord=cent_point2, color=cl, name=cent_point2, label=cent_point2)
                    mol_ind = TREE_MODEL.index(-1, 0, by_point=root)
                    TREE_MODEL.insertRow(TREE_MODEL.rowCount(parent=mol_ind), parent=mol_ind)
                    cont_ind = TREE_MODEL.index(TREE_MODEL.rowCount(parent=mol_ind) - 1, 0, parent=mol_ind)
                    ind = TREE_MODEL.index(0, 0, parent=cont_ind)
                    TREE_MODEL.attachObserver(ind, 'Sphere')
                    ind = TREE_MODEL.index(1, 0, parent=cont_ind)
                    TREE_MODEL.attachObserver(ind, 'Dashed line')
                if condition_type == 'maxMeanPlaneDiff':
                    pass

            if TREE_MODEL is not None:
                molecule = mol_sys.children[0]
                points_list = molecule.point()
                i = 0
                for cond_type in draw.conditions_d:
                    if draw.conditions_d[cond_type]:
                        visCond(cond_type, draw.conditions_d[cond_type], ret[i:i+len(draw.conditions_d[cond_type])], points_list)
                    i += len(draw.conditions_d[cond_type])
                self.organizer = SolutionOrganizer(*conditions)
                self.organizer.fillTable()
            pack = Pack.concatenatePacks(*sub_packs)
            draw.pack = pack

    def dbSearch(self, draw):
        from . import Db_viewer
        Db_viewer.show()
        Db_viewer.DB_VIEWER.list_model.populate((draw, 'substructure'))

    def exportTable(self):
        if self.organizer is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(parent=self, filter='*.csv')[0]
            try:
                out = open(filename, 'w')
                table = self.organizer.getTable()
                out.write(table)
                out.close()
            except PermissionError:
                pass


def show():
    global DRAW_WIDGET
    if DRAW_WIDGET is None:
        DRAW_WIDGET = DrawWidget()
    DRAW_WIDGET.show()


if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)
    window = DrawWidget()
    window.show()
    app.exec_()
