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
from PySide6.QtOpenGLWidgets import QOpenGLWidget

from .Db_viewer import DbWindow
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
import logging

import debug

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
        return self.conditions_d['contacts'][-1]

    def addAngle(self, *points, split=-1):
        self.conditions_value_d['angle'].append(0.0)
        self.conditions_d['angle'].append(Condition(lambda x: x >= 0.0 or x <= 0.0, lambda *x: angle(x, rad=False, split=split), [self.point_node[x] for x in points]))
        self.conditions_d['angle'][-1].split = split
        return self.conditions_d['angle'][-1]

    def addAvgDiff(self, *points, split=1):
        self.conditions_value_d['avgDiff'].append(0.0)
        self.conditions_d['avgDiff'].append(
            Condition(lambda x: x <= 0.0, lambda *x: avgDiff(x, split=split), [self.point_node[x] for x in points]))
        self.conditions_d['avgDiff'][-1].split = split
        return self.conditions_d['avgDiff'][-1]

    def addMaxMeanPlaneDiff(self, *points):
        self.conditions_value_d['maxMeanPlaneDiff'].append(0.0)
        self.conditions_d['maxMeanPlaneDiff'].append(
            Condition(lambda x: x <= 0.0, lambda *x: maxMeanPlaneDiff(*x), [self.point_node[x] for x in points]))
        return self.conditions_d['maxMeanPlaneDiff'][-1]

    def removeAngle(self, ind):
        if ind == -1:
            return
        if type(ind) is Condition:
            ind = self.conditions_d['angle'].index(ind)
        cond = self.conditions_d['angle'].pop(ind)
        self.conditions_value_d['angle'].pop(ind)
        return cond

    def removeContact(self, ind):
        if ind == -1:
            return
        if type(ind) is Condition:
            ind = self.conditions_d['contacts'].index(ind)
        cond = self.conditions_d['contacts'].pop(ind)
        node1, node2 = cond.nodes
        point1, point2 = self.node_point[node1], self.node_point[node2]
        self.conditions_value_d['contacts'].pop(ind)
        self.contacts[point1][point2].parent.destroy()
        self.contacts[point1].pop(point2)
        self.contacts[point2].pop(point1)
        return cond

    def removeAvgDiff(self, ind):
        if ind == -1:
            return
        if type(ind) is Condition:
            ind = self.conditions_d['avgDiff'].index(ind)
        cond = self.conditions_d['avgDiff'].pop(ind)
        self.conditions_value_d['avgDiff'].pop(ind)
        return cond

    def removeMaxMeanPlaneDiff(self, ind):
        if ind == -1:
            return
        if type(ind) is Condition:
            ind = self.conditions_d['maxMeanPlaneDiff'].index(ind)
        cond = self.conditions_d['maxMeanPlaneDiff'].pop(ind)
        self.conditions_value_d['maxMeanPlaneDiff'].pop(ind)
        return cond

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
        for cond_key in self.conditions_d:
            conds = self.conditions_d[cond_key]
            for cond in conds:
                if node in cond.nodes:
                    i = self.conditions_d[cond_key].index(cond)
                    self.conditions_d[cond_key].pop(i)
                    self.conditions_value_d[cond_key].pop(i)
        self.node_point.pop(node)
        node.delete()
        point.destroy()

    def removeBond(self, point1, point2):
        if point1 in self.connections:
            if point2 in self.connections[point1]:
                self.connections[point1].pop(point2)
        if point2 in self.connections:
            if point1 in self.connections[point2]:
                self.connections[point2].pop(point1)
        self.point_node[point1].removeConnect(self.point_node[point2])

    def show_cn(self):
        for point in self.point_node.keys():
            point.label = point.cn_label

    def show_at(self):
        for point in self.point_node.keys():
            point.label = point.at_label


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

    def __init__(self, remove, parent=None, model=None, widget=None):

        super().__init__(parent)
        self.__widget = widget
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
        com = self._remove(self.__widget, self.model())
        com.apply(self.currentIndex().row())
        com.appendStack(com)
        #self._remove(self.currentIndex().row())
        #self.model().removeRows(self.currentIndex().row(), 1)

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


class Command(ABC):

    stack = []
    reverse_stack = []

    def __init__(self, *args, **kwargs):
        self.applied = False
        return

    @abstractmethod
    def apply(self, *args, **kwargs):
        self.applied = True
        pass

    @abstractmethod
    def undo(self):
        pass

    @staticmethod
    def popStack():
        try:
            command = Command.stack.pop()
        except IndexError:
            return
        command.undo()
        Command.reverse_stack.append(command)

    @staticmethod
    def popReverseStack():
        try:
            command = Command.reverse_stack.pop()
        except IndexError:
            return
        command.apply()
        Command.stack.append(command)

    @staticmethod
    def appendStack(command):
        Command.stack.append(command)
        Command.reverse_stack = []


class DragCommand(Command):

    def __init__(self):
        Command.__init__(self)
        self.pos = None
        self.old_pos = None
        self.cc_point = None

    def payload(self, cc_point, pos, *args, **kwargs):
        cc_point.point.coord = pos

    def apply(self, cc_point, old_pos, pos, *args, **kwargs):
        Command.apply(self)
        logging.debug(f'Applied {self}, atom: {cc_point.point.label}')
        self.apply = lambda *args, **kwargs: self.payload(cc_point, pos, *args, **kwargs)
        self.apply()
        self.cc_point = cc_point
        self.pos = pos
        self.old_pos = old_pos

    def undo(self):
        if self.applied:
            self.payload(self.cc_point, self.old_pos)


class CreateAtomCommand(Command):

    def __init__(self, drawing, points_list, ids):
        Command.__init__(self)
        self.drawing = drawing
        self.point_list = points_list
        self.ids = ids
        self.kwargs = None
        self.point = None

    def payload(self, *args, **kwargs):
        id = self.ids.get()
        self.point = point_class.Point(parent=self.point_list,
                                       coord=kwargs['coord'],
                                       rad=self.point_list,
                                       label=kwargs['label'] + str(id),
                                       at_label=kwargs['label'] + str(id),
                                       cn_label='0-14',
                                       label_size=6,
                                       id=id,
                                       atom_type=kwargs['atom_type'],
                                       color=kwargs['color'],
                                       create_command=kwargs['create_command'],
                                       cn=[0, 14])
        self.drawing.add_point(self.point)
        kwargs['create_command'].point = self.point
        return self.point

    def apply(self, *args, **kwargs):
        Command.apply(self)
        self.kwargs = kwargs
        self.apply = lambda *args, **kwargs: self.payload(*args,  create_command=self, **self.kwargs)
        self.apply()
        logging.debug(f'Applied {self}, atom: {self.point.label}')
        return self.point

    def undo(self):
        if self.applied:
            DeleteAtomCommand(self.drawing, self.point_list, self.ids).apply(self)


class DeleteAtomCommand(Command):

    def __init__(self, drawing, points_list, ids):
        Command.__init__(self)
        self.drawing = drawing
        self.point_list = points_list
        self.ids = ids
        self.kwargs = None
        self.cc_point = None
        self.bond_commands = []

    def payload(self, cc_point, *args, **kwargs):
        for point2 in self.drawing.connections[cc_point.point].copy():
            self.bond_commands.append(DeleteBondCommand(self.drawing, self.drawing.connections[cc_point.point][point2].parent))
            self.bond_commands[-1].apply(cc_point, point2.create_command)
        self.ids.remove(cc_point.point.id)
        self.drawing.remove_point(cc_point.point)
        cc_point.point = None

    def apply(self, cc_point, *args, **kwargs):
        Command.apply(self)
        point = cc_point.point
        logging.debug(f'Applied {self}, atom: {point.label}')
        self.kwargs = {'coord': point.coord.copy(),
                       'label': ''.join([x for x in point.label if x.isalpha()]),
                       'id': point.id,
                       'atom_type': point.atom_type,
                       'color': point.color.copy(),
                       'create_command': point.create_command}
        self.cc_point = cc_point
        self.apply = lambda *args, **kwargs: self.payload(self.cc_point)
        self.apply()

    def undo(self):
        if self.applied:
            self.point = CreateAtomCommand(self.drawing, self.point_list, self.ids).payload(**self.kwargs)
            for bond_command in self.bond_commands:
                bond_command.undo()
            self.bond_commands = []


class CreateBondCommand(Command):

    def __init__(self, drawing, line_list):
        Command.__init__(self)
        self.drawing = drawing
        self.line_list = line_list
        self.cc_point1 = None
        self.cc_point2 = None

    def payload(self, cc_point1, cc_point2, *args, **kwargs):
        point1 = cc_point1.point
        point2 = cc_point2.point
        line1 = point_class.Point(parent=self.line_list, coord=point1, color=point1, rad=point1)
        line2 = point_class.Point(parent=self.line_list, coord=point2, color=point2, rad=point2)
        self.drawing.add_connection((point1, line1), (point2, line2))

    def apply(self, cc_point1, cc_point2, *args, **kwargs):
        Command.apply(self)
        self.cc_point1 = cc_point1
        self.cc_point2 = cc_point2
        logging.debug(f'Applied {self}, atoms: {self.cc_point1.point.label}, {self.cc_point2.point.label}')
        self.apply = lambda *args, **kwargs: self.payload(cc_point1, cc_point2)
        self.apply()

    def undo(self):
        if self.applied:
            DeleteBondCommand(self.drawing, self.line_list).apply(self.cc_point1, self.cc_point2)


class DeleteBondCommand(Command):

    def __init__(self, drawing, line_list):
        Command.__init__(self)
        self.drawing = drawing
        self.line_list = line_list
        self.cc_point1 = None
        self.cc_point2 = None

    def payload(self, cc_point1, cc_point2, *args, **kwargs):
        point1 = cc_point1.point
        point2 = cc_point2.point
        self.drawing.removeBond(point1, point2)
        self.line_list.children[0].destroy()
        self.line_list.children[0].destroy()

    def apply(self, cc_point1,  cc_point2, *args, **kwargs):
        Command.apply(self)
        self.cc_point1 = cc_point1
        self.cc_point2 = cc_point2
        logging.debug(f'Applied {self}, atoms: {self.cc_point1.point.label}, {self.cc_point2.point.label}')
        self.apply = lambda *args, **kwargs: self.payload(cc_point1, cc_point2)
        self.apply()

    def undo(self):
        if self.applied:
            CreateBondCommand(self.drawing, self.line_list).apply(self.cc_point1, self.cc_point2)


class ChangeTypeCommand(Command):

    def __init__(self, drawing):
        Command.__init__(self)
        self.old_type = None
        self.old_label = None
        self.old_color = None
        self.cc_point = None
        self.drawing = drawing

    def payload(self, cc_point, atom_type, label, color, *args, **kwargs):
        cc_point.point.label = label
        cc_point.point.color = color
        cc_point.point.atom_type = atom_type
        self.drawing.changeAtomType(cc_point.point)

    def apply(self, cc_point, atom_type, label, color, *args, **kwargs):
        Command.apply(self)
        self.cc_point = cc_point
        self.old_type = cc_point.point.atom_type
        self.old_label = cc_point.point.label
        self.old_color = cc_point.point.color
        logging.debug(f'Applied {self}, atoms: {self.cc_point.point.label}')
        self.apply = lambda *args, **kwargs: self.payload(self.cc_point, atom_type, label, color)
        self.apply()

    def undo(self):
        if self.applied:
            self.payload(self.cc_point, self.old_type, self.old_label, self.old_color)


class CompositeCommand(Command):

    def __init__(self, *commands):
        Command.__init__(self)
        self.commands = commands
        self.applied = True

    def apply(self, *args, **kwargs):
        Command.apply(self)
        logging.debug(f'Applied {self}')
        for command in self.commands:
            if command.applied:
                command.apply()

    def undo(self):
        if self.applied:
            for command in self.commands[::-1]:
                command.undo()


class ReverseCompositeCommand(Command):

    def __init__(self, *commands):
        Command.__init__(self)
        self.commands = commands

    def apply(self, *args, **kwargs):
        Command.apply(self)
        logging.debug(f'Applied {self}')
        for command in self.commands[::-1]:
            if command.applied:
                command.undo()

    def undo(self):
        if self.applied:
            for command in self.commands:
                command.apply()


class ConditionCommand(Command):

    def __init__(self, widget):
        Command.__init__(self)
        self.cc_points = []
        self.widget = widget
        self.condition = None

    @abstractmethod
    def payload(self, *args, **kwargs):
        pass

    def removePayload(self, *args, **kwargs):
        pass

    def apply(self, *cc_points, **kwargs):
        Command.apply(self)
        self.cc_points = cc_points
        self.apply = lambda *args, **kwargs: self.payload(*self.cc_points)
        self.apply()

    def undo(self):
        if self.applied:
            self.removePayload(self.cc_points)


class ContactsCommand(ConditionCommand):

    def __init__(self, widget, line_list):
        ConditionCommand.__init__(self, widget)
        self.line_list = line_list
        self.__model = DRAW_WIDGET.contacts_model
        self.ind = None

    def payload(self, *cc_points, **kwargs):
        pc1, pc2 = cc_points[0].point, cc_points[1].point
        dlp1 = point_class.Point(parent=self.line_list, color=pc1, coord=pc1, rad=pc1)
        dlp2 = point_class.Point(parent=self.line_list, color=pc2, coord=pc2, rad=pc2)
        self.condition = self.widget.drawing.add_contact((pc1, dlp1), (pc2, dlp2))
        self.condition.create_command = self
        self.ind = self.widget.drawing.conditions_d['contacts'].index(self.condition)
        self.widget.newContact.emit()

    def removePayload(self, cc_points, *args, **kwargs):
        if self.applied:
            self.__model.removeRows(self.ind, 1)
            self.widget.drawing.removeContact(self.ind)


class DeleteContactsCommand(ConditionCommand):

    def __init__(self, widget, model):
        ConditionCommand.__init__(self, widget)
        self.__model = model
        self.__cc_points = []
        self.__line_list = None

    def payload(self, ind, *args, **kwargs):
        cond = self.widget.drawing.removeContact(ind)
        self.__cc_points = [self.widget.drawing.node_point[x].create_command for x in cond.nodes]
        self.__line_list = self.__cc_points[0].point.parent
        self.__model.removeRows(ind, 1)

    def removePayload(self, *args, **kwargs):
        if self.applied:
            ContactsCommand(self.widget, self.__line_list).apply(*self.__cc_points)


class AngleCommand(ConditionCommand):

    def __init__(self, widget):
        ConditionCommand.__init__(self, widget)
        self.split = -1
        self.__model = DRAW_WIDGET.angle_model
        self.ind = None

    def payload(self, *cc_points, **kwargs):
        points = [x.point for x in cc_points]
        self.condition = self.widget.drawing.addAngle(*points, split=self.split)
        self.condition.create_command = self
        self.ind = self.widget.drawing.conditions_d['angle'].index(self.condition)
        self.widget.newAngle.emit()

    def apply(self, *cc_points, **kwargs):
        self.split = kwargs.get('split', -1)
        ConditionCommand.apply(self, *cc_points, **kwargs)

    def removePayload(self, cc_points, *args, **kwargs):
        if self.applied:
            self.__model.removeRows(self.ind, 1)
            self.widget.drawing.removeAngle(self.ind)


class DeleteAngleCommand(ConditionCommand):

    def __init__(self, widget, model):
        ConditionCommand.__init__(self, widget)
        self.split = -1
        self.__model = model
        self.__cc_points = []

    def payload(self, ind, *args, **kwargs):
        cond = self.widget.drawing.removeAngle(ind)
        self.split = cond.split
        self.__cc_points = [self.widget.drawing.node_point[x].create_command for x in cond.nodes]
        self.__model.removeRows(ind, 1)

    def removePayload(self, *args, **kwargs):
        if self.applied:
            AngleCommand(self.widget).apply(*self.__cc_points, split=self.split)


class AvgDiffCommand(ConditionCommand):

    def __init__(self, widget):
        ConditionCommand.__init__(self, widget)
        self.split = -1
        self.__model = DRAW_WIDGET.avgDiff_model
        self.ind = None

    def payload(self, *cc_points, **kwargs):
        points = [x.point for x in cc_points]
        self.condition = self.widget.drawing.addAvgDiff(*points, split=self.split)
        self.condition.create_command = self
        self.ind = self.widget.drawing.conditions_d['avgDiff'].index(self.condition)
        self.widget.newAvgDiff.emit()

    def apply(self, *cc_points, **kwargs):
        self.split = kwargs.get('split', -1)
        ConditionCommand.apply(self, *cc_points, **kwargs)

    def removePayload(self, cc_points, *args, **kwargs):
        if self.applied:
            self.__model.removeRows(self.ind, 1)
            self.widget.drawing.removeAvgDiff(self.ind)


class DeleteAvgDiffCommand(ConditionCommand):

    def __init__(self, widget, model):
        ConditionCommand.__init__(self, widget)
        self.split = -1
        self.__model = model
        self.__cc_points = []

    def payload(self, ind, *args, **kwargs):
        cond = self.widget.drawing.removeAvgDiff(ind)
        self.split = cond.split
        self.__cc_points = [self.widget.drawing.node_point[x].create_command for x in cond.nodes]
        self.__model.removeRows(ind, 1)

    def removePayload(self, *args, **kwargs):
        if self.applied:
            AvgDiffCommand(self.widget).apply(*self.__cc_points, split=self.split)


class MaxMeanPlaneDiffCommand(ConditionCommand):

    def __init__(self, widget):
        ConditionCommand.__init__(self, widget)
        self.__model = DRAW_WIDGET.maxMeanPlaneDiff_model
        self.ind = None

    def payload(self, *cc_points, **kwargs):
        points = [x.point for x in cc_points]
        self.condition = self.widget.drawing.addMaxMeanPlaneDiff(*points)
        self.condition.create_command = self
        self.ind = self.widget.drawing.conditions_d['maxMeanPlaneDiff'].index(self.condition)
        self.widget.newMaxMeanPlaneDiff.emit()

    def removePayload(self, cc_points, *args, **kwargs):
        if self.applied:
            self.__model.removeRows(self.ind, 1)
            self.widget.drawing.removeMaxMeanPlaneDiff(self.ind)


class DeleteMaxMeanPlaneCommand(ConditionCommand):

    def __init__(self, widget, model):
        ConditionCommand.__init__(self, widget)
        self.__model = model
        self.__cc_points = []

    def payload(self, ind, *args, **kwargs):
        cond = self.widget.drawing.removeMaxMeanPlaneDiff(ind)
        self.__cc_points = [self.widget.drawing.node_point[x].create_command for x in cond.nodes]
        self.__model.removeRows(ind, 1)

    def removePayload(self, *args, **kwargs):
        if self.applied:
            MaxMeanPlaneDiffCommand(self.widget).apply(*self.__cc_points)


class ChangeCNCommand(Command):
    def __init__(self, drawing):
        Command.__init__(self)
        self.old_cn = None
        self.old_cn_label = None
        self.cc_point = None
        self.drawing = drawing

    def payload(self, cc_point, cn, cn_label, *args, **kwargs):
        cc_point.point.cn_label = cn_label
        cc_point.point.label = cn_label
        cc_point.point.cn = cn
        self.drawing.changeAtomType(cc_point.point)

    def apply(self, cc_point, cn, cn_label, *args, **kwargs):
        Command.apply(self)
        self.cc_point = cc_point

        self.old_cn = cc_point.point.cn
        self.old_cn_label = cc_point.point.cn_label

        self.apply = lambda *args, **kwargs: self.payload(self.cc_point, cn, cn_label)
        self.apply()

    def undo(self):
        if self.applied:
            self.payload(self.cc_point, self.old_cn, self.old_cn_label)


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
        self.old_pos = None
        self.pos = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.drag_point = widget.select(event.localPos())
            logging.debug(f'{self} Mouse click, self.pos: {self.pos}')
            if self.drag_point is not None:
                self.old_pos = self.drag_point.coord
                self.pos = self.drag_point.coord
                self.timer.start()

        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.LeftButton:
            if self.timer.run():
                self.timer.stop()
            if event.modifiers() == QtCore.Qt.NoModifier and self.timer.time() > self.tol:
                logging.debug(f'{self} Mouse move, self.pos: {self.pos}')
                if self.drag_point is not None:
                    self.pos = np.array([((event.localPos().x() / self.widget.width()) * 2 - 1),
                                                      ((event.localPos().y() / self.widget.height()) * (-2) + 1) / (
                                                                  self.widget.width() / self.widget.height()),
                                                      0], dtype=np.float32)
                    self.drag_point.coord = self.pos

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            logging.debug(f'{self} Mouse release, self.pos: {self.pos}')
            if self.drag_point is not None:
                com = DragCommand()
                com.apply(self.drag_point.create_command, self.old_pos, self.pos)
                com.appendStack(com)
            self.drag_point = None
            self.pos = None
            self.old_pos = None
            self.timer.stop()

    def detach(self):
        aDrawWidgetEvent.detach(self)
        self.drag_point = None
        self.old_pos = None
        self.pos = None
        self.timer.stop()
        logging.debug(f'{self} Drag detach, self.pos: {self.pos}')


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


class UndoRedo(aDrawWidgetEvent):

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.KeyPress:
            if event.matches(QtGui.QKeySequence.Undo):
                Command.popStack()
            elif event.matches(QtGui.QKeySequence.Redo):
                Command.popReverseStack()


class Draw(aDrawWidgetEvent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.atom_type = 'C'
        self.atom_type_label = 'C'
        self.atom_color = PALETTE.getColor('C')

        self.create_atom1_command = None
        self.create_atom2_command = None
        self.create_bond_command = None
        self.change_type_command = None
        self.drag_command = None

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
                if self.create_atom1_command is None and self.point1 is None:
                    self.create_atom1_command = CreateAtomCommand(self.widget.drawing, self.widget.p_list, self.widget.ids)
                    self.point1 = self.create_atom1_command.apply(coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()),
                                                                    0], dtype=np.float32),
                                                                  rad=self.widget.p_list,
                                                                  label=self.atom_type_label,
                                                                  atom_type=self.atom_type,
                                                                  color=self.atom_color)
                if self.create_atom2_command is None:
                    self.create_atom2_command = CreateAtomCommand(self.widget.drawing, self.widget.p_list, self.widget.ids)
                    self.point2 = self.create_atom2_command.apply(coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()),
                                                                    0], dtype=np.float32),
                                                                  rad=self.widget.p_list,
                                                                  color=self.atom_color,
                                                                  atom_type=self.atom_type,
                                                                  label=self.atom_type_label)

                if self.create_bond_command is None:
                    self.line_list = point_class.PointsList(parent=self.widget.l_list)
                    self.create_bond_command = CreateBondCommand(self.widget.drawing, self.line_list)
                    self.create_bond_command.apply(self.point1.create_command, self.point2.create_command)
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
                    self.create_bond_command.undo()
                    self.create_atom2_command.undo()
                    self.create_atom2_command = None
                    self.point2 = point2
                    self.line_list = point_class.PointsList(parent=self.widget.l_list)
                    self.create_bond_command = CreateBondCommand(self.widget.drawing, self.line_list)
                    self.create_bond_command.apply(self.point1.create_command, self.point2.create_command)

                    coms = [x for x in [self.create_atom1_command, self.create_atom2_command, self.create_bond_command, self.drag_command] if x]
                    com = CompositeCommand(*coms)
                    com.appendStack(com)

                elif self.point2 is None:
                    pos = event.localPos()
                    point_a = self.widget.select(pos)
                    if point_a is None:
                        self.create_atom1_command = CreateAtomCommand(self.widget.drawing, self.widget.p_list, self.widget.ids)
                        self.create_atom1_command.apply(coord=np.array([((pos.x() / self.widget.width()) * 2 - 1),
                                                                    ((pos.y() / self.widget.height()) * (-2) + 1) / (
                                                                                self.widget.width() / self.widget.height()), 0], dtype=np.float32),
                                                        color=self.atom_color,
                                                        rad=self.widget.p_list,
                                                        atom_type=self.atom_type,
                                                        label=self.atom_type_label)
                        self.create_atom1_command.appendStack(self.create_atom1_command)
                    else:
                        self.change_type_command = ChangeTypeCommand(self.widget.drawing)
                        label = self.atom_type_label + str(point_a.id)
                        color = self.atom_color
                        atom_type = self.atom_type
                        self.change_type_command.apply(point_a.create_command, atom_type, label, color)
                        self.change_type_command.appendStack(self.change_type_command)
                else:
                    self.drag_command = DragCommand()
                    self.drag_command.apply(self.point2.create_command, self.point2.coord, self.point2.coord)
                    coms = [x for x in [self.create_atom1_command, self.create_atom2_command, self.create_bond_command, self.drag_command] if x]
                    com = CompositeCommand(*coms)
                    com.appendStack(com)

                self.point1 = None
                self.point2 = None
                self.line_list = None
                self.create_atom1_command = None
                self.create_atom2_command = None
                self.create_bond_command = None
                self.change_type_command = None
                self.drag_command = None

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
        self.delete_command = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.timer.start()
            self.point_d = self.widget.select(event.localPos())

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            if self.timer.run():
                self.timer.stop()
            if self.timer.time() > self.tol:
                if self.point_d is not None:
                    self.delete_command = DeleteAtomCommand(self.widget.drawing, self.point_d.parent, self.widget.ids)
                    self.delete_command.apply(self.point_d.create_command)
                    self.delete_command.appendStack(self.delete_command)
            else:
                self.point_d = None
            self.point_d = None
            self.delete_command = None


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
                        #dlp1 = point_class.Point(parent=dl_list, color=self.pc1, coord=self.pc1, rad=self.pc1)
                        #dlp2 = point_class.Point(parent=dl_list, color=point, coord=point, rad=point)
                        #self.widget.drawing.add_contact((self.pc1, dlp1), (point, dlp2))
                        #self.widget.newContact.emit()
                        dl_list = point_class.PointsList(parent=self.widget.dl_list)
                        com = ContactsCommand(self.widget, dl_list)
                        com.apply(self.pc1.create_command, point.create_command)
                        com.appendStack(com)
                        self.pc1.color = self.old_color
                        self.pc1 = None

            if event.buttons() == QtCore.Qt.RightButton:
                if self.pc1 is not None:
                    self.pc1.color = self.old_color
                    self.pc1 = None
                    self.old_color = None

    def detach(self):
        aDrawWidgetEvent.detach(self)
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
                    com = AngleCommand(self.widget)
                    com.apply(*[x.create_command for x in points], split=-1)
                    com.appendStack(com)
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
                    com = AngleCommand(self.widget)
                    com.apply(*[x.create_command for x in points], split=split)
                    com.appendStack(com)
                    self.p1 = []
                    self.p2 = []
                    self.old_colors = [[], []]
                    self.add = self.p1, self.old_colors[0]

    def detach(self):
        aDrawWidgetEvent.detach(self)
        colors = self.old_colors[0] + self.old_colors[1]
        points = self.p1 + self.p2
        for i, point in enumerate(points):
            point.color = colors[i]

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
                    #self.widget.drawing.addAvgDiff(*points, split=split)
                    #self.widget.newAvgDiff.emit()
                    com = AvgDiffCommand(self.widget)
                    com.apply(*[x.create_command for x in points], split=split)
                    com.appendStack(com)
                    self.p1 = []
                    self.p2 = []
                    self.old_colors = [[], []]
                    self.add = self.p1, self.old_colors[0]

    def detach(self):
        aDrawWidgetEvent.detach(self)
        colors = self.old_colors[0] + self.old_colors[1]
        points = self.p1 + self.p2
        for i, point in enumerate(points):
            point.color = colors[i]

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
                #self.widget.drawing.addMaxMeanPlaneDiff(*points)
                #self.widget.newMaxMeanPlaneDiff.emit()
                com = MaxMeanPlaneDiffCommand(self.widget)
                com.apply(*[x.create_command for x in points])
                com.appendStack(com)

                self.selected_points = []
                self.old_colors = []

    def detach(self):
        aDrawWidgetEvent.detach(self)
        for i, point in enumerate(self.selected_points):
            point.color = self.old_colors[i]

        self.selected_points = []
        self.old_colors = []


class ChangeCN(aDrawWidgetEvent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cn = [0, 14]
        self.cn_label = '0-14'

        self.change_cn_command = None

        self.point = None

    def assertEvent(self, event: QtCore.QEvent, widget):
        if event.type() == QtCore.QEvent.MouseButtonPress and event.buttons() == QtCore.Qt.LeftButton:
            self.timer.start()
            self.point = self.widget.select(event.localPos())

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            if self.timer.run():
                self.timer.stop()
            if self.timer.time() > self.tol:
                point2 = self.widget.select(event.localPos())
                if point2 is self.point and self.point is not None:
                    self.change_cn_command = ChangeCNCommand(self.widget.drawing)
                    self.change_cn_command.apply(self.point.create_command, self.cn, self.cn_label)
                    self.change_cn_command.appendStack(self.change_cn_command)
                self.point = None
                self.change_cn_command = None

    def setCN(self, cn):
        cn = [int(x) for x in cn]
        self.cn = cn
        self.cn_label = f'{cn[0]}-{cn[1]}'


class DrawerGL(QOpenGLWidget):

    newContact = Signal(name='newContact')
    newAngle = Signal(name='newAngle')
    newAvgDiff = Signal(name='newAvgDiff')
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

    def clearDraw(self):
        while len(Command.stack) != 0:
            Command.popStack()


from .ui import Drawer_model_ui


class DrawWidget(Drawer_model_ui.Ui_Dialog, QtWidgets.QDialog):

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

        self.contacts_view = TableView(remove=DeleteContactsCommand, parent=self, model=self.contacts_model, widget=self.openGl_drawer)

        self.angle_view = TableView(remove=DeleteAngleCommand, parent=self, model=self.angle_model, widget=self.openGl_drawer)

        self.avgDiff_view = TableView(remove=DeleteAvgDiffCommand, parent=self, model=self.avgDiff_model, widget=self.openGl_drawer)

        self.maxMeanPlaneDiff_view = TableView(remove=DeleteMaxMeanPlaneCommand, parent=self, model=self.maxMeanPlaneDiff_model, widget=self.openGl_drawer)

        self.openGl_drawer.newContact.connect(lambda: self.contacts_model.insertRows(self.contacts_model.rowCount(), 1))
        self.openGl_drawer.newAngle.connect(lambda: self.angle_model.insertRows(self.angle_model.rowCount(), 1))
        self.openGl_drawer.newAvgDiff.connect(lambda: self.avgDiff_model.insertRows(self.avgDiff_model.rowCount(), 1))
        self.openGl_drawer.newMaxMeanPlaneDiff.connect(lambda: self.maxMeanPlaneDiff_model.insertRows(self.maxMeanPlaneDiff_model.rowCount(), 1))

        '''self.openGl_drawer.removeContact.connect(lambda x: self.contacts_model.removeRows(x, 1))
        self.openGl_drawer.removeAngle.connect(lambda x: self.angle_model.removeRows(x, 1))
        self.openGl_drawer.removeAvgDiff.connect(lambda x: self.avgDiff_model.removeRows(x, 1))
        self.openGl_drawer.removeMaxMeanPlaneDiff.connect(lambda x: self.maxMeanPlaneDiff_model.removeRows(x, 1))'''

        self.verticalLayout_8.addWidget(self.contacts_view)
        self.verticalLayout_9.addWidget(self.angle_view)
        self.verticalLayout_9.addWidget(self.angle_view)
        self.verticalLayout_11.addWidget(self.avgDiff_view)
        self.verticalLayout_10.addWidget(self.maxMeanPlaneDiff_view)

        HighLight().attach(self.openGl_drawer)
        UndoRedo().attach(self.openGl_drawer)
        drag_event = Drag(exclusive=True, exclusive_group=1)
        drag_event.attach(self.openGl_drawer)
        self.draw_event = Draw(exclusive=True, exclusive_group=1)
        self.cn_event = ChangeCN(exclusive=True, exclusive_group=1)
        clear_event = Clear(exclusive=True, exclusive_group=1)
        contact_event = Contact(exclusive=True, exclusive_group=1)
        angle_event = Angle(exclusive=True, exclusive_group=1)
        avgDiff_event = AvgDiff(exclusive=True, exclusive_group=1)
        maxMeanPlaneDiff_event = MaxMeanPlaneDiff(exclusive=True, exclusive_group=1)

        self.pushButton_7.pressed.connect(self.exportTable)

        self.pushButton_15.pressed.connect(self.openGl_drawer.drawing.show_cn)
        self.pushButton_15.pressed.connect(self.openGl_drawer.update)
        self.pushButton_16.pressed.connect(self.openGl_drawer.drawing.show_at)
        self.pushButton_16.pressed.connect(self.openGl_drawer.update)

        self.pushButton_16.pressed.connect(self.exportTable)
        self.pushButton_10.pressed.connect(lambda: drag_event.attach(self.openGl_drawer))
        self.pushButton_14.pressed.connect(lambda: clear_event.attach(self.openGl_drawer))
        self.pushButton.pressed.connect(self.openGl_drawer.clearDraw)
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
        self.pushButton_17.pressed.connect(self.setCN)
        #self.pushButton_9.pressed.connect(lambda: self.checked_func.change_attr(func=1))
        #self.pushButton_8.pressed.connect(lambda: self.checked_func.change_attr(func=2))
        #self.dbSearchButton.pressed.connect(lambda: self.db_search(self.openGl_drawer.drawing.struct))

        #self.pushButton_6.pressed.connect(lambda: parent.substr_search(self.openGl_drawer.drawing.struct,
        #                                                               self.openGl_drawer.drawing.conditions))
        #self.pushButton_7.pressed.connect(self.openGl_drawer.save_template)
        self.gridLayout.addWidget(self.pushButton_6, 0, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_10, 1, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_14, 2, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton, 3, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_11, 7, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_17, 8, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_9, 9, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_8, 10, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_12, 11, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_13, 12, 0, 1, 2)

        #self.tableWidget_2.itemChanged.connect(lambda x: self.asd(x, tab=1))
        #self.tableWidget_3.itemChanged.connect(lambda x: self.asd(x, tab=2))

    def changeDraw(self, atom_type):
        if type(atom_type) is list and len(atom_type) == 1:
            atom_type = atom_type[0]
        self.draw_event.setType(atom_type)
        self.draw_event.attach(self.openGl_drawer)

    def changeCN(self, cn):
        self.cn_event.setCN(cn)
        self.cn_event.attach(self.openGl_drawer)

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

    def setCN(self):
        dialog = QtWidgets.QDialog(parent=self)
        dialog.setWindowTitle('C.N.')
        lineEd = QtWidgets.QLineEdit(parent=dialog)
        dialog.setLayout(QtWidgets.QVBoxLayout())
        dialog.layout().addWidget(lineEd)
        lineEd.returnPressed.connect(dialog.accept)
        func = lambda: self.getCN(lineEd)
        dialog.accepted.connect(func)
        dialog.show()

    def getCustomAtomType(self, lineEd: QtWidgets.QLineEdit):
        text = lineEd.text()
        if text == '':
            self.changeDraw('H')
        types = text.split(',')
        self.changeDraw(types)

    def getCN(self, lineEd: QtWidgets.QLineEdit):
        text = lineEd.text()
        if text == '':
            self.changeCN([0, 14])
            return
        cns = text.split('-')
        cn = [int(x) for x in cns]
        if cn[0] < 0:
            cn[0] = 0
        if cn[1] > 14:
            cn[1] = 14
        self.changeCN(cn)

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
            nodes = {}
            for point in points:
                if point not in mem:
                    if nodes.get(point, None) is not None:
                        node = nodes[point]
                    else:
                        node = Node(pack, str_atom=atoms[point])
                        nodes[point] = node
                    mem[point] = node
                    for point2 in bonds[point]:
                        if nodes.get(point2, None) is not None:
                            node2 = nodes[point2]
                        else:
                            node2 = Node(pack, str_atom=atoms[point2])
                            nodes[point2] = node2
                        node.addConnect(node2)
                else:
                    continue
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
        Db_viewer.DB_VIEWER.list_model.populate((draw, 'substructure', Db_viewer.DB_VIEWER.list_model._last_db_type, Db_viewer.DB_VIEWER.search_dialog.getDbString()))

    def exportTable(self):
        if self.organizer is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(parent=self, filter='*.csv')[0]
            if not filename:
                return
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
