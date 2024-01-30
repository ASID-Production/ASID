# Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
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


from PySide6 import QtWidgets, QtGui
from PySide6.QtCore import *
from .point_class import aPoint, Point, PointsList
from . import Observers
import typing
import numpy as np

SINGLE_OBSERVER = None


def setSingleObserver(single_observer):
    global SINGLE_OBSERVER
    SINGLE_OBSERVER = single_observer


class SimpleDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self, parent=None):
        super().__init__(parent=parent)

    def createEditor(self, parent: QtWidgets.QWidget, option: 'QStyleOptionViewItem', index: QModelIndex) -> QtWidgets.QWidget:
        if 'color' == index.data(99)[0]:
            color = QtWidgets.QColorDialog(parent=parent)
            color.setOption(QtWidgets.QColorDialog.ShowAlphaChannel, True)
            return color
        else:
            return super().createEditor(parent, option, index)

    def setEditorData(self, editor: QtWidgets.QWidget, index: QModelIndex) -> None:
        if isinstance(editor, QtWidgets.QLineEdit):
            editor.setText(index.data(Qt.ItemDataRole.DisplayRole))
        elif isinstance(editor, QtWidgets.QColorDialog):
            color = index.data(99)
            color = (color[0], color[1][0][1])
            if issubclass(type(color[1]), aPoint):
                color = (color[0], color[1].__getattribute__(color[0]))
            if color[1] is None:
                color = [0,0,0,255]
            else:
                color = [int(x*255) for x in color[1]]
            editor.setCurrentColor(QtGui.QColor(*color))

    def setModelData(self, editor: QtWidgets.QWidget, model: 'QtPointsPropertyModel', index: QModelIndex) -> None:
        if isinstance(editor, QtWidgets.QLineEdit):
            line = editor.text()
            line = line.split(':')
            if len(line) == 1:
                value = line[0]
            else:
                property, value = line
            value = value.replace(' ', '')
            if value != '' and value[0] == '[' and value[-1] == ']':
                value = value[1:-1].split(',')
                value = np.array([float(x) for x in value], dtype=np.float32)
            elif value == 'None' or value == '':
                value = None
            elif value[0] == '$':
                seq = int(value[1:])
                value = model.getBySeqData(seq)
            else:
                try:
                    for char in value:
                        if char.isdigit() or char == '.' or char == '-':
                            continue
                        else:
                            model.setData(index, value, role=99)
                            return
                    value = float(value)
                except ValueError:
                    pass
            model.setData(index, value, role=99)
        elif isinstance(editor, QtWidgets.QColorDialog):
            color = editor.currentColor()
            value = np.array(color.getRgbF(), dtype=np.float32)
            model.setData(index, value, role=99)


class UniformListModel(QAbstractListModel):

    def __init__(self, parent=None, data=None):
        super().__init__(parent)
        self.rows = 0
        self._root = None
        self.setModelData(data)

    def setModelData(self, data):
        self.beginResetModel()
        if data is None:
            if self.rowCount() > 0:
                self.removeRows(0, self.rowCount())
            self.endResetModel()
            return
        elif len(data.getInfo()) < self.rowCount():
            self.romoveRows(len(data.getInfo()), self.rowCount() - len(data.getInfo()))
        elif len(data.getInfo()) > self.rowCount():
            self.insertRows(self.rowCount(), len(data.getInfo()) - self.rowCount())
        self._root = data
        self.endResetModel()

    def rowCount(self, parent=QModelIndex(), *args, **kwargs) -> int:
        return self.rows

    def removeRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginRemoveRows(parent, row, row+count-1)
        self.rows -= count
        self.endRemoveRows()
        return True

    def insertRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginInsertRows(parent, row, row+count-1)
        self.rows += count
        self.endInsertRows()
        return True

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = ...):
        if role != Qt.ItemDataRole.DisplayRole and role != Qt.ItemDataRole.ToolTipRole:
            return None
        if self.selected is None:
            return 'Data'
        else:
            return self.selected.internalPointer().__str__()

    def getData(self, index):
        properties = list(self._root.getInfo().keys())
        property = properties[index.row()]
        return property, self._root.__getattr__(property)

    def data(self, index: QModelIndex, role: int = ...):
        if self._root is None:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            properties = list(self._root.getInfo().keys())
            try:
                property = properties[index.row()]
            except IndexError:
                return None
            line = self._root.getLine(property)
            return line
        elif role == 99:
            property = list(self._root.getInfo().keys())[index.row()]
            return property, self._root.__getattribute__(property)

    def setData(self, index, value, role=None):
        property = list(self._root.getInfo().keys())[index.row()]
        self._root.__setattr__(property, value)

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable


class QtPointsPropertyModel(QAbstractListModel):

    def __init__(self, parent=None, data=None):
        self._root = data
        self.selected = []
        self.selected_len = 0
        self.rows = 0
        self.props = {}
        super().__init__(parent)

    def addProperty(self, property):
        self.beginResetModel()
        if self.selected:
            for selected in self.selected:
                selected.internalPointer().addProperty(property, None)
        self.endResetModel()

    def removeProperty(self, index=QModelIndex()):
        if not self.selected or not index.isValid():
            return
        else:
            self.beginResetModel()
            for selected in self.selected:
                selected.internalPointer().deleteProperty(self.data(index, role=99)[0])
                self.setSelected(self.selected)
            self.endResetModel()

    def setModelData(self, data):
        self._root = data

    def setSelected(self, selected, command: QItemSelectionModel.SelectionFlags):
        self.beginResetModel()
        if isinstance(selected, QModelIndex):
            if not selected.isValid():
                if self.selected_len != 0:
                    self.selected = []
                    self.removeRows(0, self.rows, QModelIndex())
                    self.selected_len = 0
                    self.props = {}
                    return
                else:
                    return
            else:
                selected = [selected]
        elif isinstance(selected, QItemSelection):
            selected = selected.indexes()
        for sel in selected:
            if (command & QItemSelectionModel.Clear) == QItemSelectionModel.Clear:
                if self.selected_len != 0:
                    self.selected = []
                    self.removeRows(0, self.rows, QModelIndex())
                    self.selected_len = 0
                    self.props = {}
            if (command & QItemSelectionModel.Deselect) == QItemSelectionModel.Deselect:
                try:
                    self.selected.remove(sel)
                    props = sel.internalPointer().getProperties().keys()
                    for key in props:
                        self.props[key] -= 1
                        if self.props[key] == 0:
                            self.props.pop(key)
                            self.selected_len -= 1
                            self.removeRows(self.selected_len, self.rows - self.selected_len, QModelIndex())
                except ValueError:
                    pass
            if (command & QItemSelectionModel.Select) == QItemSelectionModel.Select:
                if sel not in self.selected:

                    props = sel.internalPointer().getProperties().keys()
                    old_sel_len = self.selected_len
                    for key in props:
                        num = self.props.get(key, 0)
                        if num == 0:
                            self.selected_len += 1
                        self.props[key] = num + 1

                    if old_sel_len != self.selected_len:
                        self.insertRows(self.rows, self.selected_len - self.rows, QModelIndex())
                    self.selected.append(sel)

        self.endResetModel()

    def rowCount(self, parent=QModelIndex(), *args, **kwargs) -> int:
        return self.rows

    def removeRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginRemoveRows(parent, row, row+count-1)
        self.rows -= count
        self.endRemoveRows()
        return True

    def insertRows(self, row, count, parent=QModelIndex(), *args, **kwargs) -> bool:
        self.beginInsertRows(parent, row, row+count-1)
        self.rows += count
        self.endInsertRows()
        return True

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = ...):
        if role != Qt.ItemDataRole.DisplayRole and role != Qt.ItemDataRole.ToolTipRole:
            return None
        if self.selected is None:
            return 'Data'
        else:
            return 'Data'

    def getBySeqData(self, seq):
        if self._root is None:
            return None
        else:
            return self._root.getBySeq(seq)

    def data(self, index: QModelIndex, role: int = ...):
        if not self.selected:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            try:
                property = list(self.props.keys())[index.row()]
            except IndexError:
                return None
            line_p = None
            for selected in self.selected:
                line = selected.internalPointer().getLine(property)
                if line_p is not None:
                    if line == line_p:
                        line_p = line
                        continue
                    else:
                        line = f'{property}: '
                        break
                line_p = line
            return line
        elif role == 99:
            try:
                property = list(self.props.keys())[index.row()]
            except IndexError:
                return None
            value = [(x, x.internalPointer().getProperties()[property]) for x in self.selected]
            return property, value

    def setData(self, index, value, role=None):
        if role == Qt.ItemDataRole.DisplayRole:
            line = self.data(index, role=role)
            property = line.split(':')[0]
            for selected in self.selected:
                if property in selected.internalPointer().getProperties().keys():
                    pass
                else:
                    selected.internalPointer().addProperty(property, None)
                selected.internalPointer().__setattr__(property, value)
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
        return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable


class QtPointsTreeModel(QAbstractItemModel):

    def __init__(self, parent=None, data=None):
        if data is None:
            data = PointsList()
        self._root = data
        super().__init__(parent=parent)
        a = 0

    def setModelData(self, data):
        self._root = data

    def getRoot(self):
        return self._root

    def update(self):
        self.beginResetModel()
        self.endResetModel()

    def attachObserver(self, index: QModelIndex, observer_name):
        if index.isValid() and index.internalPointer() is not None:
            observer_cls = Observers.observers.get(observer_name, None)
            if observer_cls is None:
                return
            if SINGLE_OBSERVER is not None:
                observer = SINGLE_OBSERVER.getObserver(observer_cls)
                if observer not in index.internalPointer().observers:
                    index.internalPointer().attach(observer)
                return

    def detachObserver(self, index: QModelIndex, observer_name):
        if index.isValid() and index.internalPointer() is not None:
            observer_cls = Observers.observers.get(observer_name, None)
            if observer_cls is None:
                return
            item = index.internalPointer()
            for observer in item.observers:
                if isinstance(observer, observer_cls):
                    item.detach(observer)
                    break
            return

    def index(self, row, column, parent=QModelIndex(), *args, **kwargs) -> QModelIndex:
        if kwargs.get('by_point', None) is not None:
            point = kwargs['by_point']
            child = kwargs['by_point']
            path = []
            parent = point.parent
            while parent is not self._root:
                ind = parent.children.index(child)
                child = parent
                path.append((parent, ind))
                parent = parent.parent
            path.append((QModelIndex(), self._root.children.index(child)))
            step = path.pop()
            index = self.index(step[1], column, parent=step[0])
            while path:
                step = path.pop()
                index = self.index(step[1], column, parent=index)
            return index
        if not self.hasIndex(row, column, parent):
            return QModelIndex()
        if not parent.isValid():
            parentItem = self._root
        else:
            parentItem = parent.internalPointer()
        childItem = parentItem.children[row]
        if childItem:
            index = self.createIndex(row, column, childItem)
            return index
        else:
            return QModelIndex()

    def parent(self, child: QModelIndex) -> QModelIndex:
        if not child.isValid():
            return QModelIndex()
        childItem = child.internalPointer()
        parentItem = childItem.parent
        if parentItem == self._root:
            return QModelIndex()
        if parentItem.parent is None:
            return self.createIndex(0, 0, parentItem)
        else:
            return self.createIndex(parentItem.parent.children.index(parentItem), 0, parentItem)

    def rowCount(self, parent=QModelIndex(), *args, **kwargs):
        if parent.column() > 0:
            return 0
        if parent.isValid():
            return len(parent.internalPointer().children)
        return len(self._root.children)

    def columnCount(self, parent: QModelIndex = QModelIndex(), *args, **kwargs) -> int:
        if parent.isValid():
            return 1
        return 1

    def data(self, index: QModelIndex, role: int = ...):
        if not index.isValid():
            return None
        if role != Qt.ItemDataRole.DisplayRole and role != Qt.ItemDataRole.ToolTipRole:
            return None
        item = index.internalPointer()
        if role == Qt.ItemDataRole.DisplayRole:
            if item.name is None:
                return item.__str__()
            return item.name
        return item.__str__()

    def setItemData(self, index: QModelIndex, roles, *args, **kwargs) -> bool:
        data = roles.get(Qt.ItemDataRole.DisplayRole)
        if 'Point' not in data:
            return False
        item = index.internalPointer()
        data = int(data.split(' ')[1])
        origin = self._root.getBySeq(data)
        if origin.parent == item.parent:
            item.parent.changeOrder(origin, item.parent.children.index(item))
        return True

    def setData(self, index: QModelIndex, value, role: int = ...) -> bool:
        if role == 99:
            property, value = value
            point = index.internalPointer()
            point.__setattr__(property, value)
            return True
        else:
            return super().setData(index, value, role)

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        if type(index.internalPointer()).__name__ == 'PointsList':
            return (Qt.ItemFlag.ItemIsEnabled|Qt.ItemFlag.ItemIsSelectable)
        return (Qt.ItemFlag.ItemIsEnabled|Qt.ItemFlag.ItemIsSelectable|Qt.ItemFlag.ItemIsDragEnabled)

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = ...):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return 'Points'
        return None

    def moveRows(self, sourceParent: QModelIndex, sourceRow: int, count: int, destinationParent: QModelIndex, destinationChild: int) -> bool:
        return super().moveRows(sourceParent, sourceRow, count, destinationParent, destinationChild)

    def dropMimeData(self, data: 'QMimeData', action: Qt.DropAction, row: int, column: int, parent: QModelIndex) -> bool:
        source = [int(seq) for seq in data.data('source').data().decode('utf-8').split(',')]
        dest = self.index(row, column, parent)
        parent_p = self._root.getBySeq(source[0]).parent
        if len(source) > 1:
            for seq in source[1:]:
                if parent_p is not self._root.getBySeq(seq).parent:
                    return False
        if parent.internalPointer() is not parent_p:
            return False
        for seq in source:
            if row == parent_p.children.index(parent_p.getBySeq(seq)) or row - 1 == parent_p.children.index(parent_p.getBySeq(seq)):
                continue
            self.beginMoveRows(parent, parent_p.children.index(parent_p.getBySeq(seq)), parent_p.children.index(parent_p.getBySeq(seq)),
                               parent, row)
            if row > parent_p.children.index(parent_p.getBySeq(seq)):
                parent_p.changeOrder(self._root.getBySeq(seq), row - 1)
            else:
                parent_p.changeOrder(self._root.getBySeq(seq), row)
            self.endMoveRows()
        return True

    def mimeTypes(self) -> typing.List[str]:
        return ['application/x-qabstractitemmodeldatalist', 'source']

    def mimeData(self, indexes: typing.Iterable[QModelIndex]) -> 'QMimeData':
        if len(indexes) == 0:
            return None
        types = self.mimeTypes()
        if len(types) == 0:
            return None
        data = QMimeData()
        format = types[0]
        encoded = QByteArray()
        stream = QDataStream(encoded, QIODevice.WriteOnly)
        self.encodeData(indexes, stream)
        data.setData(format, encoded)
        if 'source' in types:
            format = 'source'
            line = ','.join([str(index.internalPointer().seq) for index in indexes]).encode('utf-8')
            encoded = QByteArray(line)
            data.setData(format, encoded)
        return data

    def insertRow(self, row, parent=QModelIndex(), *args, **kwargs):
        type = kwargs.get('type', None)
        if type is None:
            self.beginResetModel()
            ret = super().insertRow(row, parent=parent)
            self.endResetModel()
            return ret
        if parent.internalPointer() is None:
            parent_item = self._root
        else:
            parent_item = parent.internalPointer()
        if isinstance(parent_item, Point):
            return False
        if type == 'point':
            Point(parent=parent_item)
        elif type == 'list':
            PointsList(parent=parent_item)
        else:
            return False
        row = self.rowCount(parent=parent)
        self.beginInsertRows(parent, row, row+1)
        super().insertRow(row, parent=parent)
        self.endInsertRows()

    def removeRow(self, row, parent=QModelIndex(), *args, **kwargs):
        if row == -1:
            return False
        if not parent.isValid():
            item = self.index(row, 0)
        else:
            item = self.index(row, 0, parent=parent)
        item = item.internalPointer()
        self.beginRemoveRows(parent, row, row+1)
        item.destroy()
        self.endRemoveRows()
        return True


class TreeView(QtWidgets.QTreeView):

    def __init__(self, parent=None):
        super(TreeView, self).__init__(parent)
        self.startDragPosition = None
        self.drag = None
        self.installEventFilter(self)
        self.setDragEnabled(True)
        self.setDropIndicatorShown(True)
        self.setAcceptDrops(True)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)

    def eventFilter(self, object: QObject, event: QEvent) -> bool:
        return super().eventFilter(object, event)

    def mousePressEvent(self, e: QtGui.QMouseEvent) -> None:
        item = self.indexAt(e.pos())
        super().mousePressEvent(e)
        if item.row() == -1 and item.column() == -1:
            self.clearSelection()
            self.selectionModel().setCurrentIndex(QModelIndex(), QItemSelectionModel.Select)

    def dragEnterEvent(self, event: QtGui.QDragEnterEvent) -> None:
        super().dragEnterEvent(event)
        if event.mimeData().hasText():
            event.acceptProposedAction()
            self.setState(self.DraggingState)

    def dragLeaveEvent(self, event: QtGui.QDragLeaveEvent) -> None:
        super().dragLeaveEvent(event)

    def dragMoveEvent(self, event: QtGui.QDragMoveEvent) -> None:
        super().dragMoveEvent(event)
        event.acceptProposedAction()

    def dropEvent(self, event: QtGui.QDropEvent) -> None:
        super().dropEvent(event)

    def showContextMenu(self, pos: QPoint):
        selected = self.currentIndex()

        menu = QtWidgets.QMenu('Context Menu', self)

        if selected.internalPointer() is not None:
            def save_state(func, *args):
                return lambda: func(*args)
            menu_attach = QtWidgets.QMenu('Attach representation')
            menu_detach = QtWidgets.QMenu('Detach representation')
            menu.addMenu(menu_attach)
            menu.addMenu(menu_detach)
            actions = []
            for i, observer_key in enumerate(Observers.observers):
                f = False
                for observer in selected.internalPointer().observers:
                    if type(observer) is Observers.observers[observer_key]:
                        action = QtGui.QAction(f'Detach {observer_key}')
                        func = save_state(self.model().detachObserver, selected, observer_key)
                        action.triggered.connect(func)
                        menu_detach.addAction(action)
                        actions.append(action)
                        f = True
                if f:
                    continue
                else:
                    action = QtGui.QAction(f'Attach {observer_key}')
                    func = save_state(self.model().attachObserver, selected, observer_key)
                    action.triggered.connect(func)
                    menu_attach.addAction(action)
                    actions.append(action)

        def delete(row, parent):
            self.model().removeRow(row=selected.row(), parent=selected.parent())
            self.selectionModel().clearSelection()
        action1 = QtGui.QAction('Add list', self)
        action2 = QtGui.QAction('Add point', self)
        action3 = QtGui.QAction('Delete', self)
        action1.triggered.connect(lambda: self.model().insertRow(row=-1, parent=selected, type='list'))
        action2.triggered.connect(lambda: self.model().insertRow(row=-1, parent=selected, type='point'))
        action3.triggered.connect(lambda: delete(row=selected.row(), parent=selected.parent()))
        menu.addAction(action1)
        menu.addAction(action2)
        menu.addAction(action3)
        menu.exec(self.mapToGlobal(pos))
        return

    def setSelection(self, rect: QRect, command) -> None:
        super().setSelection(rect, command)


class ListView(QtWidgets.QListView):

    def __init__(self, parent=None):
        super(ListView, self).__init__(parent)
        self.setItemDelegate(SimpleDelegate(self))
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)

    def mousePressEvent(self, e: QtGui.QMouseEvent) -> None:
        item = self.indexAt(e.pos())
        super().mousePressEvent(e)
        if item.row() == -1 and item.column() == -1:
            self.clearSelection()
            self.selectionModel().setCurrentIndex(QModelIndex(), QItemSelectionModel.Select)

    def addProperty(self):
        def addProp(self, prop):
            self.model().addProperty(prop)
            return
        dialog = QtWidgets.QDialog(parent=self)
        layout = QtWidgets.QHBoxLayout(dialog)
        layout.addWidget(QtWidgets.QLabel('Property name:'))
        lineEdit = QtWidgets.QLineEdit()
        lineEdit.returnPressed.connect(dialog.accept)
        dialog.accepted.connect(lambda: addProp(self, lineEdit.text()))
        layout.addWidget(lineEdit)
        dialog.open()

    def showContextMenu(self, pos: QPoint):
        menu = QtWidgets.QMenu('Context Menu', self)
        action1 = QtGui.QAction('Add property', self)
        action2 = QtGui.QAction('Delete', self)
        selected = self.currentIndex()

        action1.triggered.connect(self.addProperty)
        action2.triggered.connect(lambda: self.model().removeProperty(index=selected))
        menu.addAction(action1)
        menu.addAction(action2)
        menu.exec(self.mapToGlobal(pos))
        return


class SelectionModel(QItemSelectionModel):

    newSelection = Signal(tuple, name='newSelection') # Tuple[QModelIndex | QItemSelection, QItemSelectionModel.SelectionFlags]
    
    def __init__(self, model=None):
        super(SelectionModel, self).__init__(model=model)

    def select(self, index, command: QItemSelectionModel.SelectionFlags):
        if isinstance(self.model(), QtPointsTreeModel):
            if isinstance(index, QItemSelection):
                index_l = index.indexes()
            elif isinstance(index, QModelIndex):
                index_l = [index]
            for index_m in index_l:
                if index_m.internalPointer() is not None:
                    if index_m.internalPointer().pick is None:
                        index_m.internalPointer().addProperty('pick', 0.0)
                    if (command & QItemSelectionModel.Clear) == QItemSelectionModel.Clear:
                        for ind in self.selectedIndexes():
                            self.model().setData(ind, ('pick', 0.0), role=99)
                            pass
                    if (command & QItemSelectionModel.Deselect) == QItemSelectionModel.Deselect:
                        self.model().setData(index_m, ('pick', 0.0), role=99)
                    if (command & QItemSelectionModel.Select) == QItemSelectionModel.Select:
                        self.model().setData(index_m, ('pick', 1.0), role=99)
                    self.newSelection.emit((index, command))
        return super().select(index, command)

    def clearSelection(self) -> None:
        if isinstance(self.model(), QtPointsTreeModel):
            self.newSelection.emit((QModelIndex(), QItemSelectionModel.Clear))
            for index in self.selectedIndexes():
                self.model().setData(index, ('pick', 0.0), role=99)
        super().clearSelection()