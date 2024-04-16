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
import os.path

from PySide6 import QtWidgets
from PySide6.QtCore import *
from . import Db_bindings
import typing
import numpy as np
from ... import point_class
from ..ChemPack import PALETTE, MOLECULE_SYSTEMS, TREE_MODEL
from . import MoleculeClass

DB_VIEWER = None


class StructuresListModel(QAbstractListModel):

    _display_tags = ['id', 'refcode', 'CCDC_number', 'formula', 'temperature']

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self._data = []
        self._rows = 0
        self._display_tag = 'refcode'

    def populate(self, request):
        request, search_type = request
        self.beginResetModel()
        self._data = Db_bindings.search(request, search_type)
        self._rows = len(self._data)
        self.endResetModel()

    def setDisplayTag(self, tag):
        if tag in self._display_tags:
            self._display_tag = tag

    def displayTags(self):
        return self._display_tags

    def rowCount(self, parent: QModelIndex = ...) -> int:
        return self._rows

    def columnCount(self, parent: QModelIndex = ...) -> int:
        return 1

    def data(self, index: QModelIndex, role: int = ...) -> typing.Any:
        if not index.isValid():
            return None
        if self._data is None:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return self._data[index.row()][self._display_tag]
        if role == 99:
            return self._data[index.row()]


class InfoTableModel(QAbstractTableModel):

    info_fields = ['id', 'refcode', 'CCDC_number', 'cell', 'reduced_cells', 'compound_name', 'formula', 'authors', 'publication', 'experiment', 'refinement_info', 'coordinates', 'crystal_info']

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self._selected = None
        self._fields = []
        self._rows = 0
        self._spans = []
        self._columns = 2

    def rowCount(self, parent: QModelIndex = ...) -> int:
        return self._rows

    def columnCount(self, parent: QModelIndex = ...) -> int:
        return self._columns

    def span(self, index: QModelIndex) -> 'QSize':
        if not index.isValid():
            return QSize(1, 1)
        if index.row() in self._spans:
            return QSize(1, 2)
        else:
            return QSize(1, 1)

    def selected(self):
        return self._selected

    def setSelected(self, selected_id):
        if selected_id is None:
            self.beginResetModel()
            self._selected = None
            self._fields = []
            self._rows = 0
            self.endResetModel()

        def rec_fill(obj):
            if type(obj) is dict:
                for key in obj:
                    if type(obj[key]) is dict or type(obj[key]) is list:
                        self._fields.append([str(key), ''])
                        self._spans.append(len(self._fields)-1)
                    else:
                        if obj[key] is not None:
                            self._fields.append([str(key), str(obj[key])])
                    rec_fill(obj[key])
            elif type(obj) is list:
                for item in obj:
                    rec_fill(item)

        self.beginResetModel()
        self._fields = []
        self._spans = []
        self._selected = Db_bindings.get_full_info(selected_id)
        rec_fill(self._selected)
        self._rows = len(self._fields)
        self.endResetModel()

    def data(self, index: QModelIndex, role: int = ...) -> typing.Any:
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return self._fields[index.row()][index.column()]

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = ...) -> typing.Any:
        return None


from .ui import search_dialog


class SearchDialog(search_dialog.Ui_Dialog, QtWidgets.QDialog):

    SEARCH_TYPES = Db_bindings.search_types

    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle('Search parameters')
        self.search_type = SearchDialog.SEARCH_TYPES[0]
        self.text = self.lineEdit
        self.checkBox.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[0], self.lineEdit))
        self.checkBox_8.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[1], self.lineEdit_8))
        self.checkBox_9.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[2], self.lineEdit_9))
        self.checkBox_10.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[3], self.lineEdit_10))
        self.checkBox_11.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[4], self.lineEdit_11))
        self.checkBox_12.stateChanged.connect(lambda: self.setSearchState(SearchDialog.SEARCH_TYPES[5], self.lineEdit_12))

    def setSearchState(self, type, text):
        self.search_type = type
        self.text = text

    def getText(self):
        return self.text.text()

    def getSearchType(self):
        return self.search_type


from .ui import Settings_dialog


class DbSettings(Settings_dialog.Ui_Dialog, QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle('Db Settings')
        self.pushButton.pressed.connect(lambda: self.connectServer(self.lineEdit.text(), self.lineEdit_2.text()))
        self.pushButton_2.pressed.connect(lambda: self.startServer(self.lineEdit_2.text()))
        self.pushButton_3.pressed.connect(lambda: self.login(self.lineEdit_3.text(), self.lineEdit_4.text()))

    def connectServer(self, address, port):
        if len(address.split('.')) == 4:
            try:
                int(port)
            except ValueError:
                return
            Db_bindings.SESSION.connectServer(address, port)

    def startServer(self, port):
        try:
            int(port)
        except ValueError:
            return
        Db_bindings.SESSION.startServer(port)

    def login(self, name, passwd):
        if name and passwd:
            Db_bindings.SESSION.login(name, passwd)


from .ui import base_search_window


class DbWindow(base_search_window.Ui_Dialog, QtWidgets.QDialog):

    def __init__(self, parent=None, ):
        super().__init__()
        self.setupUi(self)

        self.setWindowTitle('DB Search')

        self.list_model = StructuresListModel()
        self.listView: QtWidgets.QListView
        self.listView.setModel(self.list_model)
        self.listView.setSelectionModel(QItemSelectionModel())
        self.listView.selectionModel().currentChanged.connect(self.newSelection)

        self.search_dialog = SearchDialog(parent=self)
        self.db_settings = DbSettings(parent=self)

        self.table_model = InfoTableModel()
        self.tableView: QtWidgets.QTableView
        self.tableView.setModel(self.table_model)
        self.table_model.modelReset.connect(self.setSpan)
        self.tableView.verticalHeader().hide()
        self.tableView.horizontalHeader().hide()
        self.tableView.resizeColumnsToContents()
        self.tableView.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        self.tableView.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)

        self.pushButton.pressed.connect(self.search_dialog.show)
        self.search_dialog.pushButton.pressed.connect(lambda: self.list_model.populate((self.search_dialog.getText(), self.search_dialog.getSearchType())))
        self.pushButton_2.pressed.connect(self.loadStruct)

        self.pushButton_3.pressed.connect(self.saveCif)

        self.pushButton_4.pressed.connect(self.uploadFile)

        self.pushButton_5.pressed.connect(self.db_settings.show)

    def setSpan(self):
        for i in range(self.table_model.rowCount()):
            ind = self.table_model.index(i, 0)
            span = self.table_model.span(ind)
            self.tableView.setSpan(ind.row(), ind.column(), span.width(), span.height())

    def newSelection(self, current, previous):
        if current.isValid():
            id = self.list_model.data(current, 99)['id']
            self.table_model.setSelected(id)
        else:
            self.table_model.setSelected(None)

    def saveCif(self):
        data = self.table_model.selected()
        if data is None:
            return
        else:
            id = data['id']
            cif = Db_bindings.getCif(id)
            filename = QtWidgets.QFileDialog.getSaveFileName()[0]
            if filename:
                out = open(filename, 'wb')
                out.write(cif)
                out.close()

    def uploadFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(filter='*.cif *.xml')[0]
        if filename:
            Db_bindings.uploadFile(filename)

    def loadStruct(self):
        from ..ChemPack import pars

        if TREE_MODEL is None:
            return

        data = self.table_model.selected()
        if data is None:
            return
        id = data['id']
        cif = Db_bindings.getCif(id)
        if not os.path.exists('temp/'):
            os.mkdir('temp')
        filename = f'temp/{id}.cif'
        out = open(filename, 'wb')
        out.write(cif)
        out.close()
        pars(filename, True, TREE_MODEL.getRoot())
        return

def show():
    global DB_VIEWER
    if DB_VIEWER is None:
        DB_VIEWER = DbWindow()
    DB_VIEWER.show()
