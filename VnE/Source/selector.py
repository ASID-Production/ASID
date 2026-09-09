import os
from PySide6.QtWidgets import (
    QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QPushButton,
    QLabel, QMessageBox,
    QGroupBox, QSizePolicy, QSplitter, QTreeView, QAbstractItemView
)
from PySide6.QtCore import Qt
from PySide6 import QtCore
from PySide6.QtGui import QStandardItemModel, QStandardItem

from .point_class import PointsList, Point
from .Observers import observers
observers = tuple(observers.keys())

try:
    from pyqcodeeditor.QCodeEditor import QCodeEditor
    from pyqcodeeditor.highlighters import QPythonHighlighter
    from pyqcodeeditor.completers import QPythonCompleter

    class PyQCodeEditor(QCodeEditor):
        def set_code(self, code):
            self.setPlainText(code)

        def get_code(self):
            return self.toPlainText()
        def setup(self):
            self.setCompleter(QPythonCompleter())
            self.setHighlighter(QPythonHighlighter())

except ImportError:
    from PySide6.QtWidgets import QPlainTextEdit
    class PyQCodeEditor(QPlainTextEdit):
        def set_code(self, code):
            self.setPlainText(code)
        def get_code(self):
            return self.toPlainText()
        def setup(self):
            return


class Executor:

    def __init__(self, model, root, selection_model=None, parent=None):
        self.selection = []
        self.root = root
        self.model = model
        self.parent = parent
        self.selection_model = selection_model

    def attach(self, point, observer):
        if observer in observers:
            index = self.model.index(by_point=point)
            self.model.attachObserver(index, observer)

    def detach(self, point, observer):
        if observer in observers:
            index = self.model.index(by_point=point)
            self.model.detachObserver(index, observer)

    def addProp(self, key, value, point):
        point.addProperty(key, value)

    def select(self, point):
        if not point:
            return
        if isinstance(point, list):
            sel = QtCore.QItemSelection()
            for i in point:
                index = self.model.index(by_point=i)
                sel.select(index, index)
            ind = sel
        else:
            ind = self.model.index(by_point=point)
        self.selection_model.select(ind, QtCore.QItemSelectionModel.Select)

    def reg(self, reg, string):
        import re
        ret = re.fullmatch(reg, string)
        if ret:
            return True
        else:
            return False

    def exec(self, code):
        attach = self.attach
        detach = self.detach
        addProp = self.addProp
        select = self.select
        regex = self.reg
        selected = []
        for ind in self.selection:
            if ind.isValid() and ind.internalPointer() is not None and ind.internalPointer()._checked:
                if  ind.internalPointer() not in selected:
                    selected.append(ind.internalPointer())
            elif ind.isValid() and ind.internalPointer() is not None:
                self.selection.remove(ind)
        root = self.root
        try:
            exec(code)
        except Exception as e:
            if self.parent:
                QMessageBox.critical(self.parent, "Error", f"Error executing macro:\n{e}")


class MacroDockWidget(QDockWidget):
    def __init__(self, macro_folder: str, parent=None, tree_model=None, selection_model=None):
        super().__init__("Macro", parent)
        self.macro_folder = macro_folder
        if tree_model is None:
            self.tree_model = QStandardItemModel()
            self.tree_model.setHorizontalHeaderLabels(["Points"])
        else:
            self.tree_model = tree_model
            self.tree_model.dataChanged.connect(self.on_check)
        self.selection_model = selection_model
        self._create_ui()
        self.load_macros()
        self.load_macros()
        self.on_macro_selected(0)
        self.executor = Executor(self.tree_model, self.tree_model.getRoot(), self.selection_model, self)

    def _create_ui(self):
        central = QWidget()
        self.setWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(6, 6, 6, 6)

        macro_layout = QHBoxLayout()
        macro_layout.addWidget(QLabel("Macro:"))
        self.macro_combo = QComboBox()
        self.macro_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.setToolTip("Refresh macros list")
        macro_layout.addWidget(self.macro_combo)
        macro_layout.addWidget(self.refresh_btn)
        main_layout.addLayout(macro_layout)

        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.setChildrenCollapsible(False)

        group_box = QGroupBox("Points selection")
        group_layout = QVBoxLayout(group_box)
        group_layout.setContentsMargins(4, 4, 4, 4)

        self.tree_view = QTreeView()
        self.tree_view.setModel(self.tree_model)
        self.tree_view.setHeaderHidden(True)
        self.tree_view.setAlternatingRowColors(True)
        self.tree_view.setDragEnabled(False)
        self.tree_view.setAcceptDrops(False)
        self.tree_view.setDragDropMode(QAbstractItemView.NoDragDrop)
        group_layout.addWidget(self.tree_view)
        self.splitter.addWidget(group_box)
        self.splitter.setObjectName('codeSplitter')

        editor_widget = QWidget()
        editor_layout = QVBoxLayout(editor_widget)
        editor_layout.setContentsMargins(0, 0, 0, 0)
        editor_layout.addWidget(QLabel("Macro editor:"))
        self.editor = PyQCodeEditor()
        self.editor.setObjectName('macroEditor')
        self.editor.setup()
        editor_layout.addWidget(self.editor)
        self.splitter.addWidget(editor_widget)

        self.splitter.setStretchFactor(0, 2)  # дерево 40%
        self.splitter.setStretchFactor(1, 3)  # редактор 60%
        main_layout.addWidget(self.splitter, stretch=1)

        button_layout = QHBoxLayout()
        self.save_btn = QPushButton("Save macro")
        self.apply_btn = QPushButton("Execute macro")
        button_layout.addStretch()
        button_layout.addWidget(self.save_btn)
        button_layout.addWidget(self.apply_btn)
        main_layout.addLayout(button_layout)

        self.macro_combo.currentIndexChanged.connect(self.on_macro_selected)
        self.refresh_btn.clicked.connect(self.load_macros)
        self.save_btn.clicked.connect(self.save_macro)
        self.apply_btn.clicked.connect(self.apply_macro)

        settings = QtCore.QSettings('ASID', 'LRSI')
        state = settings.value("macroDock/splitter")
        if state:
            self.splitter.restoreState(state)

    def load_macros(self):
        self.macro_combo.blockSignals(True)
        self.macro_combo.clear()
        if os.path.isdir(self.macro_folder):
            for fname in sorted(os.listdir(self.macro_folder)):
                if fname.endswith(".py"):
                    self.macro_combo.addItem(fname)
        self.macro_combo.blockSignals(False)
        self.on_macro_selected(self.macro_combo.currentIndex())

    def on_macro_selected(self, index):
        if index < 0:
            return
        fname = self.macro_combo.itemText(index)
        path = os.path.join(self.macro_folder, fname)
        try:
            with open(path, "r") as f:
                code = f.read()
            self.editor.set_code(code)
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Can't load macro:\n{e}")

    def save_macro(self):
        fname = self.macro_combo.currentText()
        if not fname:
            QMessageBox.warning(self, "Warning", "Choose existing macro")
            return
        if not fname.endswith(".py"):
            fname += ".py"
        path = os.path.join(self.macro_folder, fname)
        code = self.editor.get_code()
        try:
            with open(path, "w") as f:
                f.write(code)
            '''self.macro_combo.blockSignals(True)
            if self.macro_combo.findText(fname) == -1:
                self.macro_combo.addItem(fname)
                self.macro_combo.setCurrentText(fname)
            self.macro_combo.blockSignals(False)'''
            QMessageBox.information(self, "Saved", "Macro saved")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving macro:\n{e}")

    def apply_macro(self):

        code = self.editor.get_code().strip()
        self.executor.exec(code)

    def on_check(self, *args):
        _, index, roles = args
        if Qt.ItemDataRole.CheckStateRole.value in roles:
            if index.isValid() and index.internalPointer() is not None and index.internalPointer()._checked:
                self.executor.selection.append(index)
        ...