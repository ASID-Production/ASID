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

from abc import ABC, abstractmethod
from PySide6.QtWidgets import QDialog
from .ui.ui_instruments import Ui_Dialog
from PySide6 import QtCore
import time
from . import MAIN_WIDGET, MOLECULE_SYSTEMS
from PySide6.QtOpenGLWidgets import QOpenGLWidget
import numpy as np


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


class Dialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.opengl_widget = MAIN_WIDGET.findChild(QOpenGLWidget, "OpenGLWidget")
        self.old_filter = self.opengl_widget.eventFilterf
        self.mode = self.translate
        self.ui.pushButton_2.clicked.connect(self.deleteSel)

    def show(self):
        self.old_filter = self.opengl_widget.eventFilterf
        self.opengl_widget.eventFilterf = self.customEventFilter
        QDialog.show(self)

    def translate(self, dir):
        sel = self.opengl_widget.selection_model.selection()
        sel = [x.indexes()[0].internalPointer() for x in sel]
        x = dir.x()
        y = dir.y()
        x = x * 2 / (self.opengl_widget.width())
        y = -y * 2 / (self.opengl_widget.height())
        aspect_ratio = np.linalg.inv(self.opengl_widget.uniforms.aspect_ratio)
        scale = np.linalg.inv(self.opengl_widget.uniforms.scale)
        rotation = np.linalg.inv(self.opengl_widget.uniforms.rotation)
        pers = np.linalg.inv(self.opengl_widget.uniforms.perspective)
        coords = np.array([x,y,0.0,1.0])[np.newaxis].T
        coords = rotation @ scale @ aspect_ratio @ coords
        coords = np.squeeze(coords.T[0,:-1]).astype(np.float32)
        for p in sel:
            if p._atom:
                p.coord += coords
                p._atom.coord += coords
            elif p.coord is not None:
                p.coord += coords
        self.opengl_widget.update()
        ...

    def select_mol(self, pos):
        from .MoleculeClass import Bond
        from PySide6.QtCore import QItemSelectionModel
        def rec_mol(a, mem=None, mem_b=None):
            if not mem_b:
                mem_b = []
            if not mem:
                mem = [a]
            else:
                if a not in mem:
                    mem.append(a)
                else:
                    return mem, mem_b
            for bond in a.bonds():
                if bond not in mem_b:
                    mem_b.append(bond)
                rec_mol(bond.get(a), mem, mem_b)
            return mem, mem_b

        sel = self.opengl_widget.select(pos)
        mol = None

        for p in sel:
            if p._atom:
                if isinstance(p._atom, Bond):
                    mol = rec_mol(p._atom.parents()[0])
                else:
                    mol = rec_mol(p._atom)

        if mol:
            for a in mol[0]:
                index = self.opengl_widget.selection_model.model().index(0, 0, by_point=a.point())
                if a.point().pick is None:
                    a.point().addProperty('pick', 0.0)
                self.opengl_widget.selection_model.select(index, QItemSelectionModel.Select)
                self.opengl_widget.update()
            for b in mol[1]:
                p = b.point()[0].parent
                index = self.opengl_widget.selection_model.model().index(0, 0, by_point=p)
                if p.pick is None:
                    p.addProperty('pick', 0.0)
                self.opengl_widget.selection_model.select(index, QItemSelectionModel.Select)
                self.opengl_widget.update()

    def deleteSel(self):
        sel = self.opengl_widget.selection_model.selection()
        sel = [x.indexes()[0].internalPointer() for x in sel]
        for p in sel:
            if p._atom:
                p._atom.remove()
        self.opengl_widget.update()

    def customEventFilter(sl, self, obj: 'QObject', event: 'QEvent'):
        if event.type() == QtCore.QEvent.MouseButtonPress:
            self.pressed = True
            self.button = event.buttons()
            self.timer_pressed = time.perf_counter()
            self.pos = event.localPos()
        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.LeftButton:
            dir = event.localPos() - self.pos
            self.pos = event.localPos()
            if event.modifiers() == QtCore.Qt.NoModifier:
                self.rotate(dir)
                #sl.translate(dir)

        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.RightButton:
            dir = event.localPos() - self.pos
            self.pos = event.localPos()
            if event.modifiers() == QtCore.Qt.ControlModifier:
                pass
            if event.modifiers() == QtCore.Qt.NoModifier:
                self.roll(dir)

        if event.type() == QtCore.QEvent.MouseButtonRelease:
            self.pressed = False
            self.timer_pressed -= time.perf_counter()
            if self.timer_pressed >= -0.25:
                if self.button == QtCore.Qt.LeftButton and self.selection_model is not None:
                    sl.select_mol(event.localPos())
                elif self.button == QtCore.Qt.RightButton and self.selection_model is not None:
                    for index in self.selection_model.selectedIndexes():
                        if index.internalPointer().pick is None:
                            index.internalPointer().addProperty('pick', 0.0)
                        self.selection_model.model().setData(index, ('pick', 0.0), role=99)
                    self.selection_model.clearSelection()
                    self.update()
            else:
                pass
            self.button = None

        if event.type() == QtCore.QEvent.Wheel:
            if event.angleDelta().y() > 0:
                self.scale_func(1)
            else:
                self.scale_func(-1)

    def closeEvent(self, arg__1):
        self.opengl_widget.eventFilterf = self.old_filter
        QDialog.closeEvent(self, arg__1)


def execute():
    global DIALOG
    DIALOG = Dialog()
    DIALOG.show()