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


import sys
from PySide6 import QtCore
from PySide6.QtOpenGLWidgets import QOpenGLWidget
import sys
from OpenGL.GL import *
import numpy as np
from . import Scenes, UniformBuffers
from .Facade import RenderFacade
import time

from .point_class import PointsList

from PySide6 import QtWidgets, QtGui
from PySide6.QtCore import *

from . import QtModels
from .QtModels import ListView, UniformListModel, TreeView, QtPointsTreeModel, SelectionModel, QtPointsPropertyModel


class OpenGlWidget(QOpenGLWidget):

    def __init__(self, parent, facade=None, scene=None, pipeline=None, model=None, **kwargs):
        super().__init__(parent)
        self.surface_format = QtGui.QSurfaceFormat()
        self.surface_format.setSamples(4)
        self.setFormat(self.surface_format)
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update)
        self.facade = facade
        self.scene = scene
        self.pipeline = pipeline
        self._model = model
        self.uniformWidget = kwargs.get('uniformWidget', None)
        self.installEventFilter(self)
        self.selection_model = None
        self.model = None

        self.button = None
        self.timer_pressed = 0
        self.pressed = False
        self.pos = [0, 0]

    def setSelectionModel(self, selection_model: QItemSelectionModel):
        self.selection_model = selection_model
        self.selection_model.selectionChanged.connect(self.update)

    def paintGL(self) -> None:
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.facade.drawScene(self.scene)

    def initializeGL(self) -> None:
        super().initializeGL()
        self.delta_time = 0
        self.it = 100
        if self.facade is None:
            self.facade = RenderFacade(self)
            self.scene = self.facade.addScene(scene_cls=Scenes.Scene)
        from . import Observers

        SINGLE_OBSERVER = Observers.SingleObserver(self.facade, self.scene)
        QtModels.SINGLE_OBSERVER = SINGLE_OBSERVER
        self.label_observer = SINGLE_OBSERVER.getObserver(Observers.LabelObserver)
        self.label_observer.set_wh([self.width(), self.height()])
        pipeline = self.label_observer.getPipeline()
        self.facade.changePipelineUniforms(pipeline, 'const_scale', ctypes.c_float(150.0))
        self.uniforms_id = self.facade.addUniformBufferToScene(self.scene, uniform_buffer_cls=UniformBuffers.SceneUniformBuffer)
        self.uniforms = self.getUniforms()
        if self.uniformWidget is not None:
            self.uniformWidget.setUniforms(self.facade.getInst(self.uniforms_id))
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'perspective', np.array([-100, 100, 100, -100, 100, 2500]))
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'scene_shift', -500)
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'scale', np.array([5, 5, 5]))

    def resizeGL(self, w: int, h: int) -> None:
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'wh', [w, h])

    def rotate(self, dir):
        """
        Rotation scene along x,y axes function
        :param dir: x,y direction of rotation in pixels
        :type dir: QtCore.QPointF
        :return: None
        """
        x = -dir.x()
        y = dir.y()

        self.makeCurrent()
        self.uniforms.y_rotation = (x / self.width()) * 4 * np.pi
        self.uniforms.x_rotation = (y / self.height()) * 4 * np.pi
        self.update()

    def scale_func(self, zoom):
        """
        Scale scene function
        :param zoom: 1 or -1, zoom or unzoom scene
        :type zoom: int
        :return: None
        """
        self.makeCurrent()
        self.uniforms.x_scale *= 1.1 ** zoom
        self.uniforms.y_scale *= 1.1 ** zoom
        self.uniforms.z_scale *= 1.1 ** zoom
        self.update()

    def roll(self, dir):
        """
        Rotation of scene along z axis function
        :param dir: rotation direction and velocity, only along x direction of mouse movement
        :type dir: QtCore.QPointF
        :return: None
        """
        z = -dir.x()
        self.makeCurrent()
        self.uniforms.z_rotation = (z / self.width()) * 4 * np.pi
        self.update()

    def translate(self, dir):
        """
        Translate of the scene function
        :param dir: x,y direction of translation in pixels
        :type dir: QtCore.QPointF
        :return: None
        """
        x = dir.x()
        y = dir.y()
        self.makeCurrent()
        self.uniforms.x_translation += x * 2 / (self.width())
        self.uniforms.y_translation -= y * 2 / (self.height())
        self.update()

    def select(self, pos):
        if self.selection_model is None:
            return
        mod = self.uniforms.translation @ self.uniforms.perspective @ self.uniforms.aspect_ratio @ self.uniforms._rotation_point_matr @ self.uniforms.scale @ self.uniforms.rotation @ self.uniforms._r_rotation_point_matr @ self.uniforms._scene_shift
        tol = (self.uniforms.aspect_ratio @ self.uniforms._rotation_point_matr @ self.uniforms.scale @ self.uniforms._r_rotation_point_matr)[:2,:2] @ np.array([[1], [1]])
        tol = tol.flatten()
        pos = [(pos.x()/self.width()) * 2 - 1, (pos.y()/self.height()) * (-2) + 1]
        root = self.selection_model.model().getRoot()
        selected = root.select(pos, tol=tol, mod=mod)
        if selected is not None:
            selected = selected[0]
            index = self.selection_model.model().index(0, 0, by_point=selected)
            if selected.pick is None:
                selected.addProperty('pick', 0.0)
            if selected.pick == 1.0:
                self.selection_model.select(index, QItemSelectionModel.Deselect)
            else:
                self.selection_model.select(index, QItemSelectionModel.Select)
        self.update()

    def eventFilter(self, obj: 'QObject', event: 'QEvent') -> bool:
        if event.type() == QtCore.QEvent.MouseButtonPress:
            self.pressed = True
            self.button = event.buttons()
            self.timer_pressed = time.perf_counter()
            self.pos = event.localPos()
        if event.type() == QtCore.QEvent.MouseMove and event.buttons() == QtCore.Qt.LeftButton:
            dir = event.localPos() - self.pos
            self.pos = event.localPos()
            if event.modifiers() == QtCore.Qt.ControlModifier:
                self.translate(dir)
            if event.modifiers() == QtCore.Qt.NoModifier:
                self.rotate(dir)

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
                    self.select(event.localPos())
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

        return super().eventFilter(obj, event)

    def getUniforms(self):
        return self.facade.getInst(self.uniforms_id)


class UniformWid(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle('Uniforms')
        self.listView = ListView(parent=self)
        self.hlayout = QtWidgets.QHBoxLayout()
        self.hlayout.addWidget(self.listView)
        self.setLayout(self.hlayout)
        self.model = UniformListModel()

    def setUniforms(self, uniform):
        self.model.setModelData(uniform)

    def setModel(self, model):
        self.model = model
        self.listView.setModel(model)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):

        self.points_list = PointsList()
        super().__init__()
        self.setWindowTitle('ASID View & Explore')
        widget = QtWidgets.QWidget()
        frame = QtWidgets.QFrame(parent=widget)
        opengl_frame = QtWidgets.QFrame(parent=widget)
        vlayout = QtWidgets.QVBoxLayout()
        hlayout = QtWidgets.QHBoxLayout()

        self.treeView = TreeView()
        self.listView = ListView()

        self.model = QtPointsTreeModel(parent=self.treeView, data=self.points_list)
        self.selection_model = SelectionModel(model=self.model)
        self.list_model = QtPointsPropertyModel(data=self.points_list)

        self.selection_model.newSelection.connect(lambda args: self.list_model.setSelected(*args))

        self.treeView.setModel(self.model)
        self.treeView.setSelectionModel(self.selection_model)
        self.listView.setModel(self.list_model)

        self.uniformWid = UniformWid()

        self.opengl_widget = OpenGlWidget(parent=widget, model=self.points_list, uniformWidget=self.uniformWid)
        self.opengl_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        self.opengl_widget.setSelectionModel(self.selection_model)

        self.listView.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        self.listView.setMinimumWidth(300)
        self.treeView.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        self.treeView.setMinimumWidth(300)
        frame.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        frame.setMinimumWidth(300)

        vlayout.addWidget(self.treeView)
        vlayout.addWidget(self.listView)
        frame.setLayout(vlayout)
        hlayout.addWidget(frame)
        hlayout.addWidget(self.opengl_widget)

        self.uniformModel = UniformListModel()
        self.uniformWid.setModel(self.uniformModel)

        from . import Extensions

        self.menu = self.menuBar()
        self.extension_menu = Extensions.getMenu(self.model, self.uniformModel)
        self.uniformAction = self.menu.addAction('Uniforms')
        self.menu.addMenu(self.extension_menu)
        self.uniformAction.triggered.connect(self.uniformWid.show)

        widget.setLayout(hlayout)
        self.resize(1280, 720)

        self.setCentralWidget(widget)


def show():
    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()