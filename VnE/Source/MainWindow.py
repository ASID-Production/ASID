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
import sys
from PySide6 import QtCore
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtOpenGL import QOpenGLDebugLogger, QOpenGLDebugMessage
import sys
from OpenGL.GL import *
import numpy as np
from . import Scenes, UniformBuffers
from .Facade import RenderFacade
import time
import logging

from .point_class import PointsList

from PySide6 import QtWidgets, QtGui
from PySide6.QtCore import *

from . import QtModels
from .QtModels import ListView, UniformListModel, TreeView, QtPointsTreeModel, SelectionModel, QtPointsPropertyModel
from .selector import MacroDockWidget
from PIL import Image


class OpenGlWidget(QOpenGLWidget):

    def __init__(self, parent, facade=None, scene=None, pipeline=None, model=None, **kwargs):
        super().__init__(parent)
        self.surface_format = QtGui.QSurfaceFormat()
        self.surface_format.setSamples(4)
        self.surface_format.setOption(QtGui.QSurfaceFormat.DebugContext)
        self.surface_format.setRenderableType(QtGui.QSurfaceFormat.OpenGL)
        self.surface_format.setProfile(QtGui.QSurfaceFormat.CoreProfile)
        self.surface_format.setMajorVersion(4)
        self.surface_format.setMinorVersion(6)
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
        self.debug_logger = QOpenGLDebugLogger()
        self.debug_logger.messageLogged.connect(self.log)
        self.button = None
        self.timer_pressed = 0
        self.pressed = False
        self.pos = [0, 0]
        self.select_fbo = None
        self.select_crbo, self.select_dsrbo = None, None
        self.screen_fbo = None
        self.screen_crbo, self.screen_dsrbo = None, None

    def log(self, msg):
        logging.debug(f'{msg.severity()} {msg.type()} {msg.id()} {msg.source()}\n{msg.message()}')

    def setSelectionModel(self, selection_model: QItemSelectionModel):
        self.selection_model = selection_model
        self.selection_model.selectionChanged.connect(self.update)

    def paintGL(self) -> None:
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.facade.drawScene(self.scene)

    def initializeGL(self) -> None:
        super().initializeGL()

        self.select_fbo = glGenFramebuffers(1)
        self.select_crbo, self.select_dsrbo = glGenRenderbuffers(2)
        def_rbo = int(glGetIntegerv(GL_RENDERBUFFER_BINDING))
        glBindFramebuffer(GL_FRAMEBUFFER, self.select_fbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.select_crbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, self.select_crbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.select_dsrbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, self.select_dsrbo)

        self.screen_fbo = glGenFramebuffers(1)
        self.screen_crbo, self.screen_dsrbo = glGenRenderbuffers(2)
        glBindFramebuffer(GL_FRAMEBUFFER, self.screen_fbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.screen_crbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, self.screen_crbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.screen_dsrbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, self.screen_dsrbo)
        glBindRenderbuffer(GL_RENDERBUFFER, def_rbo)
        glBindFramebuffer(GL_FRAMEBUFFER, self.context().defaultFramebufferObject())

        logging.debug(f'OpenGL context profile: {self.surface_format.profile()} {self.surface_format.renderableType()} {self.surface_format.majorVersion()}.{self.surface_format.minorVersion()}')
        if self.context().hasExtension(QByteArray("GL_KHR_debug".encode())):
            logging.debug('GL_KHR_debug supported')
            self.debug_logger.initialize()
            self.debug_logger.enableMessages(QOpenGLDebugMessage.AnySource, QOpenGLDebugMessage.AnyType, QOpenGLDebugMessage.AnySeverity)
            self.debug_logger.startLogging(QOpenGLDebugLogger.SynchronousLogging)
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
        #glClearColor(1.0,0.0,0.0,1.0)

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
        self.makeCurrent()
        size = glGetIntegerv(GL_VIEWPORT)
        scale = [size[2]/self.width(), size[3]/self.height()]
        pos_new = [int(pos.x()*scale[0]), int((self.height()-int(pos.y()))*scale[1])]
        glBindFramebuffer(GL_FRAMEBUFFER, self.select_fbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.select_crbo)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32UI, *size[2:])
        glBindRenderbuffer(GL_RENDERBUFFER, self.select_dsrbo)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_STENCIL, *size[2:])

        glBindFramebuffer(GL_FRAMEBUFFER, self.select_fbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, self.select_dsrbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, self.select_crbo)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glDrawBuffers(1, GL_COLOR_ATTACHMENT0)
        self.facade.drawScene(self.scene, mode='SELECT')
        glFlush()
        glBindFramebuffer(GL_FRAMEBUFFER, self.select_fbo)
        glReadBuffer(GL_COLOR_ATTACHMENT0)

        c = np.zeros((4,), dtype=np.uint32)
        glReadPixels(*pos_new, 1, 1, GL_RGBA_INTEGER, GL_UNSIGNED_INT, c)

        pipeline_id = int.from_bytes(c[2:].tobytes(), 'little')
        point_pos = int(c[0])
        point_count = int(c[1])
        points = []
        for obs in QtModels.SINGLE_OBSERVER.obs_dict.values():
            if pipeline_id == obs._pipeline:
                points = obs._points[point_pos:point_pos + point_count]
                break
        self.makeCurrent()

        if points:
            for point in points:
                index = self.selection_model.model().index(0, 0, by_point=point)
                if point.pick is None:
                    point.addProperty('pick', 0.0)
                if point.pick == 1.0:
                    self.selection_model.select(index, QItemSelectionModel.Deselect)
                else:
                    self.selection_model.select(index, QItemSelectionModel.Select)
            self.update()
        return points

    def screenshot(self, filename, res=None):
        if not res:
            res = glGetIntegerv(GL_VIEWPORT)[2:]
        self.makeCurrent()
        old_view = glGetIntegerv(GL_VIEWPORT)
        glViewport(0,0,*res)
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'wh', res)
        glBindFramebuffer(GL_FRAMEBUFFER, self.screen_fbo)
        glBindRenderbuffer(GL_RENDERBUFFER, self.screen_crbo)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F, *res)
        glBindRenderbuffer(GL_RENDERBUFFER, self.screen_dsrbo)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_STENCIL, *res)

        glBindFramebuffer(GL_FRAMEBUFFER, self.screen_fbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, self.screen_dsrbo)
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, self.screen_crbo)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glDrawBuffers(1, GL_COLOR_ATTACHMENT0)
        self.facade.drawScene(self.scene)
        glFlush()
        glBindFramebuffer(GL_FRAMEBUFFER, self.screen_fbo)
        glReadBuffer(GL_COLOR_ATTACHMENT0)

        c = np.zeros((res[1], res[0], 4), dtype=np.uint8)
        glReadPixels(0, 0, *res, GL_RGBA, GL_UNSIGNED_BYTE, c)
        c = c[::-1]
        im = Image.fromarray(c, 'RGBA')
        im.save(filename)
        glViewport(*old_view)
        self.facade.changeUniformBufferProperty(self.uniforms_id, 'wh', old_view[2:])

    def eventFilter(self, obj: 'QObject', event: 'QEvent') -> bool:
        self.eventFilterf(self, obj, event)

        return super().eventFilter(obj, event)

    @staticmethod
    def eventFilterf(self, obj, event):
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

    def getUniforms(self):
        return self.facade.getInst(self.uniforms_id)


class UniformWid(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setWindowIcon(QtGui.QIcon('Source/ico.svg'))
        self.setWindowTitle('Uniforms')
        self.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint, True)
        self.listView = ListView(parent=self)
        self.hlayout = QtWidgets.QHBoxLayout()
        self.hlayout.addWidget(self.listView)
        self.setLayout(self.hlayout)
        self.hlayout.setSpacing(0)
        self.hlayout.setContentsMargins(0,0,0,0)
        self.model = UniformListModel(self)

    def setUniforms(self, uniform):
        self.model.setModelData(uniform)

    def setModel(self, model):
        self.model = model
        self.listView.setModel(model)


class SaveScreenDialog(QtWidgets.QDialog):

    def __init__(self, function, parent=None):
        super().__init__(parent)
        self.function = function
        self.setWindowTitle("Save screen")
        self.setMinimumWidth(400)

        self.file_path = ""
        self.width = 800
        self.height = 600

        self._create_file_selection()
        self._create_resolution_inputs()
        self._create_buttons()

        self._setup_layout()

    def _create_file_selection(self):
        self.file_path_edit = QtWidgets.QLineEdit()
        self.file_path_edit.setPlaceholderText("screen.png")

        self.browse_button = QtWidgets.QPushButton("...")
        self.browse_button.clicked.connect(self._browse_file)

    def _create_resolution_inputs(self):
        self.width_spin = QtWidgets.QSpinBox()
        self.width_spin.setRange(1, 10000)
        self.width_spin.setValue(self.width)

        self.height_spin = QtWidgets.QSpinBox()
        self.height_spin.setRange(1, 10000)
        self.height_spin.setValue(self.height)

    def _create_buttons(self):
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self._on_accept)
        self.button_box.rejected.connect(self.reject)

    def _setup_layout(self):
        main_layout = QtWidgets.QVBoxLayout(self)

        file_layout = QtWidgets.QHBoxLayout()
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.browse_button)
        main_layout.addLayout(file_layout)

        res_layout = QtWidgets.QHBoxLayout()
        res_layout.addWidget(QtWidgets.QLabel("Width:"))
        res_layout.addWidget(self.width_spin)
        res_layout.addWidget(QtWidgets.QLabel("Height:"))
        res_layout.addWidget(self.height_spin)
        res_layout.addStretch()
        main_layout.addLayout(res_layout)

        main_layout.addWidget(self.button_box)

    def _browse_file(self):
        """Открывает диалог сохранения файла и обновляет поле пути."""
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save image",
            "",
            "Image (*.png *.tiff)"
        )
        if file_path:
            self.file_path_edit.setText(file_path)

    def _on_accept(self):
        if not self.file_path_edit.text():
            return

        self.file_path = self.file_path_edit.text()
        self.width = self.width_spin.value()
        self.height = self.height_spin.value()
        if self.file_path:
            self.function(self.file_path, [self.width, self.height])
        self.accept()


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):

        self.points_list = PointsList()
        super().__init__()
        self.setWindowIcon(QtGui.QIcon('Source/ico.ico'))
        self.setWindowTitle('ASID View & Explore')
        widget = QtWidgets.QWidget()
        #self.frame = QtWidgets.QFrame(parent=widget)
        #self.frame.setObjectName('lists_frame')
        self.opengl_frame = QtWidgets.QFrame(parent=widget)
        self.opengl_frame.setObjectName('opengl_frame')
        self.opengl_frame.setLayout(QtWidgets.QVBoxLayout())
        self.opengl_frame.layout().setSpacing(0)
        #self.setCentralWidget(self.opengl_frame)
        self.setCentralWidget(widget)
        #vlayout = QtWidgets.QVBoxLayout()
        hlayout = QtWidgets.QHBoxLayout()
        hlayout.setSpacing(0)

        self.treeView = TreeView()
        self.listView = ListView()
        self.treeDock = QtWidgets.QDockWidget('Points list', self)
        self.treeDock.setObjectName('treeDock')
        self.treeDock.setWidget(self.treeView)
        self.treeDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetFloatable)
        self.listDock = QtWidgets.QDockWidget('Props list', self)
        self.listDock.setObjectName('listDock')
        self.listDock.setWidget(self.listView)
        self.listDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetFloatable)

        self.model = QtPointsTreeModel(parent=self.treeView, data=self.points_list, main_window=self)
        self.selection_model = SelectionModel(model=self.model)
        self.list_model = QtPointsPropertyModel(data=self.points_list, main_window=self)

        self.selection_model.newSelection.connect(lambda args: self.list_model.setSelected(*args))

        self.treeView.setModel(self.model)
        self.treeView.setSelectionModel(self.selection_model)
        self.listView.setModel(self.list_model)

        self.uniformWid = UniformWid()
        self.uniformDock = QtWidgets.QDockWidget('Uniforms', self)
        self.uniformDock.setObjectName('uniformDock')
        self.uniformDock.setWidget(self.uniformWid)
        self.uniformDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetFloatable | QtWidgets.QDockWidget.DockWidgetClosable)
        self.uniformDock.hide()

        macro_dir = os.path.normpath(f'{os.path.dirname(__file__)}/macros')
        os.makedirs(macro_dir, exist_ok=True)
        self.macroDock = MacroDockWidget(macro_dir, self, self.model, self.selection_model)
        self.macroDock.setObjectName('macroDock')
        self.macroDock.hide()

        self.addDockWidget(Qt.LeftDockWidgetArea, self.uniformDock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.treeDock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.listDock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.macroDock)

        self.opengl_widget = OpenGlWidget(parent=widget, model=self.points_list, uniformWidget=self.uniformWid)
        self.opengl_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        self.opengl_widget.setSelectionModel(self.selection_model)
        self.opengl_frame.layout().addWidget(self.opengl_widget)

        #self.listView.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        #self.listView.setMinimumWidth(300)
        #self.treeView.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        #self.treeView.setMinimumWidth(300)
        #self.frame.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Expanding)
        #self.frame.setMinimumWidth(300)

        #vlayout.addWidget(self.treeView)
        #vlayout.addWidget(self.listView)
        #self.frame.setLayout(vlayout)
        #hlayout.addWidget(self.frame)
        hlayout.addWidget(self.opengl_frame)

        self.uniformModel = UniformListModel()
        self.uniformWid.setModel(self.uniformModel)

        widget.setLayout(hlayout)
        self.resize(1280, 720)

        from . import Extensions

        self.menu = self.menuBar()
        self.opengl_widget.setObjectName('OpenGLWidget')
        self.menu.setObjectName('MenuBar')
        self.extension_menu = Extensions.getMenu(self.model, self.uniformModel, main_widget=widget, main_menu=self.menu)
        self.uniformAction = self.menu.addAction('Uniforms')
        self.macroAction = self.menu.addAction('Macros')
        self.screenshotAction = self.menu.addAction('Screenshot')
        self.saveStateAction = self.menu.addAction('Save state')
        self.loadStateAction = self.menu.addAction('Load state')
        self.menu.addMenu(self.extension_menu)
        self.uniformAction.triggered.connect(lambda: self.toggleWidget(self.uniformDock))
        self.macroAction.triggered.connect(lambda: self.toggleWidget(self.macroDock))
        self.screenshotAction.triggered.connect(self.screenshot)
        self.saveStateAction.triggered.connect(self.saveSceneState)
        self.loadStateAction.triggered.connect(self.loadState)

        self.about = self.menu.addAction('About')

        self.about_dialog = QtWidgets.QDialog(self)
        self.about_dialog.setWindowTitle('About')
        self.about_dialog.setLayout(QtWidgets.QVBoxLayout())
        self.about_dialog.layout().addWidget(QtWidgets.QLabel("ASID 1.0.1\nAuthors:\nAlexander A. Korlyukov (head),\nAlexander D. Volodin (author of cpplib),\nPetr A. Buikin (author of api_database),\nAlexander R. Romanenko (author of VnE)"))

        self.about.triggered.connect(self.about_dialog.show)

        settings = QtCore.QSettings('ASID', 'LRSI')
        state = settings.value("windowState")
        geom = settings.value("windowGeom")
        if state:
            self.restoreState(state, version=0)
        if geom:
            self.restoreGeometry(geom)

    def toggleWidget(self, widget):
        widget.setVisible(not widget.isVisible())

    def closeEvent(self, event, *args, **kwargs):
        settings = QtCore.QSettings('ASID', 'LRSI')
        settings.setValue("windowState", self.saveState(version=0))
        settings.setValue("windowGeom", self.saveGeometry())
        settings.setValue("macroDock/splitter", self.macroDock.splitter.saveState())
        QtWidgets.QApplication.quit()

    def screenshot(self):
        self.screen_dialog = SaveScreenDialog(self.opengl_widget.screenshot)
        self.screen_dialog.show()

    def saveSceneState(self):
        import json
        d = []
        ps = []
        def recDict(p):
            if p not in ps:
                d.append(p.toDict())
                ps.append(p)
            else:
                return
            for c in p.children:
                recDict(c)
        recDict(self.points_list)
        uniforms = self.opengl_widget.uniforms
        scene_state = {k: getattr(uniforms, k) for k in uniforms.getInfo()}
        def rec(l):
            if isinstance(l, np.ndarray):
                l = list(l)
            for i, v in enumerate(l):
                if isinstance(v, np.ndarray):
                    l[i] = rec(v)
                elif isinstance(v, np.int32) or isinstance(v, np.intc):
                    l[i] = int(v)
                elif isinstance(v, np.float32):
                    l[i] = float(v)
            return l

        for k in scene_state:
            if isinstance(scene_state[k], np.ndarray):
                scene_state[k] = scene_state[k].tolist()
            elif isinstance(scene_state[k], np.int32) or isinstance(scene_state[k], np.intc):
                scene_state[k] = int(scene_state[k])
            elif isinstance(scene_state[k], np.float32):
                scene_state[k] = float(scene_state[k])
        d = {'scene_state': scene_state,
             'points': d}
        j = json.dumps(d, indent=2)
        f, _ = QtWidgets.QFileDialog.getSaveFileName(caption='Scene state', filter='*.json')
        if f:
            f = open(f, 'w')
            f.write(j)
            f.close()
        ...

    def loadState(self):
        import json
        from . import point_class
        file, _ = QtWidgets.QFileDialog.getOpenFileName(caption='Scene state', filter='*.json')
        if not file:
            return
        f = open(file, 'r').read()
        d = json.loads(f)
        t = {'Point': point_class.Point, 'PointsList': point_class.PointsList}
        points = {}
        obs = {}
        for p in d['points']:
            obj = t[p['type']]()
            points[p['seq']] = obj
            for prop, v in p['static_props'].items():
                obj.addProperty(prop, v)
            for ob in p['obs']:
                l = obs.get(ob, None)
                if not l:
                    obs[ob] = [obj]
                else:
                    l.append(obj)
        for p in d['points']:
            obj = points[p['seq']]
            parent = None if p['parent'] is None else points[p['parent']]
            if parent is None:
                root = obj
            else:
                parent.addChild(obj)

            children = [points[x] for x in p['children']]
            if children:
                for ch in children:
                    obj.addChild(ch)

            for prop, seq in p['dynamic_props'].items():
                obj.addProperty(prop, points[seq])

        root.addProperty('name', os.path.basename(file))
        self.points_list.addChild(root)
        self.model.update()

        for ob, ps in obs.items():
            for p in ps:
                ind = self.model.index(by_point=p)
                self.model.attachObserver(ind, ob)




        for k, v in d['scene_state'].items():
            if isinstance(v, list):
                v = np.array(v, dtype=np.float32)
            self.opengl_widget.uniforms.__setattr__(k, v)

    def updateOpenGL(self):
        self.opengl_widget.update()


def show():
    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()