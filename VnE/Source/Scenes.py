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


from OpenGL.GL import *
from abc import ABC, abstractmethod
import numpy as np

from .ShaderPipelines import *
from .UniformBuffers import SceneUniformBuffer

import debug

class aScene(ABC):
    def __init__(self, context):
        self.shader_pipelines = []
        self.uniform_buffers = []
        self._context = context

    def makeCurrent(self):
        self._context.makeCurrent()

    def addShaderPipeline(self, pipeline_cls=None, pipeline_inst=None):
        self._context.makeCurrent()
        if pipeline_inst:
            pipe_line = pipeline_inst
        elif pipeline_cls:
            pipe_line = pipeline_cls()
        if pipe_line not in self.shader_pipelines:
            self.shader_pipelines.append(pipe_line)
        return pipe_line

    def addUniformBuffer(self, uniform_buffer_cls=None, uniform_buffer_inst=None):
        self._context.makeCurrent()
        if uniform_buffer_inst:
            buffer = uniform_buffer_inst
        elif uniform_buffer_cls:
            buffer = uniform_buffer_cls()
        if buffer not in self.uniform_buffers:
            self.uniform_buffers.append(buffer)
        return buffer

    def getInfo(self):
        return

    @abstractmethod
    def draw(self):
        self._context.makeCurrent()
        for buffer in self.uniform_buffers:
            buffer.bind()
        for pipeline in self.shader_pipelines:
            pipeline.draw()


class Scene(aScene):

    def __init__(self, context):
        super().__init__(context=context)
        self.conf = False

    def draw(self):
        self.makeCurrent()
        if self.conf:
            for buffer in self.uniform_buffers:
                buffer.bind()
            for pipe_line in self.shader_pipelines:
                pipe_line.draw()
        else:
            glEnable(GL_DEPTH_TEST)
            glEnable(GL_MULTISAMPLE)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glEnable(GL_BLEND)
            glClearColor(1.0, 1.0, 1.0, 1.0)
            glUseProgram(0)
            self.conf = True
            self.draw()


class TestScene(aScene):

    def __init__(self, context):
        self._context = context
        self.shaderProgram = PlaneShaderPipeline()
        self.uniform_buffers = [SceneUniformBuffer()]

    def draw(self):
        self._context.makeCurrent()
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_MULTISAMPLE)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_BLEND)
        glClearColor(1.0, 1.0, 1.0, 1.0)
        glUseProgram(0)
        for buffer in self.uniform_buffers:
            buffer.bind()
        self.shaderProgram.draw()


if __name__ == '__main__':
    from PySide6 import QtCore, QtWidgets
    import sys
    class OpenGlWidget(QtWidgets.QOpenGLWidget):

        def __init__(self, parent):
            super().__init__(parent)
            self.timer = QtCore.QTimer()
            self.timer.timeout.connect(self.update)

        def paintGL(self) -> None:
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            self.scene.draw()
            self.delta_time += 1
            self.scene.uniform_buffers[0].y_rotation = np.pi / 1000
            self.scene.uniform_buffers[0].x_rotation = np.pi / 1000
            self.scene.uniform_buffers[0].z_rotation = np.pi / 1000

            self.timer.start(1)

        def initializeGL(self) -> None:
            super().initializeGL()
            self.delta_time = 0
            self.it = 100
            self.dec = True
            self.scene = TestScene(self)
            self.scene.uniform_buffers[0].perspective = np.array([-10, 10, 10, -10, 25, 100])
            self.scene.uniform_buffers[0].scale = np.array([1, 1, 1])
            self.scene.shaderProgram.add_shader_data()
            sl = [-5, 5, 1]
            self.scene.shaderProgram.shader_data[0].addData(0, np.array([[x, y, z] for x in range(*sl) for y in range(*sl) for z in range(-50, -40, 1)], dtype=np.float32))
            self.scene.shaderProgram.shader_data[0].addData(1, np.array([[1, 0, 0, 1] for x in range(*sl) for y in range(*sl) for z in range(-50, -40, 1)], dtype=np.float32))

        def resizeGL(self, w: int, h: int) -> None:
            self.scene.uniform_buffers[0].aspect_ratio = w / h


    class MainWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super(QtWidgets.QMainWindow, self).__init__()
            layout = QtWidgets.QVBoxLayout()
            self.opengl = OpenGlWidget(parent=self)
            layout.addWidget(self.opengl)
            widget = QtWidgets.QWidget()
            widget.setLayout(layout)
            self.setCentralWidget(widget)


    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()
