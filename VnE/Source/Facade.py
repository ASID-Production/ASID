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


from . import ShaderPipelines
from . import ShaderDataObjects
import debug


class RenderFacade:

    def __init__(self, context):
        self.ids = {}
        self.scenes = {}
        self.pipe_lines = {}
        self.uniform_buffers = {}
        self.data_buffers = {}
        self.context = context

    def getInfo(self, id):
        return self.ids[id].getInfo()

    def makeCurrent(self):
        if self.context is None:
            raise Exception
        self.context.makeCurrent()
        ShaderPipelines.setCreator(self.context)
        ShaderDataObjects.setCreator(self.context)

    def addScene(self, scene_cls):
        self.makeCurrent()
        scene = scene_cls(self.context)

        self.scenes[id(scene)] = scene
        self.ids[id(scene)] = scene
        return id(scene)

    def drawScene(self, scene_id):
        self.makeCurrent()
        scene = self.scenes[scene_id]
        scene.draw()

    def addPipelineToScene(self, scene_id, pipeline_cls=None, pipeline_inst=None):
        self.makeCurrent()
        scene = self.scenes[scene_id]

        if pipeline_inst:
            scene.addShaderPipeline(pipeline_inst=pipeline_inst)
            self.pipe_lines[id(pipeline_inst)] = pipeline_inst
            self.ids[id(pipeline_inst)] = pipeline_inst
            return id(pipeline_inst)
        elif pipeline_cls:
            pipeline_inst = pipeline_cls()
            scene.addShaderPipeline(pipeline_inst=pipeline_inst)
            self.pipe_lines[id(pipeline_inst)] = pipeline_inst
            self.ids[id(pipeline_inst)] = pipeline_inst
            return id(pipeline_inst)

    def changePipelineUniforms(self, pipeline_id, property, value):
        self.makeCurrent()
        pipeline = self.pipe_lines[pipeline_id]
        pipeline.changeUniform(property, value)

    def addUniformBufferToScene(self, scene_id, uniform_buffer_inst=None, uniform_buffer_cls=None):
        self.makeCurrent()
        scene = self.scenes[scene_id]

        if uniform_buffer_inst:
            scene.addUniformBuffer(uniform_buffer_inst=uniform_buffer_inst)
            self.uniform_buffers[id(uniform_buffer_inst)] = uniform_buffer_inst
            self.ids[id(uniform_buffer_inst)] = uniform_buffer_inst
            return id(uniform_buffer_inst)
        elif uniform_buffer_cls:
            uniform_buffer_inst = uniform_buffer_cls(self.makeCurrent)
            scene.addUniformBuffer(uniform_buffer_inst=uniform_buffer_inst)
            self.uniform_buffers[id(uniform_buffer_inst)] = uniform_buffer_inst
            self.ids[id(uniform_buffer_inst)] = uniform_buffer_inst
            return id(uniform_buffer_inst)

    def changeUniformBufferProperty(self, uniform_buffer_id, property, data):
        self.makeCurrent()
        buffer = self.uniform_buffers[uniform_buffer_id]
        buffer.__setattr__(property, data)

    def addDataBufferToPipeline(self, pipeline_id, shader_data_inst=None):
        self.makeCurrent()
        pipeline = self.pipe_lines[pipeline_id]
        if shader_data_inst:
            pipeline.add_shader_data(shader_data_inst=shader_data_inst)
        else:
            shader_data_inst = pipeline.add_shader_data()
        self.data_buffers[id(shader_data_inst)] = shader_data_inst
        self.ids[id(shader_data_inst)] = shader_data_inst
        return id(shader_data_inst)

    def addDataToShaderData(self, shader_data_id, attr, data):
        self.makeCurrent()
        shader_data = self.data_buffers[shader_data_id]
        shader_data.addData(attr, data)

    def deleteDataInShaderData(self, shader_data_id, attr, offset, bytes):
        self.makeCurrent()
        shader_data = self.data_buffers[shader_data_id]
        shader_data.deleteData(attr, offset, bytes)

    def replaceDataInShaderData(self, shader_data_id, attr, data, offset):
        self.makeCurrent()
        shader_data = self.data_buffers[shader_data_id]
        shader_data.replaceData(attr, data, offset)

    def getInst(self, id):
        return self.ids[id]

