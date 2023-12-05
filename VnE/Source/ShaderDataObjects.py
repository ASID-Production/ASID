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
import numpy as np
from typing import Tuple, Dict, List


class VAOCreator:
    types = {np.byte: GL_BYTE,
             np.ubyte: GL_UNSIGNED_BYTE,
             np.intc: GL_INT,
             np.uintc: GL_UNSIGNED_INT,
             np.float32: GL_FLOAT}

    def __init__(self):
        self.VAOs: {Tuple[int]: Dict[Tuple[str]: int]} = {}

    def create_VAO(self, VAOFormat):
        shape = tuple([x[0] for x in VAOFormat])
        data_types = tuple([x[1] for x in VAOFormat])
        if self.VAOs.get(shape, None):
            if self.VAOs[shape].get(data_types, None):
                return self.VAOs[shape].get(data_types)
        else:
            self.VAOs[shape] = {}

        VAO = glGenVertexArrays(1)
        glBindVertexArray(VAO)
        self.VAOs[shape][data_types] = VAO
        for i, attr in enumerate(shape):
            glVertexAttribFormat(i, attr, ShaderData.types[data_types[i]], GL_FALSE, 0)
            glVertexAttribBinding(i, i)
            glEnableVertexAttribArray(i)
        glBindVertexArray(0)
        return VAO


class ShaderData:
    types = {np.byte: GL_BYTE,
             np.ubyte: GL_UNSIGNED_BYTE,
             np.intc: GL_INT,
             np.uintc: GL_UNSIGNED_INT,
             np.float32: GL_FLOAT}

    def __init__(self, VAOFormat, allocation_size):

        self.shape = [x[0] for x in VAOFormat]
        self.shape_bytes = [x[0] * np.dtype(x[1]).itemsize for x in VAOFormat]
        self.data_types = [x[1] for x in VAOFormat]

        self.VAO = VAO_CREATOR.create_VAO(VAOFormat)

        self._VBOSize = {}
        self.VBOs = []
        for i,attr in enumerate(self.shape):
            id = glGenBuffers(1)
            self._VBOSize[id] = 0
            VBO = {'id': id,
                   'size': 0,
                   'size_bytes': 0 * np.dtype(self.data_types[i]).itemsize * attr,
                   'data': None,
                   'dtype': ShaderData.types[self.data_types[i]],
                   'dtype_bytes': np.dtype(self.data_types[i]).itemsize,
                   'shape': attr,
                   'stride': attr * np.dtype(self.data_types[i]).itemsize,
                   'allocated_size': allocation_size}
            self.VBOs.append(VBO)
            glBindBuffer(GL_ARRAY_BUFFER, VBO['id'])
            glBufferData(GL_ARRAY_BUFFER, VBO['allocated_size'], VBO['data'], GL_STATIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def getInfo(self):
        return self.VBOs

    @staticmethod
    def _validateData(VBO, data):
        if isinstance(data, np.ndarray):
            if ShaderData.types[data.dtype.type] is not VBO['dtype']:
                raise TypeError(f"Invalid data type, expected {VBO['dtype']}, got {ShaderData.types[data.dtype]}")

    @staticmethod
    def _validateVBORange(VBO, data, offset):
        if offset >= VBO['size_bytes']:
            raise ValueError(f'Offset out of range, expected offset < {VBO["size_bytes"]}, got {offset}')
        last_pos = offset + data.nbytes
        if last_pos > VBO['size_bytes']:
            raise ValueError(f"Replace location out of VBO range, expected from 0 to {VBO['size_bytes']} bytes, got {last_pos} bytes")
        return

    @staticmethod
    def _checkForMemoryAllocation(VBO, data, offset):
        last_pos = offset + data.nbytes
        if last_pos <= VBO['allocated_size']:
            return
        if last_pos - VBO['allocated_size'] >= VBO['allocated_size']:
            glNamedBufferData(VBO['id'], last_pos, VBO['data'], GL_STATIC_DRAW)
            VBO['allocated_size'] = last_pos
        else:
            glNamedBufferData(VBO['id'], VBO['allocated_size'] * 2, None, GL_STATIC_DRAW)
            glNamedBufferSubData(VBO['id'], 0, VBO['data'].nbytes, VBO['data'])
            VBO['allocated_size'] = VBO['allocated_size'] * 2

    @property
    def size(self):
        """

        :return: number of points to render
        """
        return min(self._VBOSize.values())

    @size.setter
    def size(self, value: List[int]):
        """

        :param value: [attr, size in points] attr - VBO index in self.VBOs
        :return:
        """
        if type(value) is not List and len(value) != 2:
            raise ValueError('Value is not list or not list of size 2')
        if value[0] in [VBO['id'] for VBO in self.VBOs]:
            self._VBOSize[value[0]] = int(value[1])
        else:
            raise ValueError('Invalid VBO id')

    def draw(self, mode):
        """

        :param mode: GL_LINES, GL_PATCH, ...
        :return:
        """
        if glGetIntegerv(GL_VERTEX_ARRAY_BINDING) is not self.VAO:
            glBindVertexArray(self.VAO)
        for i, VBO in enumerate(self.VBOs):
            glBindVertexBuffer(i, VBO['id'], 0, VBO['stride'])
        glDrawArrays(mode, 0, self.size)

    def addData(self, attr, data):
        """

        :param attr: VBO index in self.VBOs
        :param data: numpy array
        :return:
        """
        VBO = self.VBOs[attr]
        offset = VBO['size_bytes']
        self._validateData(VBO, data)
        self._checkForMemoryAllocation(VBO, data, offset)
        glNamedBufferSubData(VBO['id'], offset, data.nbytes, data)

        if VBO['data'] is not None:
            VBO['data'] = np.append(VBO['data'], data)
        else:
            VBO['data'] = data.flatten()
        VBO['size'] = VBO['data'].size
        self.size = [VBO['id'], VBO['size']/VBO['shape']]
        VBO['size_bytes'] = VBO['data'].nbytes

    def deleteData(self, attr, offset, bytes):
        """

        :param attr: VBO index in self.VBOs
        :param offset: offset in bytes
        :param bytes: num of bytes to delete
        :return:
        """
        VBO = self.VBOs[attr]
        data = np.empty(int(bytes/VBO['dtype_bytes']), dtype=VBO['data'].dtype.type)
        self._validateVBORange(VBO, data, offset)
        VBO['data'] = np.delete(VBO['data'], slice(int(offset/VBO['dtype_bytes']), int((offset+bytes)/VBO['dtype_bytes'])))
        VBO['size'] = VBO['data'].size
        VBO['size_bytes'] = VBO['data'].nbytes
        glNamedBufferSubData(VBO['id'], 0, VBO['data'].nbytes, VBO['data'])
        self.size = [VBO['id'], VBO['size']/VBO['shape']]

    def clearData(self, attr):
        VBO = self.VBOs[attr]
        VBO['data'] = None
        VBO['size'] = 0
        VBO['size_bytes'] = 0
        self.size = 0

    def replaceData(self, attr, data, offset):
        """

        :param attr: VBO index in self.VBOs
        :param data: numpy array
        :param offset: offset in bytes
        :return:
        """
        VBO = self.VBOs[attr]
        self._validateData(VBO, data)
        self._validateVBORange(VBO, data, offset)
        glNamedBufferSubData(VBO['id'], offset, data.nbytes, data)
        VBO['data'][int(offset/VBO['dtype_bytes']):int(offset/VBO['dtype_bytes']+len(data))] = data


class ShaderDataText(ShaderData):

    def draw(self, mode, textures=None):
        glActiveTexture(GL_TEXTURE0)
        if glGetIntegerv(GL_VERTEX_ARRAY_BINDING) is not self.VAO:
            glBindVertexArray(self.VAO)
        for i, VBO in enumerate(self.VBOs):
            glBindVertexBuffer(i, VBO['id'], 0, VBO['stride'])
        for i in range(0, self.size, 6):
            if textures is not None:
                glBindTexture(GL_TEXTURE_2D, textures[int(i/6)])
            glDrawArrays(mode, i, 6)


class ShaderDataCreator:
    """
    VAOFormat: List[List[count, type]]
    """

    DEFAULT_VBO_ALLOCATION_SIZE = 1024

    @classmethod
    def createShaderData(cls, VAOFormat, VBOAllocationSize=None, shader_data=ShaderData):
        if VBOAllocationSize is None:
            VBOAllocationSize = cls.DEFAULT_VBO_ALLOCATION_SIZE
        return shader_data(VAOFormat, VBOAllocationSize)


VAO_CREATORS = {}
VAO_CREATOR = VAOCreator()


def setCreator(context):
    global VAO_CREATORS
    creator = VAO_CREATORS.get(context, None)
    if creator is None:
        creator = VAOCreator()
    global VAO_CREATOR
    VAO_CREATOR = creator
