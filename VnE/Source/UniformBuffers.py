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
import numpy as np
from typing import List
from OpenGL.GL import *


class aUniformBuffer(ABC):

    types = {np.byte: GL_BYTE,
             np.ubyte: GL_UNSIGNED_BYTE,
             np.intc: GL_INT,
             np.uintc: GL_UNSIGNED_INT,
             np.float32: GL_FLOAT}

    @abstractmethod
    def __init__(self):
        self.id: int
        self.size: int
        pass

    @abstractmethod
    def getInfo(self):
        pass

    def _validateVBORange(self, data, offset):
        if offset >= self.size:
            raise ValueError(f'Offset out of range, expected offset < {self.size}, got {offset}')
        last_pos = offset + data.nbytes
        if last_pos > self.size:
            raise ValueError(
                f"Replace location out of VBO range, expected from 0 to {self.size} bytes, got {last_pos} bytes")
        return

    def replaceData(self, data, offset):
        """

        :param data: numpy array
        :param offset: offset in bytes
        :return:
        """
        self._validateVBORange(data, offset)
        glNamedBufferSubData(self.id, offset, data.nbytes, data)

    def bind(self, binding_point=0):
        glBindBufferBase(GL_UNIFORM_BUFFER, binding_point, self.id)

    def getLine(self, property):
        value = self.__getattribute__(property)
        if issubclass(type(value), list) or issubclass(type(value), np.ndarray):
            data = [str(x) for x in value]
            value_line = '[' + ', '.join(data) + ']'
        else:
            value_line = str(value)
        return f'{property}: {value_line}'


class SceneUniformBuffer(aUniformBuffer):

    def __init__(self):
        self.properties = {'scale': 0,
                           'translation': 64,
                           'rotation': 128,
                           'aspect_ratio': 192,
                           'clip_distance': 256,
                           'perspective': 320,
                           'scene_shift': 384,
                           'wh': 448}

        self.info = {'scale': List[float],
                     'x_scale': float,
                     'y_scale': float,
                     'z_scale': float,
                     'translation': List[float],
                     'x_translation': float,
                     'y_translation': float,
                     'z_translation': float,
                     'rotation': List[float],
                     'x_rotation': float,
                     'y_rotation': float,
                     'z_rotation': float,
                     'rotation_point': List[float],
                     'rotation_point_x': float,
                     'rotation_point_y': float,
                     'rotation_point_z': float,
                     'aspect_ratio': List[float],
                     'perspective': List[float],
                     'left': float,
                     'right': float,
                     'top': float,
                     'bottom': float,
                     'near': float,
                     'far': float,
                     'wh': List[float],
                     'width': float,
                     'height': float,
                     'scene_shift': float}

        self.size = 456
        self._scale = np.array([[1.0, 0.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0, 0.0],
                               [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)

        self._translation = np.array([[1.0, 0.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0, 0.0],
                                     [0.0, 0.0, 1.0, 0.0],
                                     [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)

        self._rotation = np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]], dtype=np.float32)
        self._x_rot, self._y_rot, self._z_rot = np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]]), np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]]), np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]])
        self._x_rot_rad, self._y_rot_rad, self._z_rot_rad = 0, 0, 0
        self._rotation_point = np.array([0,0,0], dtype=np.float32)
        self._rotation_point_matr = np.array([[1.0, 0.0, 0.0, 0.0],
                                              [0.0, 1.0, 0.0, 0.0],
                                              [0.0, 0.0, 1.0, 0.0],
                                              [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)
        self._r_rotation_point_matr = np.array([[1.0, 0.0, 0.0, 0.0],
                                                [0.0, 1.0, 0.0, 0.0],
                                                [0.0, 0.0, 1.0, 0.0],
                                                [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)
        self._aspect_ratio = np.array([[1.0, 0.0, 0.0, 0.0],
                                  [0.0, 1.0, 0.0, 0.0],
                                  [0.0, 0.0, 1.0, 0.0],
                                  [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)

        self._clip_distance = np.array([[1.0, 0.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0, 0.0],
                               [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)

        l, r, t, b, n, f = -50.0, 50.0, 50.0, -50.0, 5, 10000
        self._perspective = np.array([[2 * n / (r - l), 0, (r + l) / (r - l), 0],
                                      [0, 2 * n / (t - b), (t + b) / (t - b), 0],
                                      [0, 0, -(f + n) / (f - n), -2 * f * n / (f - n)],
                                      [0, 0, -1, 0]]).astype(np.float32).transpose().astype(dtype=np.float32)
        self._left, self._right, self._top, self._bottom, self._near, self._far = l, r, t, b, n, f

        self._width, self._height = 240, 240
        self._wh = np.array([self._width, self._height], dtype=np.float32)

        self._scene_shift = np.array([[1.0, 0.0, 0.0, 0.0],
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, 1.0, 0.0],
                                      [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)


        self.id = glGenBuffers(1)
        glBindBuffer(GL_UNIFORM_BUFFER, self.id)
        glBufferData(GL_UNIFORM_BUFFER, self.size, np.array([self._scale.transpose(),
                                                             self._translation.transpose(),
                                                             self._rotation.transpose(),
                                                             self._aspect_ratio.transpose(),
                                                             self._clip_distance.transpose(),
                                                             self._perspective.transpose(),
                                                             self._scene_shift.transpose()], dtype=np.float32), GL_DYNAMIC_DRAW)
        glBufferSubData(GL_UNIFORM_BUFFER, self.properties['wh'], self._wh.nbytes, self._wh)
        glBindBuffer(GL_UNIFORM_BUFFER, 0)
        pass

    def replaceData(self, data, offset):
        data = data.astype(dtype=np.float32)
        super().replaceData(data, offset)

    def getInfo(self):
        return self.info

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, scale):
        if len(scale.shape) == 2:
            if scale.shape[0] == 4 and scale.shape[1] == 4:
                self._scale = scale
        elif len(scale) == 3:
            self._scale = np.array([[scale[0], 0.0, 0.0, 0.0],
                                    [0.0, scale[1], 0.0, 0.0],
                                    [0.0, 0.0, scale[2], 0.0],
                                    [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)
        scale = self._rotation_point_matr @ self._scale @ self._r_rotation_point_matr
        self.replaceData(scale.astype(dtype=np.float32).transpose(), self.properties['scale'])

    @property
    def x_scale(self):
        return self._scale[0, 0]

    @x_scale.setter
    def x_scale(self, scale):
        self._scale[0, 0] = scale
        scale = self._rotation_point_matr @ self._scale @ self._r_rotation_point_matr
        self.replaceData(scale.astype(dtype=np.float32).transpose(), self.properties['scale'])

    @property
    def y_scale(self):
        return self._scale[1, 1]

    @y_scale.setter
    def y_scale(self, scale):
        self._scale[1, 1] = scale
        scale = self._rotation_point_matr @ self._scale @ self._r_rotation_point_matr
        self.replaceData(scale.astype(dtype=np.float32).transpose(), self.properties['scale'])

    @property
    def z_scale(self):
        return self._scale[2, 2]

    @z_scale.setter
    def z_scale(self, scale):
        self._scale[2, 2] = scale
        scale = self._rotation_point_matr @ self._scale @ self._r_rotation_point_matr
        self.replaceData(scale.astype(dtype=np.float32).transpose(), self.properties['scale'])

    @property
    def translation(self):
        return self._translation

    @translation.setter
    def translation(self, transl):
        if len(transl.shape) == 2:
            if transl.shape[0] == 4 and transl.shape[1] == 4:
                self._translation = transl
        elif len(transl) == 3:
            self._translation = np.array([[1.0, 0.0, 0.0, transl[0]],
                                     [0.0, 1.0, 0.0, transl[1]],
                                     [0.0, 0.0, 1.0, transl[2]],
                                     [0.0, 0.0, 0.0, 1.0]]).astype(dtype=np.float32)
        self.replaceData(self._translation.transpose(), self.properties['translation'])

    @property
    def x_translation(self):
        return self._translation[0, 3]

    @x_translation.setter
    def x_translation(self, transl):
        self._translation[0, 3] = transl
        self.replaceData(self._translation.transpose(), self.properties['translation'])

    @property
    def y_translation(self):
        return self._translation[1, 3]

    @y_translation.setter
    def y_translation(self, transl):
        self._translation[1, 3] = transl
        self.replaceData(self._translation.transpose(), self.properties['translation'])

    @property
    def z_translation(self):
        return self._translation[2, 3]

    @z_translation.setter
    def z_translation(self, transl):
        self._translation[2, 3] = transl
        self.replaceData(self._translation.transpose(), self.properties['translation'])

    @property
    def rotation(self):
        return self._rotation

    @rotation.setter
    def rotation(self, rotation):
        if len(rotation.shape) == 2:
            if rotation.shape[0] == 4 and rotation.shape[1] == 4:
                self._rotation = rotation
        elif len(rotation) == 3:
            x, y, z = self._rotation_point
            self._x_rot = np.array([[1, 0, 0, 0],
                                    [0, np.cos(rotation[0]), -np.sin(rotation[0]), 0],
                                    [0, np.sin(rotation[0]), np.cos(rotation[0]), 0],
                                    [0, 0, 0, 1]])
            self._x_rot_rad = rotation[0] % (np.pi * 2)
            self._y_rot = np.array([[np.cos(rotation[1]), 0, -np.sin(rotation[1]), 0],
                                    [0, 1, 0, 0],
                                    [np.sin(rotation[1]), 0, np.cos(rotation[1]), 0],
                                    [0, 0, 0, 1]])
            self._y_rot_rad = rotation[1] % (np.pi * 2)
            self._z_rot = np.array([[np.cos(rotation[2]), -np.sin(rotation[2]), 0, 0],
                                    [np.sin(rotation[2]), np.cos(rotation[2]), 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 1]])
            self._z_rot_rad = rotation[2] % (np.pi * 2)
            self._rotation = self._z_rot @ self._y_rot @ self._x_rot @ self._rotation
            rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr
            self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])

    @property
    def x_rotation(self):
        return self._x_rot_rad

    @x_rotation.setter
    def x_rotation(self, rot):
        self._x_rot = np.array([[1, 0, 0, 0],
                                [0, np.cos(rot), -np.sin(rot), 0],
                                [0, np.sin(rot), np.cos(rot), 0],
                                [0, 0, 0, 1]])
        self._x_rot_rad = rot % (np.pi * 2)
        self._rotation = self._x_rot @ self._rotation
        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr
        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])

    @property
    def y_rotation(self):
        return self._y_rot_rad

    @y_rotation.setter
    def y_rotation(self, rot):
        self._y_rot = np.array(
            [[np.cos(rot), 0, -np.sin(rot), 0],
             [0, 1, 0, 0],
             [np.sin(rot), 0, np.cos(rot), 0],
             [0, 0, 0, 1]])
        self._y_rot_rad = rot % (np.pi * 2)
        self._rotation = self._y_rot @ self._rotation
        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr
        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])

    @property
    def z_rotation(self):
        return self._z_rot_rad

    @z_rotation.setter
    def z_rotation(self, rot):
        self._z_rot = np.array([[np.cos(rot), -np.sin(rot), 0, 0],
                                [np.sin(rot), np.cos(rot), 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]])
        self._z_rot_rad = rot % (np.pi * 2)
        self._rotation = self._z_rot @ self._rotation
        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr
        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])

    @property
    def rotation_point(self):
        return self._rotation_point

    @rotation_point.setter
    def rotation_point(self, coord):
        coord = np.array(coord, dtype=np.float32)
        self._rotation_point = coord
        self._rotation_point_matr[0, 3], self._r_rotation_point_matr[0, 3] = coord[0], -coord[0]
        self._rotation_point_matr[1, 3], self._r_rotation_point_matr[1, 3] = coord[1], -coord[1]
        self._rotation_point_matr[2, 3], self._r_rotation_point_matr[2, 3] = coord[2] + self.scene_shift, -(coord[2] + self.scene_shift)

        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr

        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])
        self.x_scale = self.x_scale

    @property
    def rotation_point_x(self):
        return self._rotation_point[0]

    @rotation_point_x.setter
    def rotation_point_x(self, x):
        self._rotation_point[0] = x
        self._rotation_point_matr[0, 3], self._r_rotation_point_matr[0, 3] = x, -x

        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr

        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])
        self.x_scale = self.x_scale

    @property
    def rotation_point_y(self):
        return self._rotation_point[1]

    @rotation_point_y.setter
    def rotation_point_y(self, y):
        self._rotation_point[1] = y

        self._rotation_point_matr[1, 3], self._r_rotation_point_matr[1, 3] = y, -y

        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr

        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])
        self.x_scale = self.x_scale

    @property
    def rotation_point_z(self):
        return self._rotation_point[2]

    @rotation_point_z.setter
    def rotation_point_z(self, z):
        self._rotation_point[2] = z

        self._rotation_point_matr[2, 3], self._r_rotation_point_matr[2, 3] = z + self.scene_shift, -(z + self.scene_shift)

        rot = self._rotation_point_matr @ self._rotation @ self._r_rotation_point_matr

        self.replaceData(rot.astype(dtype=np.float32).transpose(), self.properties['rotation'])
        self.x_scale = self.x_scale

    @property
    def aspect_ratio(self):
        return self._aspect_ratio

    @aspect_ratio.setter
    def aspect_ratio(self, wh):
        try:
            iter(wh)
            if len(wh) == 2:
                self._aspect_ratio[1, 1] = wh[0]/wh[1]
        except TypeError:
            self._aspect_ratio[1, 1] = wh
        self.replaceData(self._aspect_ratio, self.properties['aspect_ratio'])

    @property
    def perspective(self):
        return self._perspective

    @perspective.setter
    def perspective(self, persp):
        if len(persp.shape) == 2:
            if persp.shape[0] == 4 and persp.shape[1] == 4:
                self._perspective = persp
        elif len(persp) == 6:

            l, r, t, b, n, f = persp
            self._left, self._right, self._top, self._bottom, self._near, self._far = l, r, t, b, n, f
            self._perspective = np.array([[2*n / (r-l), 0, (r+l)/(r-l), 0],
                                          [0, 2*n / (t-b), (t+b)/(t-b), 0],
                                          [0, 0, -(f + n) / (f - n), -2 * f * n / (f - n)],
                                          [0, 0, -1, 0]]).astype(np.float32).astype(dtype=np.float32)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def left(self):
        return self._left

    @left.setter
    def left(self, l):
        self._left = l
        self._perspective[0, 0] = 2*self._near / (self._right - self._left)
        self._perspective[0, 2] = (self._right + self._left)/(self._right - self._left)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, r):
        self._right = r
        self._perspective[0, 0] = 2*self._near / (self._right - self._left)
        self._perspective[0, 2] = (self._right + self._left) / (self._right - self._left)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, t):
        self._top = t
        self._perspective[1, 1] = 2*self._near / (self._top - self._bottom)
        self._perspective[1, 2] = (self._top + self._bottom)/(self._top - self._bottom)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def bottom(self):
        return self._bottom

    @bottom.setter
    def bottom(self, b):
        self._bottom = b
        self._perspective[1, 1] = 2*self._near / (self._top - self._bottom)
        self._perspective[1, 2] = (self._top + self._bottom) / (self._top - self._bottom)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def near(self):
        return self._near

    @near.setter
    def near(self, n):
        self._near = n
        self._perspective[0, 0] = 2*self._near / (self._right - self._left)
        self._perspective[1, 1] = 2*self._near / (self._top - self._bottom)
        self._perspective[2, 2] = -(self._far + self._near) / (self._far - self._near)
        self._perspective[2, 3] = -2 * self._far * self._near / (self._far - self._near)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def far(self):
        return self._far

    @far.setter
    def far(self, f):
        self._far = f
        self._perspective[2, 2] = -(self._far + self._near) / (self._far - self._near)
        self._perspective[2, 3] = -2 * self._far * self._near / (self._far - self._near)
        self.replaceData(self._perspective.transpose(), self.properties['perspective'])

    @property
    def wh(self):
        return self._wh

    @wh.setter
    def wh(self, wh):
        self._width = wh[0]
        self._height = wh[1]
        self.aspect_ratio = self._width / self._height
        self._wh = np.array(wh, dtype=np.float32)
        self.replaceData(self._wh, self.properties['wh'])

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, width):
        self._width = width
        self.aspect_ratio = self._width / self._height
        self._wh = np.array([width, self._height], dtype=np.float32)
        self.replaceData(self._wh, self.properties['wh'])

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, height):
        self._height = height
        self.aspect_ratio = self._width / self._height
        self._wh = np.array([self._width, height], dtype=np.float32)
        self.replaceData(self._wh, self.properties['wh'])

    @property
    def scene_shift(self):
        return self._scene_shift[2, 3]

    @scene_shift.setter
    def scene_shift(self, shift):
        self._scene_shift[2, 3] = shift
        self.rotation_point = self.rotation_point
        self.replaceData(self._scene_shift.transpose(), self.properties['scene_shift'])