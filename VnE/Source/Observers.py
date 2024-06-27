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
from . import text_render
from .point_class import Point, PointsList

import debug


class SingleObserver:

    def __init__(self, facade, scene):
        self.scene = scene
        self.facade = facade
        self.obs_dict = {}

    def getObserver(self, observer_cls):
        observer = self.obs_dict.get(observer_cls, None)
        if observer:
            return observer
        else:
            observer = observer_cls(self.facade, self.scene)
            self.obs_dict[observer_cls] = observer
            return observer


class aObserver(ABC):

    @abstractmethod
    def __init__(self, facade, scene):
        self._pipeline = None
        pass

    @abstractmethod
    def update(self, object, property, value):
        pass

    @abstractmethod
    def add(self, object, *args, **kwargs):
        pass

    @abstractmethod
    def remove(self, object, *args, **kwargs):
        pass

    def getPipeline(self):
        return self._pipeline


class SphereObserver(aObserver):

    def __init__(self, facade, scene):

        from .ShaderPipelines import BallsShaderPipeline

        self._scene = scene
        self._facade = facade
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=BallsShaderPipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._points = []
        self._properties_list = ['coord', 'color', 'rad', 'pick']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'color': np.array([0, 0, 0, 1], dtype=np.float32),
                            'rad': np.array([1], dtype=np.float32),
                            'pick': np.array([0], dtype=np.float32)}

    def update(self, object, property, value):
        if type(object).__name__ is PointsList.__name__:
            return
        if object not in self._points:
            return
        if property not in self._properties_list:
            return
        if value is None:
            value = self._properties[property]
        if isinstance(value, float) or isinstance(value, int):
            value = np.array([value], dtype=np.float32)
        self._facade.replaceDataInShaderData(self._shader_data, self._properties_list.index(property), value, self._points.index(object) * self._properties[property].nbytes)

    def add(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            if object.children is None:
                return
            for point in object.children:
                self.add(point)
        elif type(object).__name__ is Point.__name__:
            data = {}

            for property in self._properties:
                try:
                    data[property] = object.__getattribute__(property)
                    if isinstance(data[property], float) or isinstance(data[property], int):
                        data[property] = np.array([data[property]], dtype=np.float32)
                except AttributeError:
                    data[property] = self._properties[property]
                if data[property] is None:
                    data[property] = self._properties[property]
            for property in self._properties:
                self._facade.addDataToShaderData(self._shader_data, self._properties_list.index(property), data[property])
            self._points.append(object)

    def remove(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            for point in object.children:
                self.remove(point)
        elif type(object).__name__ is Point.__name__:
            try:
                pos = self._points.index(object)
            except ValueError:
                return
            self._points.remove(object)
            for property in self._properties_list:
                self._facade.deleteDataInShaderData(self._shader_data, self._properties_list.index(property), pos * self._properties[property].nbytes, self._properties[property].nbytes)


class BondsObserver(aObserver):

    def __init__(self, facade, scene):

        from .ShaderPipelines import BondShaderPipeline

        self._scene = scene
        self._facade = facade
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=BondShaderPipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._points = []
        self._properties_list = ['coord', 'color', 'rad', 'pick']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'color': np.array([0, 0, 0, 1], dtype=np.float32),
                            'rad': np.array([1], dtype=np.float32),
                            'pick': np.array([0], dtype=np.float32)}

    def changeOrder(self, object, order):
        old_pos = self._points.index(object)
        new_pos = order.index(object)
        shader_data = self._facade.getInst(self._shader_data)
        if new_pos > old_pos:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(old_pos + 1) * VBO['shape']:(new_pos + 1) * VBO['shape']]
                data = np.append(data, shift_point_data)
                VBO['data'][old_pos * VBO['shape']:old_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, old_pos * VBO['stride'])
        else:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(new_pos) * VBO['shape']:(old_pos) * VBO['shape']]
                data = np.append(shift_point_data, data)
                VBO['data'][new_pos * VBO['shape']:new_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, new_pos * VBO['stride'])
        self._points.remove(object)
        self._points.insert(new_pos, object)

    def update(self, object, property, value):
        if property == 'order':
            self.changeOrder(object, value)
        if type(object).__name__ is PointsList.__name__:
            return
        if object not in self._points:
            return
        if property not in self._properties_list:
            return
        if value is None:
            value = self._properties[property]
        if isinstance(value, float) or isinstance(value, int):
            value = np.array([value], dtype=np.float32)
        self._facade.replaceDataInShaderData(self._shader_data, self._properties_list.index(property), value, self._points.index(object) * self._properties[property].nbytes)

    def add(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            if object.children is None:
                return
            for point in object.children:
                self.add(point)
        elif type(object).__name__ is Point.__name__:
            data = {}

            for property in self._properties:
                try:
                    data[property] = object.__getattribute__(property)
                    if isinstance(data[property], float) or isinstance(data[property], int):
                        data[property] = np.array([data[property]], dtype=np.float32)
                except AttributeError:
                    data[property] = self._properties[property]
                if data[property] is None:
                    data[property] = self._properties[property]
            for property in self._properties:
                self._facade.addDataToShaderData(self._shader_data, self._properties_list.index(property), data[property])
            self._points.append(object)

    def remove(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            for point in object.children:
                self.remove(point)
        elif type(object).__name__ is Point.__name__:
            try:
                pos = self._points.index(object)
            except ValueError:
                return
            self._points.remove(object)
            for property in self._properties_list:
                self._facade.deleteDataInShaderData(self._shader_data, self._properties_list.index(property), pos * self._properties[property].nbytes, self._properties[property].nbytes)


class LabelObserver(aObserver):

    def __init__(self, facade, scene, wh=None):
        if wh is None:
            wh = [1, 1]
        from .ShaderPipelines import TextShaderPipeline
        self._font = text_render.Font()
        self._wh = wh

        self._facade = facade

        self._last_pos = 0
        self._scene = scene
        self._points = []
        self._points_list = []
        self._textures = []
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=TextShaderPipeline)
        self._pipeline_inst = self._facade.getInst(self._pipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._properties_list = ['coord', 'label_color', 'label_size', 'label']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'label_color': np.array([0, 0, 0, 1], dtype=np.float32),
                            'label_size': np.array([6], dtype=np.float32),
                            'label': None}
        self._point = None

    def indexPoint(self, point):
        for i, p in enumerate(self._points):
            if p['point'] is point:
                return i
        raise ValueError(f'{point} not in list')

    def set_wh(self, wh):
        self._wh = wh

    def add(self, object, *args, **kwargs):

        if type(object).__name__ is PointsList.__name__:
            if object.children is None:
                return
            for point in object.children:
                self.add(point)
        elif type(object).__name__ is Point.__name__:
            if object.label is None or object.label == '':
                self._points.append({'point': object, 'pos': -1, 'word': object.label})
                return
            else:
                self._points.append({'point': object, 'pos': self._last_pos, 'word': object.label})
                self._last_pos += len(object.label) * 6
                self._points_list.append(object)
                coord = object.coord
                word = object.label

                size = object.label_size
                if size is None:
                    size = int(self._properties['label_size'][0])
                word_obj = text_render.Word(coord, self._wh, self._font, word, size)
                data = word_obj.gen_buffer()
                textures = word_obj.get_textures()
                self._textures += textures
                self._pipeline_inst.textures = self._textures
                coord = np.array([x[:3] for x in data], dtype=np.float32).flatten()
                texCoord = np.array([x[3:5] for x in data], dtype=np.float32).flatten()
                size = np.array([x[5:7] for x in data], dtype=np.float32).flatten()
                shift = np.array([x[7:] for x in data], dtype=np.float32).flatten()
                self._facade.addDataToShaderData(self._shader_data, 0, coord)
                self._facade.addDataToShaderData(self._shader_data, 1, texCoord)
                self._facade.addDataToShaderData(self._shader_data, 2, size)
                self._facade.addDataToShaderData(self._shader_data, 3, shift)

    def remove(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            for point in object.children:
                self.remove(point)
        elif type(object).__name__ is Point.__name__:
            try:
                ind = self.indexPoint(object)
            except ValueError:
                return
            if self._points[ind]['word'] is None:
                self._points.pop(ind)
                return
            pos = self._points[ind]['pos']
            for point_ind in range(ind, len(self._points)):
                self._points[point_ind]['pos'] -= len(self._points[ind]['word']) * 6
            del self._textures[int(pos/6):int(pos/6 + len(self._points[ind]['word']))]
            self._pipeline_inst.textures = self._textures
            self._facade.deleteDataInShaderData(self._shader_data, 0, pos * self._properties['coord'].nbytes, self._properties['coord'].nbytes * len(self._points[ind]['word']) * 6)
            self._facade.deleteDataInShaderData(self._shader_data, 1, pos * 2 * 4, 2 * 4 * len(self._points[ind]['word']) * 6)
            self._facade.deleteDataInShaderData(self._shader_data, 2, pos * 2 * 4, 2 * 4 * len(self._points[ind]['word']) * 6)
            self._facade.deleteDataInShaderData(self._shader_data, 3, pos * 2 * 4, 2 * 4 * len(self._points[ind]['word']) * 6)
            self._last_pos -= len(self._points[ind]['word']) * 6
            self._points.pop(ind)

    def update(self, object, property, value):
        if type(object).__name__ is PointsList.__name__:
            return
        try:
            ind = self.indexPoint(object)
        except ValueError:
            return
        if property not in self._properties_list:
            return
        if value is None:
            value = self._properties[property]
        if isinstance(value, float) or isinstance(value, int):
            value = np.array([value], dtype=np.float32)
        if property == 'label_size':
            self._font.change_size(int(value[0]))
        self.remove(object)
        self.add(object)


class PlaneObserver(aObserver):
    def __init__(self, facade, scene):

        from .ShaderPipelines import PlaneShaderPipeline

        self._scene = scene
        self._facade = facade
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=PlaneShaderPipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._points = []
        self._properties_list = ['coord', 'color']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'color': np.array([0, 0, 0, 1], dtype=np.float32)}

    def changeOrder(self, object, order):
        old_pos = self._points.index(object)
        new_pos = order.index(object)
        shader_data = self._facade.getInst(self._shader_data)
        if new_pos > old_pos:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(old_pos + 1) * VBO['shape']:(new_pos + 1) * VBO['shape']]
                data = np.append(data, shift_point_data)
                VBO['data'][old_pos * VBO['shape']:old_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, old_pos * VBO['stride'])
        else:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(new_pos) * VBO['shape']:(old_pos) * VBO['shape']]
                data = np.append(shift_point_data, data)
                VBO['data'][new_pos * VBO['shape']:new_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, new_pos * VBO['stride'])
        self._points.remove(object)
        self._points.insert(new_pos, object)

    def update(self, object, property, value):
        if property == 'order':
            self.changeOrder(object, value)
        if type(object).__name__ is PointsList.__name__:
            return
        if object not in self._points:
            return
        if property not in self._properties_list:
            return
        if value is None:
            value = self._properties[property]
        if isinstance(value, float) or isinstance(value, int):
            value = np.array([value], dtype=np.float32)
        self._facade.replaceDataInShaderData(self._shader_data, self._properties_list.index(property), value, self._points.index(object) * self._properties[property].nbytes)

    def add(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            if object.children is None:
                return
            for point in object.children:
                self.add(point)
        elif type(object).__name__ is Point.__name__:
            data = {}

            for property in self._properties:
                try:
                    data[property] = object.__getattribute__(property)
                except AttributeError:
                    data[property] = self._properties[property]
                if data[property] is None:
                    data[property] = self._properties[property]
            for property in self._properties:
                self._facade.addDataToShaderData(self._shader_data, self._properties_list.index(property), data[property])
            self._points.append(object)

    def remove(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            for point in object.children:
                self.remove(point)
        elif type(object).__name__ is Point.__name__:
            if object in self._points:
                pos = self._points.index(object)
                self._points.remove(object)
                for property in self._properties_list:
                    self._facade.deleteDataInShaderData(self._shader_data, self._properties_list.index(property), pos * self._properties[property].nbytes, self._properties[property].nbytes)


class LineObserver(aObserver):

    def __init__(self, facade, scene):
        from .ShaderPipelines import LinesShaderPipeline

        self._scene = scene
        self._facade = facade
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=LinesShaderPipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._points = []
        self._properties_list = ['coord', 'color', 'pick']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'color': np.array([0, 0, 0, 1], dtype=np.float32),
                            'pick': np.array([0], dtype=np.float32)}

    def changeOrder(self, object, order):
        old_pos = self._points.index(object)
        new_pos = order.index(object)
        shader_data = self._facade.getInst(self._shader_data)
        if new_pos > old_pos:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(old_pos + 1) * VBO['shape']:(new_pos + 1) * VBO['shape']]
                data = np.append(data, shift_point_data)
                VBO['data'][old_pos * VBO['shape']:old_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, old_pos * VBO['stride'])
        else:
            for attr, property in enumerate(self._properties_list):
                VBO = shader_data.VBOs[attr]
                shift_point_data = VBO['data'][old_pos * VBO['shape']:(old_pos + 1) * VBO['shape']]
                data = VBO['data'][(new_pos) * VBO['shape']:(old_pos) * VBO['shape']]
                data = np.append(shift_point_data, data)
                VBO['data'][new_pos * VBO['shape']:new_pos * VBO['shape'] + data.size] = data
                self._facade.replaceDataInShaderData(self._shader_data, attr, data, new_pos * VBO['stride'])
        self._points.remove(object)
        self._points.insert(new_pos, object)

    def update(self, object, property, value):
        if property == 'order':
            self.changeOrder(object, value)
        if type(object).__name__ is PointsList.__name__:
            return
        if object not in self._points:
            return
        if property not in self._properties_list:
            return
        if value is None:
            value = self._properties[property]
        if isinstance(value, float) or isinstance(value, int):
            value = np.array([value], dtype=np.float32)
        self._facade.replaceDataInShaderData(self._shader_data, self._properties_list.index(property), value, self._points.index(object) * self._properties[property].nbytes)

    def add(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            if object.children is None:
                return
            for point in object.children:
                self.add(point)
        elif type(object).__name__ is Point.__name__:
            data = {}

            for property in self._properties:
                try:
                    data[property] = object.__getattribute__(property)
                    if isinstance(data[property], float) or isinstance(data[property], int):
                        data[property] = np.array([data[property]], dtype=np.float32)
                except AttributeError:
                    data[property] = self._properties[property]
                if data[property] is None:
                    data[property] = self._properties[property]
            for property in self._properties:
                self._facade.addDataToShaderData(self._shader_data, self._properties_list.index(property), data[property])
            self._points.append(object)

    def remove(self, object, *args, **kwargs):
        if type(object).__name__ is PointsList.__name__:
            for point in object.children:
                self.remove(point)
        elif type(object).__name__ is Point.__name__:
            try:
                pos = self._points.index(object)
            except ValueError:
                return
            self._points.remove(object)
            for property in self._properties_list:
                self._facade.deleteDataInShaderData(self._shader_data, self._properties_list.index(property), pos * self._properties[property].nbytes, self._properties[property].nbytes)


class DashedLineObserver(LineObserver):

    def __init__(self, facade, scene):
        from .ShaderPipelines import DashedLineShaderPipeline

        self._scene = scene
        self._facade = facade
        self._pipeline = self._facade.addPipelineToScene(self._scene, pipeline_cls=DashedLineShaderPipeline)
        self._shader_data = self._facade.addDataBufferToPipeline(self._pipeline)
        self._points = []
        self._properties_list = ['coord', 'color', 'pick']
        self._properties = {'coord': np.array([0, 0, 0], dtype=np.float32),
                            'color': np.array([0, 0, 0, 1], dtype=np.float32),
                            'pick': np.array([0], dtype=np.float32)}


observers = {'Sphere': SphereObserver,
             'Bond': BondsObserver,
             'Label': LabelObserver,
             'Plane': PlaneObserver,
             'Line': LineObserver,
             'Dashed line': DashedLineObserver}

