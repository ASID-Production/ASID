# Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
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
from typing import List, Dict, Any
import numpy as np

import debug


class CopyPointObserver:

    def __init__(self, owner):
        self._owner = owner
        self._subscribers = {}

    def update(self, object, property, value):
        if object is not self._owner:
            return
        subscribers = self._subscribers.get(property, None)
        if subscribers is not None:
            for subscriber in subscribers:
                subscriber.__setattr__(property, object)

    def addSubscriber(self, object, property):
        if self._subscribers.get(property, None) is not None:
            self._subscribers[property].append(object)
        else:
            self._subscribers[property] = [object]

    def removeSubscriber(self, object, property):
        try:
            self._subscribers[property].remove(object)
        except ValueError:
            return

    def add(self, object, *args, **kwargs):
        pass

    def remove(self, object, *args, **kwargs):
        for key in self._subscribers.keys():
            try:
                self.removeSubscriber(object, key)
            except ValueError:
                pass


class Signal(ABC):

    def __init__(self):
        self.methods = []

    @abstractmethod
    def emit(self, *args, **kwargs):
        for method in self.methods:
            method.__call__()

    def connect(self, method):
        self.methods.append(method)


class RemoveSeq(Signal):
    def emit(self, seq):
        for method in self.methods:
            method.__call__(seq)


class Sequence:

    def __init__(self):
        self._max_index = 0
        self._empty_index = []
        self.seqRemoved = RemoveSeq()

    def getIndex(self):
        if len(self._empty_index) > 0:
            index = self._empty_index.pop()
            return index
        else:
            index = self._max_index + 1
            self._max_index = index
            return index

    def removeIndex(self, index):
        if index not in self._empty_index:
            self._empty_index.append(index)
        self.seqRemoved.emit(index)


class aPoint(ABC):
    
    @abstractmethod
    def __init__(self, parent=None):

        super().__setattr__('_destroyed', False)
        self._not_observed_properties: List[str, ...] = ['parent', 'children', 'observers', 'dependency_observer']
        self._dependent_properties = {}
        self._properties: Dict[str: Any] = {'seq': self.seq}
        self._dependency_observer = CopyPointObserver(self)
        self._seq = sequence.getIndex()
        self.parent = parent
        if self.parent is not None:
            self.parent.addChild(self)
        self.children: List = []
        self.observers: List = [self._dependency_observer]

    @abstractmethod
    def copy(self):
        pass

    @property
    def dependency_observer(self):
        return self._dependency_observer

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, value):
        return

    def isValid(self):
        return not self._destroyed

    def getBySeq(self, seq):
        if self._seq == seq:
            return self
        else:
            for child in self.children:
                ret = child.getBySeq(seq)
                if ret is not None:
                    return ret
        return None

    def destroy(self, point=None):
        if point:
            for observer in self.observers.copy():
                observer.remove(point)
            try:
                self.children.remove(point)
            except ValueError:
                pass
            if self.parent is not None:
                self.parent.destroy(point=point)
        else:
            for child in self.children.copy():
                child.destroy()
            for observer in self.observers.copy():
                self.detach(observer)
            if self.parent is not None:
                self.parent.destroy(point=self)
            sequence.removeIndex(self.seq)
            self._destroyed = True

    def attach(self, observer):
        if observer not in self.observers:
            self.observers.append(observer)
            observer.add(self)

    def detach(self, observer):
        if observer in self.observers:
            self.observers.remove(observer)
            observer.remove(self)

    def notify(self, property, value, **kwargs):
        point = kwargs.get('point', self)
        for observer in self.observers:
            observer.update(point, property, value)
        if self.parent is not None:
            self.parent.notify(property, value, point=point)

    def addChild(self, child, branch=False):
        if not branch:
            if child not in self.children:
                for observer in self.observers:
                    if observer is not self._dependency_observer:
                        observer.add(child)
                self.children.append(child)
                if self.parent is not None:
                    self.parent.addChild(child, branch=True)
        else:
            for observer in self.observers:
                if observer is not self._dependency_observer:
                    observer.add(child)
            if self.parent is not None:
                self.parent.addChild(child, branch=True)

    def setParent(self, parent):
        self.parent = parent

    def addProperty(self, property, value, observed=True):
        if not observed:
            if property not in self._not_observed_properties:
                self._not_observed_properties.append(property)
        self._properties[property] = value
        self.__setattr__(property, value)

    def deleteProperty(self, property):
        if property == 'seq':
            return
        if property in self._not_observed_properties:
            self._not_observed_properties.remove(property)
        del self._properties[property]
        self.__delattr__(property)

    def getProperties(self):
        return self._properties

    def getLine(self, property):
        value = self.__getattribute__(property)
        if property in self._dependent_properties:
            value_line = str(self._dependent_properties[property])
        elif issubclass(type(value), list) or issubclass(type(value), np.ndarray):
            data = [str(x) for x in value]
            value_line = '[' + ', '.join(data) + ']'
        else:
            value_line = str(value)
        return f'{property}: {value_line}'

    def select(self, pos, tol=1, selection=None, mod=None, obs=False):
        if (len(self.observers) > 1 and self.coord is not None) or (obs and self.coord is not None):
            if mod is not None:
                coord = mod @ np.append(self.coord, 1)[:, np.newaxis]
                coord = coord.flatten()
                tol_curr = tol/coord[3]
                coord = coord/coord[3]
            else:
                coord = self.coord
            if self.rad is None:
                rad1 = 1
            else:
                rad1 = self.rad
            if np.all(np.abs(coord[:2]-pos) < tol_curr * rad1):
                if selection is None:
                    selection = [self, coord[2]]
                else:
                    if coord[2] < selection[1]:
                        selection = [self, coord[2]]
                    elif coord[2] == selection[1]:
                        if selection[0].rad is None:
                            rad2 = 1
                        else:
                            rad2 = selection[0].rad
                        if rad1 > rad2:
                            selection = [self, coord[2]]

        for child in self.children:
            selection = child.select(pos, tol, selection, mod)
        return selection

    def __del__(self):
        self.destroy()

    def __str__(self):
        return f'Point {self._seq}'

    def __delattr__(self, item):
        self.__setattr__(item, None)
        super().__delattr__(item)

    def __setattr__(self, key, value):
        if issubclass(type(value), aPoint) and key != 'parent':
            if self._dependent_properties.get(key, None) is value:
                pass
            elif self._dependent_properties.get(key, None) is None:
                value.dependency_observer.addSubscriber(self, key)
                self._dependent_properties[key] = value
            else:
                value.dependency_observer.addSubscriber(self, key)
                self._dependent_properties[key].dependency_observer.removeSubscriber(self, key)
                self._dependent_properties[key] = value
            if self._properties is not None and key in self._properties.keys():
                self._properties[key] = value
            value = value.__getattribute__(key)
        elif (key[0] != '_' and key not in self._not_observed_properties) and key in self._dependent_properties and key != 'parent':
            self._dependent_properties[key].dependency_observer.removeSubscriber(self, key)
            self._dependent_properties.pop(key)
        if key[0] != '_' and key not in self._not_observed_properties:
            super().__setattr__(key, value)
            if self._properties is not None and key in self._properties.keys() and key not in self._dependent_properties:
                self._properties[key] = value
            self.notify(key, self.__getattribute__(key))
        else:
            if self._properties is not None and key in self._properties.keys() and key not in self._dependent_properties:
                self._properties[key] = value
            super().__setattr__(key, value)

    def __getattribute__(self, item):
        try:
            value = super().__getattribute__(item)
        except AttributeError:
            value = None
        if issubclass(type(value), aPoint) and item != 'parent':
            value = value.__getattribute__(item)
        return value


class Point(aPoint):

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent=parent)

        properties = {key: kwargs[key] for key in kwargs}
        self._properties.update(properties)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])

    def copy(self, obs_clone=True):
        clone, linked_props, clones = self._copy(obs_clone)
        for source in linked_props.keys():
            source_clone = clones.get(source, None)
            if source_clone:
                for consumer in linked_props[source].keys():
                    for prop in linked_props[source][consumer]:
                        consumer.__setattr__(prop, source_clone)
        return clone

    def _copy(self, obs_clone=True, linked_props=None, parent=None, clones=None):
        if clones is None:
            clones = {}
        if linked_props is None:
            linked_props = {}
        props = {key: val for key, val in self._properties.items() if key != 'seq'}
        for key, val in props.items():
            if issubclass(type(val), aPoint):
                if linked_props.get(val, None) is None:
                    linked_props[val] = {self: [key]}
                else:
                    if linked_props[val].get(self, None) is None:
                        linked_props[val][self] = [key]
                    else:
                        if key not in linked_props[val][self]:
                            linked_props[val][self].append(key)

        clone = type(self)(parent=parent, **props)
        clones[self] = clone

        if obs_clone:
            for observer in self.observers:
                if observer is not self.dependency_observer:
                    clone.attach(observer)
        return clone, linked_props, clones

    def linkCopy(self):
        pass

    @property
    def coord(self):
        return self._coord

    @coord.setter
    def coord(self, coord):
        if coord is None:
            self._coord = None
            return
        self._coord = np.array(coord, dtype=np.float32)

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, color):
        if color is None:
            self._color = None
            return
        self._color = np.array(color, dtype=np.float32)

    @property
    def rad(self):
        return self._rad

    @rad.setter
    def rad(self, rad):
        if rad is None:
            self._rad = None
            return
        try:
            iter(rad)
        except TypeError:
            rad = np.array([rad], dtype=np.float32)

        self._rad = rad


class PointsList(aPoint):
    
    def __init__(self, parent=None, **kwargs):
        super().__init__(parent=parent)

        properties = {key: kwargs[key] for key in kwargs}
        self._properties.update(properties)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])

    def addChild(self, child, branch=False):
        if issubclass(type(child), aPoint):
            super().addChild(child, branch)
            child.setParent(self)

    def copy(self, obs_clone=True):
        clone, linked_props, clones = self._copy(obs_clone)
        for source in linked_props.keys():
            source_clone = clones.get(source, None)
            if source_clone:
                for consumer in linked_props[source].keys():
                    clone_consumer = clones[consumer]
                    for prop in linked_props[source][consumer]:
                        clone_consumer.__setattr__(prop, source_clone)
        return clone

    def _copy(self, obs_clone=True, linked_props=None, parent=None, clones=None):
        if clones is None:
            clones = {}
        if linked_props is None:
            linked_props = {}
        props = {key: val for key, val in self._properties.items() if key != 'seq'}
        for key, val in props.items():
            if issubclass(type(val), aPoint):
                if linked_props.get(val, None) is None:
                    linked_props[val] = {self: [key]}
                else:
                    if linked_props[val].get(self, None) is None:
                        linked_props[val][self] = [key]
                    else:
                        if key not in linked_props[val][self]:
                            linked_props[val][self].append(key)

        clone = type(self)(parent=parent, **props)
        clones[self] = clone

        if obs_clone:
            for observer in self.observers:
                if observer is not self.dependency_observer:
                    clone.attach(observer)

        for child in self.children:
            _, linked_props, clones = child._copy(obs_clone, linked_props, parent=clone, clones=clones)
        return clone, linked_props, clones

    def linkCopy(self):
        pass

    def changeOrder(self, object, index):
        if object not in self.children:
            return
        self.children.remove(object)
        if index > len(self.children):
            self.children.append(object)
        else:
            self.children.insert(index, object)
        self.notify('order', self.children, point=object)

    def select(self, pos, tol=1, selection=None, mod=None, obs=False):
        if len(self.observers) > 1:
            obs = True
        for child in self.children:
            selection = child.select(pos, tol, selection, mod, obs=obs)
        return selection

    def __str__(self):
        return f'PointsList {self._seq}'


sequence = Sequence()