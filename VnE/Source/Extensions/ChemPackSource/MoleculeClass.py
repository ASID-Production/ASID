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
from ctypes import *
from ..ChemPack import PALETTE


class aIter(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __next__(self):
        pass


class aEntity(ABC):

    def __init__(self, parent=None):
        self.children = []
        self._parent = parent
        if self._parent is not None:
            self._parent.addChild(self)

    def __iter__(self):
        return iter(self.children)

    def __del__(self):
        return

    def __len__(self):
        return len(self.children)

    def addChild(self, child):
        if child not in self.children:
            self.children.append(child)

    def parent(self):
        return self._parent
        pass

    def remove(self):
        if self.parent() is not None:
            self.parent().removeChild(self)
        self.__del__()

    def removeChild(self, child):
        if child in self.children:
            self.children.pop(child)


class Bond:
    def __init__(self, atom1, atom2):
        self._parents = [atom1, atom2]

    def changeParent(self, pos, parent):
        self._parents[pos] = parent

    def parents(self):
        return self._parents


class Atom:
    def __init__(self, coord, atom_type: int, parent=None, **kwargs):
        self._parent = parent
        if parent is not None:
            parent.addChild(self)
        self._point = None
        self.coord = coord
        if type(atom_type) is list:
            self.atom_type = []
            for atype in atom_type:
                if type(atype) is str:
                    self.atom_type.append(PALETTE.getName(atype))
                else:
                    self.atom_type.append(atype)
        elif type(atom_type) is str:
            self.atom_type = PALETTE.getName(atom_type)
        else:
            self.atom_type = atom_type
        self._bonds = []
        for arg in kwargs:
            if arg not in self.__dict__:
                self.__setattr__(arg, kwargs[arg])

    def __del__(self):
        return

    def bonds(self):
        return self._bonds

    def addBond(self, bond):
        if bond not in self._bonds:
            self._bonds.append(bond)

    def parent(self):
        return self._parent

    def remove(self):
        if self.parent() is not None:
            self.parent().removeChild(self)
        self.__del__()

    def assignPoint(self, point):
        self._point = point

    def point(self):
        return self._point


class Molecule(aEntity):
    """Molecule class"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._points_list = None

    def addChild(self, child):
        if type(child) is Molecule:
            for atom in child:
                self.addChild(atom)
        elif type(child) is Atom:
            if child not in self.children:
                self.children.append(child)

    def assignPoint(self, point):
        self._points_list = point

    def point(self):
        return self._points_list

    def genBonds(self):
        if self.__len__() < 1:
            return
        b = [(c_float * 4)(*([c_float(a.atom_type)].__add__([c_float(coord) for coord in a.coord]))) for a in self.children]
        p = (POINTER(c_float) * len(b))(*b)
        direct = '\\'.join(__file__.split('\\')[:-1] + ['GenBonds.dll'])
        lib = CDLL(direct)
        lib.genBonds.restype = c_char_p
        lib.genBonds.argtypes = (c_uint, POINTER(POINTER(c_float)))
        line = lib.genBonds(c_uint(len(p)), p)
        line = line.decode().split('\n')
        line = [x for x in line if x]
        for i, pair in enumerate(line):
            pair = [self.children[int(x)] for x in pair.split(':')]
            bond = Bond(*pair)
            pair[0].addBond(bond)
            pair[1].addBond(bond)


class MoleculeSystem(aEntity):
    """Molecule system class, main class"""

    def __init__(self, parent=None):
        super().__init__(parent)

    def genBonds(self):
        for child in self.children:
            child.genBonds()