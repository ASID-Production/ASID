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

import debug


class DefaultData:

    __default_data = {
        'name': 'X',
        'atom_type': 999,
        'coord': np.array([0,0,0], dtype=np.float32),
        'pdb_flag': 'HETATM',
        'pdb_atom_seq': 0,
        'pdb_name': 'X',
        'pdb_alt_loc': ' ',
        'pdb_res_name': 'UND',
        'pdb_chain': 'X',
        'pdb_res_seq': 0,
        'pdb_iCode': ' ',
        'pdb_occupancy': 1.00,
        'pdb_tempFactor': 0.00,
        'pdb_charge': '  ',
        'cif_space_group': 'P 1',
        'cif_sym_codes': [[0, 'x,y,z']],
        'cif_cell_a': -1,
        'cif_cell_b': -1,
        'cif_cell_c': -1,
        'cif_cell_al': 90,
        'cif_cell_be': 90,
        'cif_cell_ga': 90,
        'cif_frac_coords': np.array([0,0,0], dtype=np.float32),
    }

    def __getattr__(self, item):
        if item in DefaultData.__default_data:
            return DefaultData.__default_data[item]
        elif item.startswith('__') and item.endswith('__'):
            raise AttributeError
        else:
            return None


class aIter(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __next__(self):
        pass


class aEntity(ABC):

    def __init__(self, parent=None, name=None):
        self.children = []
        self._parent = parent
        self.name = name
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
            self.children.remove(child)

    def findProp(self, prop, value, res=None):
        if res is None:
            res = []
        try:
            if self.__getattribute__(prop) == value:
                res.append(self)
        except AttributeError:
            pass
        for child in self.children:
            res = child.findProp(prop, value, res=res)
        return res

    def updateProps(self, props: dict):
        for arg in props:
            self.__setattr__(arg, props[arg])


class Bond:
    def __init__(self, atom1, atom2):
        self._parents = [atom1, atom2]
        self._parents.sort(key=id)

    def __contains__(self, item):
        return self._parents.__contains__(item)

    def get(self, atom):
        return self._parents[self._parents.index(atom)-1]

    def changeParent(self, old_parent, parent):
        self._parents[self._parents.index(old_parent)] = parent
        self._parents.sort(key=id)

    def parents(self):
        return self._parents

    def __eq__(self, other):
        if isinstance(other, Bond):
            if self._parents == other.parents():
                return True
        else:
            return False


class Atom(DefaultData):
    def __init__(self, coord, atom_type: int|str|list, parent=None, **kwargs):
        self._parent = parent
        if parent is not None:
            parent.addChild(self)
        self._point = None
        if isinstance(coord, list) or isinstance(coord, tuple):
            coord = np.array(coord, dtype=np.float32)
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
        if any([bond == x for x in self._bonds]):
            return
        else:
            self._bonds.append(bond)

    def parent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def remove(self):
        if self.parent() is not None:
            self.parent().removeChild(self)
        self.__del__()

    def assignPoint(self, point):
        self._point = point

    def point(self):
        return self._point

    def findProp(self, prop, value, res=None):
        if res is None:
            res = []
        try:
            if self.__getattribute__(prop) == value:
                res.append(self)
        except AttributeError:
            pass
        return res

    def updateProps(self, props: dict):
        for arg in props:
            self.__setattr__(arg, props[arg])


class Molecule(aEntity):
    """Molecule class"""

    def __init__(self, parent=None):
        aEntity.__init__(self, parent=parent)
        self._points_list = None

    def addChild(self, child):
        if type(child) is Molecule:
            for atom in child:
                self.addChild(atom)
        elif type(child) is Atom:
            if child not in self.children:
                self.children.append(child)
                child.setParent(self)

    def assignPoint(self, point):
        self._points_list = point

    def point(self):
        return self._points_list

    def genBonds(self):
        import os
        import cpplib
        if self.__len__() < 1:
            return
        b = [(int(x.atom_type), *x.coord) for x in self.children]
        res = cpplib.GenBonds(b)['bonds']
        for ind in res:
            pair = [self.children[ind[0]], self.children[ind[1]]]
            bond = Bond(*pair)
            pair[0].addBond(bond)
            pair[1].addBond(bond)


class MoleculeSystem(aEntity):
    """Molecule system class, main class"""

    def __init__(self, parent=None):
        aEntity.__init__(self, parent=parent)

    def genBonds(self):
        for child in self.children:
            child.genBonds()