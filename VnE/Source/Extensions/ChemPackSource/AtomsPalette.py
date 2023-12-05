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


import json

import numpy as np


class Palette:

    def __init__(self):
        path = '\\'.join(__file__.split('\\')[:-1])
        self._palette = json.load(open(f'{path}\\palette.json', 'r'))
        self._palette_num = {}
        nums = []
        for key in self._palette:
            if self._palette[key][0] not in nums:
                self._palette_num[self._palette[key][0]] = [key] + self._palette[key][1:]
                nums.append(self._palette[key][0])
        self._point_dict = {}

    @property
    def point_dict(self):
        return self._point_dict

    @property
    def palette(self):
        return self._palette

    def addPoint(self, name, point):
        self._point_dict[name] = point

    def getColor(self, par):
        if type(par) is int:
            color = self._palette_num.get(par, [0, [0, 0, 0], '000000'])[1]
        elif type(par) is str and par.isdigit():
            color = self._palette_num.get(par, [0, [0, 0, 0], '000000'])[1]
        elif type(par) is list:
            color = np.zeros((3))
            for atom_type in par:
                color += self.getColor(atom_type)
            color = list(color/len(par))
        else:
            color = self._palette.get(par, [0, [0, 0, 0], '000000'])[1]
        return color

    def getName(self, par):
        if type(par) is int:
            color = self._palette_num.get(par, ['NON', [0, 0, 0], '000000'])[0]
        elif type(par) is str and par.isdigit():
            color = self._palette_num.get(int(par), ['NON', [0, 0, 0], '000000'])[0]
        else:
            color = self._palette.get(par, [999, [0, 0, 0], '000000'])[0]
        return color