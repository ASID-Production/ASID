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
# *****************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# *****************************************************************************************

from rest_framework import serializers
import ast


class NodesListField(serializers.ListField):
    """
    Serialize such lists:
    ["(1, {'type': 'C', 'Hnum': 0, 'cord_min': 1, 'cord_max': 4})", "(2, {'type': 'Br', 'Hnum': 0, 'cord_min': 1, 'cord_max': 4})", "..."]
    to
    [(1, {'type': 'C', 'Hnum': 2, 'cord_min': 1, 'cord_max': 4}), (2, {'type': 'N', 'Hnum': 1, 'cord_min': 1, 'cord_max': 4}), ...]
    """
    def to_representation(self, data):
        """
        List of object instances -> List of dicts of primitive datatypes.
        """
        temp_data = []
        for element in data:
            element = list(ast.literal_eval(element))
            element[0] = int(element[0])
            element[1] = ast.literal_eval(str(element[1]))
            for v_key, v_value in element[1].items():
                if v_key in ['Hnum', 'cord_min', 'cord_max']:
                    element[1][v_key] = int(v_value)
            temp_data.append(tuple(element))
        data = temp_data.copy()
        return [item if item is not None else None for item in data]


class EdgesListField(serializers.ListField):
    """
    Serialize such lists:
    ["(1, 2, {'distance': 1.250,})", "(2, 3, {'distance': 1.340,})", ...]
    to
    [(1, 2, {'distance': 1.250,}), (2, 3, {'distance': 1.340,}), ...]
    """
    def to_representation(self, data):
        """
        List of object instances -> List of dicts of primitive datatypes.
        """
        temp_data = []
        for element in data:
            element = list(ast.literal_eval(element))
            element[0] = int(element[0])
            element[1] = int(element[1])
            element[2] = ast.literal_eval(str(element[2]))
            for v_key, v_value in element[2].items():
                if v_key == 'distance':
                    element[2][v_key] = float(v_value)
            temp_data.append(tuple(element))
        data = temp_data.copy()
        return [item if item is not None else None for item in data]
