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

from structure.models import CoordinatesBlock
from django_project.loggers import add_graphs_to_db_logger


def add_string_graph_to_db(graphs, refcode):
    add_graphs_to_db_logger.info(f'Write a graph string to the database')
    coord_block = CoordinatesBlock.objects.get(refcode=refcode)
    graphs = str(refcode.id) + ' ' + graphs
    coord_block.graph = graphs
    coord_block.save()


def upload_graphs_to_db(graph, structure_in_db):
    add_string_graph_to_db(graph['graph_str'], structure_in_db)
    add_graphs_to_db_logger.info(f'Structure {structure_in_db} successfully added to the database')
