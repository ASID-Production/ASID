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

import logging
import sys
import os
from django.conf import settings

level = logging.WARNING
if settings.DEBUG:
    level = logging.INFO


def set_logs():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s: %(levelname)s: %(message)s: %(name)s',
        handlers=[logging.StreamHandler(stream=sys.stdout)]
    )


def set_prm_log(prm):
    add_graphs_logger = logging.getLogger('add_graphs_logger')
    add_graphs_logger.setLevel(level)
    add_graphs_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', f'process_{str(prm)}.log'), mode="w")
    add_graphs_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    add_graphs_handler.setFormatter(add_graphs_formatter)
    add_graphs_logger.addHandler(add_graphs_handler)
    return add_graphs_logger


substructure_logger = logging.getLogger('substructure_logger')
substructure_logger.setLevel(level)
substructure_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', 'add_substructure_filtration.log'), mode='w')
substructure_formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(message)s: %(name)s')
substructure_handler.setFormatter(substructure_formatter)
substructure_logger.addHandler(substructure_handler)

all_cif_data_logger = logging.getLogger('all_cif_data_logger')
all_cif_data_logger.setLevel(level)
all_cif_data_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', 'add_all_cif_data.log'), mode='w')
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
all_cif_data_handler.setFormatter(formatter)
all_cif_data_logger.addHandler(all_cif_data_handler)

add_graphs_to_db_logger = logging.getLogger('add_graphs_to_db_logger')
add_graphs_to_db_logger.setLevel(level)
add_graphs_to_db_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', 'add_graphs_to_db.log'), mode='w')
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
add_graphs_to_db_handler.setFormatter(formatter)
add_graphs_to_db_logger.addHandler(add_graphs_to_db_handler)

vasp_logger = logging.getLogger('vasp_logger')
vasp_logger.setLevel(level)
vasp_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', 'vasp.log'), mode='w')
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
vasp_handler.setFormatter(formatter)
vasp_logger.addHandler(vasp_handler)

cif_db_update_main_logger = logging.getLogger('cif_db_update_main')
cif_db_update_main_logger.setLevel(level)
cif_db_update_main_handler = logging.FileHandler(os.path.join(settings.BASE_DIR, 'logs', 'cif_db_update.log'), mode='a')
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
cif_db_update_main_handler.setFormatter(formatter)
cif_db_update_main_logger.addHandler(cif_db_update_main_handler)
