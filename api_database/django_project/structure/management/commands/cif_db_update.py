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

import time

from django.core.management.base import BaseCommand
import os
from .cif_db_update_modules._cifparser import add_coords, add_cell_parms_with_error, add_other_info
from .cif_db_update_modules._make_graphs_c import add_graphs_c
from .cif_db_update_modules._add_graphs_to_db import upload_graphs_to_db
from .cif_db_update_modules._add_substructure_filtration import add_substructure_filters
from .cif_db_update_modules._add_all_cif_data import add_all_cif_data
from CifFile import ReadCif
from structure.models import StructureCode
import multiprocessing
from django_project.loggers import cif_db_update_main_logger as logger_main
import codecs

NUM_OF_PROC = int(multiprocessing.cpu_count() / 2)  # number of physical processors
MAX_TIME_WAIT = 600  # maximum time to wait for process completion (sec)
CHUNK_SIZE = 1000  # size of the processed part of the array


def collect_cif_data(file: str, cif_blocks: dict, user_refcode=''):
    # read cif file
    try:
        cif = ReadCif(file)
    except:
        with codecs.open(file, 'r', encoding='ansi') as f1:
            lines = f1.read()
        with codecs.open(file, 'w', encoding='utf8') as f2:
            f2.write(lines)
        try:
            cif = ReadCif(file)
        except Exception as err:
            logger_main.error(f"Exception: Failed to read cif file!", exc_info=True)
            return cif_blocks
    # look through each structural block in cif file
    for block in cif.items():
        db = ''
        if user_refcode:
            refcode = user_refcode
        elif '_database_code_icsd' in block[1].keys():
            refcode = 'ICSD_' + str(block[1]['_database_code_icsd'])
            db = 'ICSD'
        elif '_database_code_csd' in block[1].keys():
            refcode = block[1]['_database_code_csd']
        elif '_cod_database_code' in block[1].keys():
            refcode = 'COD_' + str(block[1]['_cod_database_code'])
            db = 'COD'
        else:
            raise Exception(f'No refcode was found in cif file or in input parameters!')
        # if it belongs to another database, then we record the information
        if db:
            str_obj, created = StructureCode.objects.get_or_create(refcode=refcode)
            setattr(str_obj, db, True)
            str_obj.save()
        cif_blocks[refcode] = block
        return cif_blocks


def get_files(args):
    files = []
    for arg in args:
        if arg == 'all_data':
            continue
        if os.path.exists(arg):
            if arg.endswith('.cif'):
                files.append(arg)
            else:
                # we get a list of cif files in the directory
                all_files = os.listdir(arg)
                for file in all_files:
                    if file.endswith('.cif'):
                        files.append(os.path.join(arg, file))
        else:
            raise FileNotFoundError(f'File {arg} does not exist!')
    return files


def manager_collect_cifs(files, user_refcodes: dict):
    cif_blocks: dict = {}
    for file in files:
        if user_refcodes and file in user_refcodes.keys():
            user_refcode = user_refcodes[file]
        else:
            user_refcode = ''
        collect_cif_data(file, cif_blocks, user_refcode)
    return cif_blocks


def manager_add_coords_and_params_to_db(cif_blocks):
    for refcode, cif_block in cif_blocks.items():
        structure, created = StructureCode.objects.get_or_create(refcode=refcode)
        atoms = add_coords(cif_block, structure)
        add_cell_parms_with_error(cif_block, structure)
        add_other_info(atoms, structure)


def create_queue(cif_blocks: dict):
    queue = multiprocessing.Queue()
    refcodes_to_graph = []
    for refcode, cif_block in cif_blocks.items():
        if {'_atom_site_label', '_atom_site_type_symbol',
            '_atom_site_fract_x', '_atom_site_fract_y',
            '_atom_site_fract_z'}.issubset(cif_block[1].keys()
                                           ):
            queue.put([refcode, cif_block])
            refcodes_to_graph.append(refcode)
    return queue, refcodes_to_graph


def create_graph_c(queue):
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    processes = []
    procs = NUM_OF_PROC
    # if the queue is small, then we process it in one thread
    if queue.qsize() < 20:
        procs = 1
    for i in range(procs):
        process = multiprocessing.Process(target=add_graphs_c, args=(queue, return_dict, i + 1))
        process.start()
        processes.append(process)
    # monitor the progress of the task
    while not queue.empty():
        time.sleep(2)

    # waiting for all processes to complete
    start = time.time()
    while True:
        time.sleep(1)
        count = 0
        for process in processes:
            if not process.is_alive():
                count += 1
        # if all processes have already completed
        if count == len(processes):
            break
        # if the waiting time is exceeded, we interrupt the processes
        if time.time() - start > 180:
            for process in processes:
                for i in range(1, 101, 1):
                    if process.is_alive():
                        # wait for the process to complete for another 1 minute
                        time.sleep(60)
                    else:
                        break
                    if i == 100:
                        logger_main.warning(f'Thread {process.name} was force terminated: thread termination timed out!')
                        process.terminate()
            break
    return return_dict


def manager_upload_graphs_to_db(graphs):
    for refcode, graphs in graphs.items():
        graphs_in_cif, bonds, angles = graphs
        structure = StructureCode.objects.get(refcode=refcode)
        upload_graphs_to_db(graphs_in_cif, structure, bonds, angles)


def main(args, all_data=False, user_refcodes=''):
    """
    user_refcodes: {'path_file': 'user_refcode', ...}
    example: {'C:\dev\cifs\my1.cif': 'SDFIREJS'}
    """
    if not user_refcodes:
        user_refcodes = dict()
    if 'all_data' in args:
        all_data = True
    logger_main.info(f"Counting cif files in a specified directory")
    cif_files = get_files(args)
    # split an array of cif files in parts of CHUNK_SIZE size
    for i in range(0, len(cif_files), CHUNK_SIZE):
        logger_main.info(f"Start reading cif files")
        cif_blocks = manager_collect_cifs(cif_files[i:i + CHUNK_SIZE], user_refcodes)
        # Add all data from cif file
        if all_data:
            logger_main.info(f"Start adding all information from the cif to the database")
            add_all_cif_data(cif_blocks)
        # Adding coordinates and parameters with deviations
        logger_main.info(f"Start adding coordinates and cell parameters with deviations")
        manager_add_coords_and_params_to_db(cif_blocks)
        # Create a queue for multi-threaded processing and get a list of structures that should be added
        logger_main.info(f"Start creating a queue for multi-threaded processing of cif files")
        queue, refcodes_to_graph = create_queue(cif_blocks)
        # Creating molecule graphs
        logger_main.info(f"Start generating molecule graphs in multi-threaded mode")
        graphs: dict = create_graph_c(queue)
        # Checking which structures were not processed
        not_added_structures = set(refcodes_to_graph).difference(set(graphs.keys()))
        logger_main.info(f"Graph Generation Results:\n"
                         f"\tTotal structures: {len(cif_blocks)}\n"
                         f"\tStructures without coordinates: {len(cif_blocks) - len(refcodes_to_graph)}\n"
                         f"\tAdded {len(graphs.keys())} structures of {len(refcodes_to_graph)}\n"
                         f"\tNot added {len(not_added_structures)} structures (addition error in {round(len(not_added_structures) / len(refcodes_to_graph) * 100, 2)} % cases)\n"
                         f"\tList of unadded structures:\n"
                         f"\t\t{', '.join(not_added_structures)}")
        # Adding graphs to the database
        logger_main.info(f"Start adding graphs into the database")
        manager_upload_graphs_to_db(graphs)
        # Adding substructure info
        logger_main.info(f"Start adding substructure information")
        add_substructure_filters(graphs.keys(), NUM_OF_PROC)
    logger_main.info(f"Script was finished successfully!")
    return 0


class Command(BaseCommand):
    help = 'Add new data to database from cif files.'

    def handle(self, *args, all_data=False, user_refcodes='', **options):
        main(args, all_data=False, user_refcodes='')

    def add_arguments(self, parser):
        parser.add_argument(
            nargs='+',
            type=str,
            help='Path to cif file(s) or directory path with cif files',
            dest='args'
        )

