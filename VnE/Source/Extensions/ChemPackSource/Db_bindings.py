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


import requests
import json

search_types = ['refcode', 'name', 'elements', 'doi', 'authors', 'cell']

SETUP = False
SESSION = None
SERVER_PROC = None


class Session:
    def __init__(self):
        self.url_base = None
        self.user_token = None

    def connectServer(self, address, port):
        url_base = f'http://{address}:{port}'
        if url_base != self.url_base:
            self.user_token = None
            self.url_base = url_base

    def login(self, name, passwd):
        if self.url_base is not None:
            b = requests.post(f'{self.url_base}/api/auth/users/', data={
                "email": "user@example.com",
                "username": f"{name}",
                "first_name": "none",
                "last_name": "none",
                "password": f"{passwd}"
            })
            token = requests.post(f'{self.url_base}/api/auth/token/login/', data={
                "username": f"{name}",
                "password": f"{passwd}"
            })
            self.user_token = json.loads(token.text)['auth_token']

    def startServer(self, port):
        global SERVER_PROC
        if SERVER_PROC is None:
            import subprocess
            import os
            if os.name == 'nt':
                proc_cmd1 = '../api_database/venv/Scripts/python.exe ../api_database/django_project/manage.py migrate'.split(' ')
                proc_cmd2 = '../api_database/venv/Scripts/python.exe ../api_database/django_project/manage.py runserver'.split(' ')
            elif os.name == 'posix':
                proc_cmd1 = '../api_database/venv/bin/python3 ../api_database/django_project/manage.py migrate'.split(' ')
                proc_cmd2 = '../api_database/venv/bin/python3 ../api_database/django_project/manage.py runserver'.split(' ')
            proc = subprocess.Popen(proc_cmd1)
            while proc.poll() is None:
                pass
            proc = subprocess.Popen(proc_cmd2)
            SERVER_PROC = proc
            try:
                proc.wait(timeout=5)
                raise Exception('Server Error')
            except subprocess.TimeoutExpired:
                pass
            self.url_base = f'http://localhost:{port}'
            self.user_token = None



def search(text, search_type):
    if search_type == 'substructure':
        data_out = structureSearch(text)
        return data_out
    data = requests.get(f'{SESSION.url_base}/api/v1/structures/?{search_type}={text}&limit=1000')
    data = data.content.decode(data.apparent_encoding)
    data = json.loads(data)
    data_out = data['results']
    while data['next'] is not None:
        data = requests.get(data['next'])
        data = data.content.decode(data.apparent_encoding)
        data = json.loads(data)
        data_out += data['results']
    return data_out


SESSION: Session


def get_full_info(id):
    data = requests.get(f'{SESSION.url_base}/api/v1/structures/{id}')
    data = data.content.decode(data.apparent_encoding)
    data_out = json.loads(data)
    return data_out


def getCif(id):
    data = requests.get(f'{SESSION.url_base}/api/v1/structures/{id}/download/')
    return data.content


def uploadFile(file):
    import os

    global TOKEN
    if SESSION.user_token is not None:
        url = f'{SESSION.url_base}/api/v1/structures/upload/'
        files = [('file', (os.path.basename(file), open(file, 'rb'), 'application/octet-stream'))]
        headers = {'Authorization': f'Token {SESSION.user_token}'}
        resp = requests.request('POST', url, headers=headers, data={}, files=files)


def structureSearch(struct):
    from .DrawerWidget import Drawing
    struct: Drawing

    def getEdges(struct):
        edges = []
        mem = []
        for point in struct.connections:
            if point.atom_type == 'H':
                continue
            mem.append(point)
            for point2 in struct.connections[point]:
                if point2.atom_type == 'H':
                    continue
                if point2 not in mem:
                    edges += [f"[{point.seq}, {point2.seq}, {{}}]"]
                else:
                    continue
        return edges
    body = {'search_type': 'substructure', 'nodes': None, 'edges': None}
    h_num = lambda atom: len([x for x in struct.connections[atom] if x.atom_type == 'H'])
    nodes = [f"[{point.seq}, {{'type': '{point.atom_type}', 'Hnum': {h_num(point)}}}]" for point in struct.connections if point.atom_type != 'H']
    edges = getEdges(struct)
    body['nodes'] = nodes
    body['edges'] = edges
    data = requests.get(f'{SESSION.url_base}/api/v1/structures/search/?limit=1000', data=body)
    data = data.content.decode(data.apparent_encoding)
    data = json.loads(data)
    data_out = data['results']
    while data['next'] is not None:
        data = requests.get(data['next'], data=body)
        data = data.content.decode(data.apparent_encoding)
        data = json.loads(data)
        data_out += data['results']
    return data_out


def create_table(data):
    text_dict = {}
    text_dict['RefCode'] = data['refcode']
    text_dict['Author(s)'] = ', '.join([f'{x["family_name"]} {x["initials"]}' for x in data['authors']])
    text_dict['Reference'] = f'{data["publication"]["journal"]} ({data["publication"]["year"]}), {data["publication"]["volume"]}, {data["publication"]["page"]}, DOI: {data["publication"]["doi"]}'
    text_dict['CCDC number'] = data['CCDC_number']

    text_dict['Formula'] = data['formula']['formula_moiety']
    text_dict['Compound'] = data['compound_name']['systematic_name']
    text_dict['Spacegroup'] = f'Name: {data["cell"]["spacegroup"]["name"]:<10}Number: {data["cell"]["spacegroup"]["number"]}'
    text_dict['Cell'] = f'a:          {data["cell"]["a"]:<10.3f}b:      {data["cell"]["b"]:<10.3f}c:            {data["cell"]["c"]:<10.3f}\n' \
                        f'alpha:   {data["cell"]["al"]:<10.2f}beta:  {data["cell"]["be"]:<10.2f}gamma: {data["cell"]["ga"]:<10.2f}\n' \
                        f'Volume:  {data["reduced_cells"][0]["volume"]}'

    return text_dict


def setup():
    global SESSION
    SESSION = Session()
    return


if SETUP is False:
    setup()


if __name__ == '__main__':
    get_full_info(search('ZUDGAJ01', 'refcode')[0]['id'])
    'ABOPUG'
