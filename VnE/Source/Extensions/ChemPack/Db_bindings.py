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

import sys
import requests
import json
from PySide6.QtWidgets import QDialog, QLabel, QVBoxLayout, QTextEdit
from PySide6.QtCore import QProcess
import os.path as opath
import os
import base64
import atexit
import subprocess

import debug

search_types = ['refcode', 'name', 'elements', 'doi', 'authors', 'cell', 'formula']

SETUP = False
SESSION = None
SERVER_PROC = None
TEST = None
requests.packages.urllib3.util.connection.HAS_IPV6 = False


class ErrorDialog(QDialog):

    def __init__(self, error, title='Error'):
        QDialog.__init__(self)
        self.setWindowTitle(title)
        text = [f'{x}: {str(y)}' for x, y in error.items()]
        text = '\n'.join(text)
        self.label = QTextEdit(text)
        self.vblayout = QVBoxLayout()
        self.vblayout.addWidget(self.label)
        self.setLayout(self.vblayout)
        self.finished.connect(lambda x: SESSION.deleteErrorDialog(self))


class Session:
    def __init__(self):
        self.url_base = None
        self.user_token = None
        self.usr_login = None
        self.usr_pass = None
        self.ready = False
        self.last_connect_start_op = lambda: True
        self.error_dialog = []

    def connectServer(self, address, port):
        self.last_connect_start_op = lambda: self.connectServer(address, port)
        url_base = f'http://{address}:{port}'
        if url_base != self.url_base:
            self.user_token = None
            self.url_base = url_base
            self.ready = True

    def login(self, name, passwd):
        if not self.ready:
            return
        if name and passwd:
            if self.url_base is not None:
                self.usr_login, self.usr_pass = name, passwd
                b = requests.post(f'{self.url_base}/api/auth/users/', data={
                    "email": f"{self.usr_login}@example.com",
                    "username": f"{name}",
                    "first_name": "none",
                    "last_name": "none",
                    "password": f"{passwd}"
                })
                sh = b.ok
                if not b.ok:
                    error = json.loads(b.text)
                    try:
                        err = error['username']
                        if 'A user with that username already exists.' in err:
                            sh = True
                    except KeyError:
                        pass
                if not sh:
                    SESSION.error_dialog.append(ErrorDialog(error))
                    SESSION.error_dialog[-1].show()
                token = requests.post(f'{self.url_base}/api/auth/token/login/', data={
                    "username": f"{name}",
                    "password": f"{passwd}"
                })
                if not token.ok:
                    error = json.loads(token.text)
                    SESSION.error_dialog.append(ErrorDialog(error))
                    SESSION.error_dialog[-1].show()
                    return
                self.user_token = json.loads(token.text)['auth_token']

    def startServer(self, port):
        global SERVER_PROC
        self.last_connect_start_op = lambda: self.startServer(port)
        if SERVER_PROC is None:
            root = opath.normpath(f'{opath.dirname(__file__)}/../../../..')
            path = opath.normpath(f'{root}/api_database/django_project/manage.py')
            if os.name == 'nt':
                proc_cmd1 = (f'{root}\\venv\\Scripts\\python.exe', f'{path}', 'migrate')
                proc_cmd2 = (f'{root}\\venv\\Scripts\\python.exe', f'{path}', 'runserver')
            elif os.name == 'posix':
                proc_cmd1 = (f'{root}/venv/bin/python', f'{path}', 'migrate')
                proc_cmd2 = (f'{root}/venv/bin/python', f'{path}', 'runserver')
            proc = subprocess.Popen(proc_cmd1)
            while proc.poll() is None:
                pass
            proc.kill()
            proc = subprocess.Popen(proc_cmd2)
            SERVER_PROC = proc
            try:
                proc.wait(timeout=5)
                raise Exception('Server Error')
            except subprocess.TimeoutExpired:
                pass
            self.url_base = f'http://localhost:{port}'
            self.user_token = None
            self.ready = True
        else:
            SERVER_PROC.kill()
            SERVER_PROC = None
            self.startServer(port)
            self.login(self.usr_login, self.usr_pass)

    def triggerLastConnectOp(self):
        self.last_connect_start_op()

    def deleteErrorDialog(self, dialog):
        self.error_dialog.remove(dialog)

SESSION: Session


def search(text, search_type, db_type='cryst', exact=None, process=None, CSD=True, COD=True, ICSD=True, user_db=True):
    url_mods = {'cryst': 'api/v1/structures',
                'qm': 'api/v1/qc_structures'}
    dbs = {'CSD': CSD,
           'COD': COD,
           'ICSD': ICSD,
           'user_db': user_db}
    db_string = '&'.join([f'{x}=False' for x in dbs if not dbs[x]])
    if process is None:
        process = QProcess()
    if not SESSION.ready:
        process.kill()
        return
    url_mod = url_mods.get(db_type, 'api/v1/structures')
    if search_type == 'substructure':
        process = structureSearch(text, url_mod, process, db_string=db_string)
        return process
    token = SESSION.user_token
    if exact is None:
        req = f'{SESSION.url_base}/{url_mod}/?{search_type}={text}&limit=10000000'
    else:
        req = f'{SESSION.url_base}/{url_mod}/?{search_type}={text}&exact={exact}&limit=10000000'
    if db_string:
        req = '&'.join([req, db_string])

    root = opath.normpath(f'{opath.dirname(__file__)}/../../../..')
    path = opath.normpath(f'{root}/VnE/Source/Extensions/ChemPack/searchProcess.py')
    if os.name == 'nt':
        prog = opath.normpath(f'{root}\\venv\\Scripts\\python.exe')
    elif os.name == 'posix':
        prog = opath.normpath(f'{root}/venv/bin/python3')
    process.startCommand(f'{prog} {path} -f search -r "{req}" -t {token}')
    process.waitForStarted()
    return process


def get_full_info(id, db_type='cryst'):
    url_mods = {'cryst': 'api/v1/structures',
                'qm': 'api/v1/qc_structures'}
    url_mod = url_mods.get(db_type, 'api/v1/structures')
    if SESSION.user_token is not None:
        headers = {'Authorization': f'Token {SESSION.user_token}'}
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}', headers=headers)
    else:
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}')
    if not data.ok:
        error = json.loads(data.text)
        SESSION.error_dialog.append(ErrorDialog(error))
        SESSION.error_dialog[-1].show()
        return {}
    data = data.content.decode(data.apparent_encoding)
    data_out = json.loads(data)
    return data_out


def get_image_temp(id, db_type='cryst', w=250, h=250):
    root = opath.normpath(f'{opath.dirname(__file__)}/../../..')
    path = opath.join(root, 'temp')
    o_file = opath.join(path, f'{id}.svg')
    return get_image(id, o_file, db_type, w, h)


def get_image(id, o_file_path, db_type='cryst', w=250, h=250):
    url_mods = {'cryst': 'api/v1/structures',
                'qm': 'api/v1/qc_structures'}
    url_mod = url_mods.get(db_type, 'api/v1/structures')
    if SESSION.user_token is not None:
        headers = {'Authorization': f'Token {SESSION.user_token}'}
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}/export/2d/?h={h}&w={w}&file=0&f=svg', headers=headers)
    else:
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}/export/2d/?h={h}&w={w}&file=0&f=svg')
    if not data.ok:
        error = json.loads(data.text)
        SESSION.error_dialog.append(ErrorDialog(error))
        SESSION.error_dialog[-1].show()
        return None
    data = data.content.decode(data.apparent_encoding)
    if data == '0':
        return None
    image_str = data
    image_data = image_str
    path = o_file_path
    image = open(path, 'w')
    image.write(image_data)
    image.close()
    return o_file_path


def getImageFromFile(file_path, o_file_path, format='gif', w=250, h=250):
    formats = ['gif', 'cml', 'svg']
    if format not in formats:
        format = formats[0]
    url_mod = 'api/v1/generate'
    data_i = {'file_format ': opath.splitext(file_path)[1][1:],
            'output_format ': format,
            'return_type': 'file',
            'h_size': f'{h}',
            'w_size': f'{w}',
            'name': opath.splitext(opath.basename(o_file_path))[0]}
    files = [('file', (opath.basename(file_path), open(file_path, 'rb'), 'application/octet-stream'))]
    if SESSION.user_token is not None:
        headers = {'Authorization': f'Token {SESSION.user_token}'}
        data = requests.request('GET', f'{SESSION.url_base}/{url_mod}/2d/', headers=headers, data=data_i, files=files)
    else:
        data = requests.request('GET', f'{SESSION.url_base}/{url_mod}/2d/', headers={}, data=data_i, files=files)
    if data.ok:
        image_data = data.content
        image = open('.'.join([opath.splitext(o_file_path)[0], format]), 'wb')
        image.write(image_data)
        image.close()
        return '.'.join([opath.splitext(o_file_path)[0], format])
    else:
        dialog = ErrorDialog(data.text)
        dialog.show()
        return None


def getCif(id, db_type='cryst'):
    url_mods = {'cryst': 'api/v1/structures',
                'qm': 'api/v1/qc_structures'}
    url_mod = url_mods.get(db_type, 'api/v1/structures')
    if db_type == 'qm':
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}/export/cif')
    else:
        data = requests.get(f'{SESSION.url_base}/{url_mod}/{id}/download/')
    if not data.ok:
        error = json.loads(data.text)
        SESSION.error_dialog.append(ErrorDialog(error))
        SESSION.error_dialog[-1].show()
        return b''
    return data.content


def uploadFile(file, ext):
    import os

    def uploadCif(files, headers):
        url = f'{SESSION.url_base}/api/v1/structures/upload/'
        resp = requests.request('POST', url, headers=headers, data={}, files=files)
        return resp

    def uploadXml(files, headers):
        url = f'{SESSION.url_base}/api/v1/qc_structures/upload/vasp/'
        resp = requests.request('POST', url, headers=headers, data={}, files=files)
        return resp

    ext_method = {'.cif': uploadCif,
                  '.xml': uploadXml}

    if SESSION.user_token is not None:
        method = ext_method.get(ext, None)
        if method:
            files = [('file', (os.path.basename(file), open(file, 'rb'), 'application/octet-stream'))]
            headers = {'Authorization': f'Token {SESSION.user_token}'}
            resp = method(files, headers)
            if not resp.ok:
                error = json.loads(resp.text)
                SESSION.error_dialog.append(ErrorDialog(error))
                SESSION.error_dialog[-1].show()
        else:
            SESSION.error_dialog = ErrorDialog('Unsupported file format')
            SESSION.error_dialog[-1].show()


def structureSearch(struct, url_mod, process=None, db_string=''):
    from .DrawerWidget import Drawing
    struct: Drawing
    if process is None:
        process = QProcess()

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
    nodes = [f"[{point.seq}, {{'type': '{' '.join(point.atom_type)}', 'Hnum': {h_num(point)}, 'cord_min': {point.cn[0]}, 'cord_max': {point.cn[1]}}}]" if type(point.atom_type) is list else f"[{point.seq}, {{'type': '{point.atom_type}', 'Hnum': {h_num(point)}, 'cord_min': {point.cn[0]}, 'cord_max': {point.cn[1]}}}]" for point in struct.connections if point.atom_type != 'H']
    edges = getEdges(struct)
    body['nodes'] = nodes
    body['edges'] = edges
    body['chunk_size'] = 10000
    body['iter_num'] = 0
    token = SESSION.user_token
    req = f'{SESSION.url_base}/{url_mod}/search/?limit=10000'
    if db_string:
        req = '&'.join([req, db_string])
    root = opath.normpath(f'{opath.dirname(__file__)}/../../../..')
    path = opath.normpath(f'{root}/VnE/Source/Extensions/ChemPack/searchProcess.py')
    if os.name == 'nt':
        prog = opath.normpath(f'{root}\\venv\\Scripts\\python.exe')
    elif os.name == 'posix':
        prog = opath.normpath(f'{root}/venv/bin/python3')
    body = json.dumps(body)+'\n'
    process.startCommand(f'{prog} {path} -f structureSearch -r "{req}" -t {token} -b True')
    process.waitForStarted()
    process.write(bytes(body, encoding='utf-8'))
    return process


def create_table(data, img=None):
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


def deleteEntry(id, db_type='cryst'):
    url_mods = {'cryst': 'api/v1/structures',
                'qm': 'api/v1/qc_structures'}
    url_mod = url_mods.get(db_type, 'api/v1/structures')
    if SESSION.user_token is not None:
        headers = {'Authorization': f'Token {SESSION.user_token}'}
        data = requests.request('DELETE', f'{SESSION.url_base}/{url_mod}/{id}/', headers=headers)
    else:
        data = requests.request('DELETE', f'{SESSION.url_base}/{url_mod}/{id}')
    if not data.ok:
        error = json.loads(data.text)
        SESSION.error_dialog.append(ErrorDialog(error))
        SESSION.error_dialog[-1].show()


def setup():
    global SESSION
    SESSION = Session()
    return


def term():
    global SERVER_PROC
    if SERVER_PROC:
        if os.name == 'nt':
            from signal import CTRL_BREAK_EVENT
            ret = subprocess.Popen('netstat -ano | findstr 127.0.0.1:8000', shell=True, stdout=subprocess.PIPE).stdout.read()
            ret = ret.decode()
            procs = ret.split('\n')
            pids = [int(x.split(' ')[-1]) for x in procs if x]
            for pid in pids:
                try:
                    os.kill(pid, CTRL_BREAK_EVENT)
                except SystemError as e:
                    pass
            SERVER_PROC = None
            sys.exit()
        if os.name == 'posix':
            from signal import SIGTERM
            from signal import CTRL_BREAK_EVENT
            ret = subprocess.Popen('ps ax | grep runserver', shell=True, stdout=subprocess.PIPE).stdout.read()
            ret = ret.decode()
            procs = ret.split('\n')
            pids = [int(x.split(' ')[0]) for x in procs if x]
            for pid in pids:
                try:
                    os.kill(pid, SIGTERM)
                except SystemError as e:
                    pass
            SERVER_PROC = None
            sys.exit()

if SETUP is False:
    atexit.register(term)
    setup()


if __name__ == '__main__':
    get_full_info(search('ZUDGAJ01', 'refcode')[0]['id'])
    'ABOPUG'
