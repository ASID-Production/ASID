import requests
import json
import sys
import argparse


def structureSearch(req, body, token=None):
    body = json.loads(body)

    if token is not None and token != 'None':
        headers = {'Authorization': f'Token {token}'}
        reqf = lambda: requests.get(req, data=body, headers=headers)
    else:
        reqf = lambda: requests.get(req, data=body)
    file.write('\n'+str(req)+'\n'+str(body)+'\n'+str(headers)+'\n')
    data = reqf()
    data = data.content.decode(data.apparent_encoding)
    file.write(data + '\n')
    sys.stdout.write(data+'\n')
    data = json.loads(data)
    for iter in range(1, data['max_iter_num']+1):
        body['iter_num'] = iter
        data = reqf()
        data = data.content.decode(data.apparent_encoding)
        file.write(data + '\n')
        sys.stdout.write(data+'\n')
        # data = json.loads(data)


def search(req, token=None):
    if token and token != 'None':
        headers = {'Authorization': f'Token {token}'}
        data = requests.get(req, headers=headers)
    else:
        data = requests.get(req)
    data = data.content.decode(data.apparent_encoding)
    file.write(data + '\n')
    sys.stdout.write(data+'\n')
    #data = json.loads(data)
    while data['next'] is not None:
        data = requests.get(data['next'])
        data = data.content.decode(data.apparent_encoding)
        file.write(data + '\n')
        sys.stdout.write(data+'\n')
        #data = json.loads(data)


FUNCTIONS = {'structureSearch': (structureSearch, ('request', 'body', 'token')),
             'search': (search, ('request', 'token'))
             }

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('-r', '--request')
parser.add_argument('-b', '--body')
parser.add_argument('-f', '--function')
parser.add_argument('-t', '--token')


if __name__ == '__main__':
    file = open('debug.txt', 'w')
    file.write('2\n')
    for arg in sys.argv:
        file.write(str(arg)+'\n')
    file.write('2 done\n')
    file.write(str(sys.argv) + '\n')
    args = parser.parse_args(sys.argv)
    if args.body == 'True':
        args.body = sys.stdin.readline()
        file.write(args.body+'\n')
    func = FUNCTIONS.get(args.function, None)
    if func:
        f = func[0]
        req_args = func[1]
        p_args = []
        for req_arg in req_args:
            p_args.append(args.__getattribute__(req_arg))
        file.write(str(p_args) + '\n')
        f(*p_args)
        file.close()
