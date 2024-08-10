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


from copy import copy
import numpy as np
from .MoleculeClass import Atom, Bond, Molecule, MoleculeSystem

import debug


class Pack:
    def __init__(self):
        self.nodes = []
        self.mapped = 0
        self.solutions = []

    def find_node(self, node):
        if isinstance(node, Node):
            try:
                return self.nodes[self.nodes.index(node)]
            except ValueError:
                return None
        elif isinstance(node, Atom):
            atom = node
            for node in self.nodes:
                if atom is node.struct_atom or atom is node.assigned_node:
                    return node
        return None

    def addNode(self, node):
        if node not in self.nodes:
            self.nodes.append(node)
            node.pack = self
            if node.assigned_node is None:
                self.mapped += 1

    def checkMapping(self):
        if self.mapped == 0:
            return True
        else:
            return False

    def addSolution(self, pack):
        ch = [x for x in range(len(self.solutions))]
        if pack in self.solutions:
            return
        for sol_i, psol in enumerate(self.solutions):
            for i, node in enumerate(self.nodes):
                node_s = pack.nodes[i]
                if psol.nodes[i].assigned_node == node_s.assigned_node:
                    pass
                else:
                    ch.remove(sol_i)
                    break
            if sol_i in ch:
                break
        if not ch:
            self.solutions.append(pack)
            for i, node in enumerate(self.nodes):
                node_s = pack.nodes[i]
                node.addRealAtom(node_s.assigned_node.struct_atom)

    @staticmethod
    def concatenatePacks(*packs):
        ret_pack = Pack()
        for pack in packs:
            for node in pack.nodes:
                ret_pack.addNode(node)
        return ret_pack

    def splitPack(self):
        def recSplit(node, nodes, mem=None, pack=None):
            if pack is None:
                pack = Pack()
            if mem is None:
                mem = {}
            if node in mem:
                return pack
            else:
                mem[node] = True
                nodes.remove(node)
                pack.addNode(node)
                for node_c in node.connect:
                    pack = recSplit(node_c, nodes, mem, pack)
                return pack

        packs = []
        nodes = copy(self.nodes)
        while nodes:
            packs.append(recSplit(nodes[0], nodes))
        return packs

    def __deepcopy__(self, memodict=None):
        memodict_l = {}
        new_pack = type(self)()
        for node in self.nodes:
            new_node = Node(new_pack)
            memodict_l[node] = new_node
            new_node.assignNode(node.assigned_node)
            new_node.setStructAtom(node.struct_atom)
            new_pack.addNode(new_node)
        for node in memodict_l:
            new_node = memodict_l[node]
            for connect_node in node.connect:
                new_node.addConnect(memodict_l[connect_node])
        return new_pack


class Node:

    def __init__(self, pack, connect=None, str_atom=None, assigned_node=None):
        self.pack = pack
        if connect is None:
            self.connect = []
        else:
            self.connect = copy(connect)
        self.struct_atom = str_atom
        self.assigned_node = assigned_node
        self.real_atoms = []
        if self.struct_atom:
            if type(self.struct_atom.atom_type) is list:
                self.atom_type = self.struct_atom.atom_type
            else:
                self.atom_type = [self.struct_atom.atom_type]
        elif self.assigned_node:
            if type(self.assigned_node.atom_type) is list:
                self.atom_type = self.assigned_node.atom_type
            else:
                self.atom_type = [self.assigned_node.atom_type]
        else:
            self.atom_type = None

        self.pack.addNode(self)

    def addConnect(self, node):
        if node not in self.connect:
            self.connect.append(node)
            node.addConnect(self)

    def removeConnect(self, node):
        try:
            self.connect.remove(node)
            node.removeConnect(self)
        except ValueError:
            pass

    def delete(self):
        for node in self.connect:
            node.removeConnect(self)
        self.pack.nodes.remove(self)
        if self.assigned_node is None:
            self.pack.mapped -= 1

    def setStructAtom(self, struct_atom: Atom):
        if struct_atom is None:
            if self.assigned_node is None:
                self.atom_type = None
            self.struct_atom = None
            return
        if self.assigned_node:
            if type(struct_atom.atom_type) is list:
                z = True
                for atom_type in struct_atom.atom_type:
                    if atom_type in self.atom_type:
                        z = False
                        break
                if z:
                    raise ValueError('Incompatible atom type')
            elif struct_atom.atom_type not in self.atom_type:
                raise ValueError('Incompatible atom type')
        self.struct_atom = struct_atom
        if struct_atom.atom_type != self.atom_type:
            if type(struct_atom.atom_type) is list:
                self.atom_type = self.struct_atom.atom_type
            else:
                self.atom_type = [self.struct_atom.atom_type]

    def assignNode(self, node):
        if node is None:
            if self.struct_atom is None:
                self.atom_type = None
            self.assigned_node = None
            if self.assigned_node is not None:
                self.pack.mapped += 1
            return
        if self.struct_atom:
            z = True
            for type in node.atom_type:
                if type in self.atom_type:
                    z = False
                    break
            if z:
                raise ValueError('Incompatible atom type')
        elif self.atom_type is None:
            self.atom_type = node.atom_type
        self.assigned_node = node
        self.pack.mapped -= 1

    def addRealAtom(self, atom: Atom):
        if atom.atom_type in self.atom_type:
            self.real_atoms.append(atom)
        else:
            raise ValueError('Incompatible atom type')

    def deleteRealAtom(self, atom: Atom):
        if atom in self.real_atoms:
            self.real_atoms.remove(atom)

    @classmethod
    def typeSubseteq(cls, lnode, rnode):
        for atom_type in lnode.atom_type:
            if atom_type not in rnode.atom_type:
                return False
        return True

    def __deepcopy__(self, memodict=None):
        pack = self.pack.__deepcopy__()
        i = self.pack.nodes.index(self)
        return pack.nodes[i]

    def __str__(self):
        return f'{self.pack.nodes.index(self)}:{self.atom_type}'


class Condition:

    def __init__(self, comp_func, value_func, nodes):
        self.nodes = nodes
        self.packs = {}
        for node in self.nodes:
            pack = node.pack
            if pack in self.packs:
                self.packs[pack].append(node)
            else:
                self.packs[pack] = [node]
        self.packs_list = list(self.packs.keys())
        self.comp_func = comp_func
        self.value_func = value_func
        self.last_eval = []

    def eval(self):

        self.packs = {}
        for node in self.nodes:
            pack = node.pack
            if pack in self.packs:
                self.packs[pack].append(node)
            else:
                self.packs[pack] = [node]
            self.packs_list = list(self.packs.keys())

        def rec_eval(packs, ind=None, result=None):
            if result is None:
                result = []
            if ind is None:
                ind = []
            if packs:
                packs_c = packs.copy()
                ind_c = ind.copy()
                pack = packs_c.pop()
                ind_c.append(-1)
                for i in range(len(pack.solutions)):
                    ind_c[-1] = i
                    result = rec_eval(packs_c, ind_c, result)
                return result
            else:
                ind = ind[::-1]
                packs = list(self.packs_list)
                atoms = []
                for node in self.nodes:
                    atoms.append(node.real_atoms[ind[packs.index(node.pack)]])
                value = self.value_func(*atoms)
                if self.comp_func(value):
                    result.append(ind.copy())
                return result
        for node in self.nodes:
            if node.real_atoms:
                continue
            else:
                return []

        self.last_eval = rec_eval(self.packs_list)
        return self.last_eval

    def setCompFunc(self, func):
        self.comp_func = func

    def setValueFunc(self, func):
        self.value_func = func

    @staticmethod
    def intersection(*conditions, eval=False):
        if eval or not all([True if x.last_eval else False for x in conditions]):
            for condition in conditions:
                condition.eval()
        solutions = [{pack: [solution[i] for solution in condition.last_eval] for i, pack in enumerate(condition.packs)} for condition in conditions]
        packs = {}
        ret = solutions[0]
        for i, cond1 in enumerate(solutions[:-1]):
            for z, cond2 in enumerate(solutions[i+1:]):
                packs = list(set(cond1.keys()) & set(cond2.keys()))
                if packs:
                    sol1 = [[cond1[pack][x] for pack in packs] for x in range(len(cond1[packs[0]]))]
                    sol2 = [[cond2[pack][x] for pack in packs] for x in range(len(cond2[packs[0]]))]
                    soln1 = [True if x in sol2 else False for x in sol1]
                    soln2 = [True if x in sol1 else False for x in sol2]
                    cond1 = {key: [cond1[key][i] for i, x in enumerate(soln1) if x] for key in cond1.keys()}
                    cond2 = {key: [cond2[key][i] for i, x in enumerate(soln2) if x] for key in cond2.keys()}
                    a = 0
                    solutions[i] = cond1
                    solutions[z+i+1] = cond2
                    a = 0
        return solutions


def avgDiff(atoms, split=1):
    coord1 = np.sum([x.coord for x in atoms[:split]], axis=0)/split
    coord2 = np.sum([x.coord for x in atoms[split:]], axis=0)/(len(atoms)-split)
    return np.linalg.norm(coord2-coord1)


def dist(atom1: Atom, atom2: Atom):
    return np.linalg.norm(atom1.coord-atom2.coord)


def plane_eq(atom1: Atom, atom2: Atom, atom3: Atom):
    relV1 = atom1.coord - atom2.coord
    relV2 = atom3.coord - atom2.coord
    normal = np.cross(relV1, relV2)
    if np.count_nonzero(normal) == 0:
        raise ValueError("Vectors are co-incident")
    d = -np.dot(normal, atom2.coord)
    return np.array([normal[0], normal[1], normal[2], d])


def meanPlaneEq(*atoms):
    points = np.array([atom.coord for atom in atoms])
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid
    _, _, V = np.linalg.svd(centered_points)
    normal_vector = V[-1]
    d = -np.dot(normal_vector, centroid)
    plane = np.append(normal_vector, d)

    return plane


def pointPlaneDist(plane, atom):
    a, b, c, d = plane
    distance = np.abs(a * atom.coord[0] + b * atom.coord[1] + c * atom.coord[2] + d) / np.sqrt(a ** 2 + b ** 2 + c ** 2)
    return distance


def maxMeanPlaneDiff(*atoms):
    mean_plane = meanPlaneEq(*atoms)

    max_distance = -1

    for atom in atoms:
        distance = pointPlaneDist(mean_plane, atom)

        if distance > max_distance:
            max_distance = distance

    return max_distance


def angle(atoms, rad=True, split=-1):
    """
    Angle, torsion, improper, plane angle calculation
    :param atoms: 3-6 atoms
    :type atoms: Atom
    :param rad: True - return in radian, False - in degree
    :type rad: bool
    :return: Angle, float
    """
    if split == -1:
        if len(atoms) == 3:
            v1 = atoms[1].coord - atoms[0].coord
            v2 = atoms[1].coord - atoms[2].coord
            cosC = (v1 @ v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
            angleC = np.arccos(cosC)
            if rad:
                return angleC
            else:
                return (angleC / np.pi) * 180
        else:
            plane1 = plane_eq(*atoms[:3])
            plane2 = plane_eq(*atoms[-3:])
            np.sqrt(sum([x**2 for x in plane1[:3]])) * np.sqrt(sum([x**2 for x in plane2[:3]]))
            cos = (sum([x ** 2 for x in plane1[:3]]) + sum([x**2 for x in plane2[:3]]) - sum([x**2 for x in plane1[:3] - plane2[:3]])) /\
                  (2 * np.sqrt(sum([x ** 2 for x in plane1[:3]])) * np.sqrt(sum([x**2 for x in plane2[:3]])))
            angle = np.arccos(cos)
            if rad:
                return angle
            else:
                return (angle / np.pi) * 180
    if split != -1:
        plane_eq1 = meanPlaneEq(*atoms[:split])
        plane_eq2 = meanPlaneEq(*atoms[split:])
        cosC = (plane_eq1[:3] @ plane_eq2[:3]) / (np.linalg.norm(plane_eq1[:3]) * np.linalg.norm(plane_eq2[:3]))
        angleC = np.arccos(cosC)
        if rad:
            return angleC
        else:
            return (angleC / np.pi) * 180


def findSubGraph(pack: Pack, sub_pack: Pack):

    def check(sub_node, node, comp_mat):
        sub_nodes_ind = [sub_node.pack.nodes.index(x) for x in sub_node.connect]
        nodes_ind = [pack.nodes.index(x) for x in node.connect]
        mat = np.array([[comp_mat[sn][n] for n in nodes_ind] for sn in sub_nodes_ind])
        mapped = np.zeros(mat.shape[0], dtype=int)

        for i, csn in enumerate(sub_node.connect):
            if csn.assigned_node is not None:
                try:
                    mat[i, :] = False
                    mat[:, node.connect.index(csn.assigned_node)] = False
                    mat[i, node.connect.index(csn.assigned_node)] = True
                    mapped[i] = 1
                except ValueError:
                    return False
        while max([sum(x) for x in mat]) != 1 and max([sum(z) for z in mat.T]) != 1:
            x = [sum(x) for x in mat]
            if 0 in x:
                return False
            x_i = min([(x[z],z) for z in range(len(x)) if x[z]!=0 and (not mapped[z])], key=lambda v: v[0])[1]
            y = [(sum(mat.T[z]), z) for z in range(len(mat.T)) if mat[x_i,z]]
            y_i = min(y, key=lambda v: v[0])[1]
            mat[x_i,:] = False
            mat[:,y_i] = False
            mat[x_i,y_i] = True
            mapped[x_i] = 1
        if 0 in [sum(x) for x in mat]:
            return False
        else:
            return True

    def check_prob(node: Node, sub_node: Node):

        def rec_find(mat, mapped, pos_pos):
            x = [sum(x) for x in mat]
            x_i = [a for a in range(len(x)) if x[a]!=0 and not mapped[a]]
            if 0 in x:
                return pos_pos
            if not all(mapped):
                for x_curr in x_i:
                    y = [(sum(z), i) for i,z in enumerate(mat.T) if mat[x_curr, i]]
                    def asd(x):
                        a = y.index(x)
                        y[a] = -1
                        return a
                    y_i = [z[1] for z in y]
                    for y_curr in y_i:
                        mat_copy = np.copy(mat)
                        mapped_copy = np.copy(mapped)
                        mat_copy[x_curr, :] = False
                        mat_copy[:, y_curr] = False
                        mat_copy[x_curr, y_curr] = True
                        mapped_copy[x_curr] = 1
                        pos_pos = rec_find(mat_copy, mapped_copy, pos_pos)
                return pos_pos
            else:
                if not pos_pos or not np.any([(mat == x).all() for x in pos_pos]):
                    pos_pos.append(mat)
                return pos_pos

        sub_nodes_ind = [sub_node.pack.nodes.index(x) for x in sub_node.connect]
        nodes_ind = [pack.nodes.index(x) for x in node.connect]
        if len(sub_nodes_ind) > len(nodes_ind):
            return []
        mat = np.array([[comp_mat[sn][n] for n in nodes_ind] for sn in sub_nodes_ind])
        mapped = np.zeros(mat.shape[0], dtype=int)
        for i, csn in enumerate(sub_node.connect):
            csn: Node
            if csn.assigned_node is not None:
                try:
                    mat[i, :] = False
                    mat[:, node.connect.index(csn.assigned_node)] = False
                    mat[i, node.connect.index(csn.assigned_node)] = True
                except ValueError:
                    return []

        while (max([sum(x) for x in mat]) != 1 and max([sum(z) for z in mat.T]) != 1) and (1 in [sum(x) for i, x in enumerate(mat) if mapped[i] == 0]):
            x = [sum(x) for x in mat]
            x_i = min([(x[z],z) for z in range(len(x)) if x[z]!=0 and (not mapped[z])], key=lambda v: v[0])[1]
            y = [(sum(mat.T[z]), z) for z in range(len(mat.T)) if mat[x_i,z]]
            y_i = min(y, key=lambda v: v[0])[1]
            mat[x_i, :] = False
            mat[:, y_i] = False
            mat[x_i, y_i] = True
            mapped[x_i] = 1

        pos_pos = rec_find(mat, mapped, [])
        return pos_pos

    def recurse(node: Node, sub_node: Node, memory=None, result=None):
        sub_node_copy_init: Node = sub_node.__deepcopy__()
        if result is None:
            result = []
        if memory is None:
            memory = {'Pack_p': [node], 'SubPack_p': [sub_node]}
        else:
            memory = {'Pack_p': copy(memory['Pack_p']), 'SubPack_p': copy(memory['SubPack_p'])}
        if sub_node_copy_init.assigned_node is not None:
            return True, sub_node_copy_init, result
        else:
            sub_node_copy_init.assignNode(node)
            memory['Pack_p'].append(node)
            memory['SubPack_p'].append(sub_node_copy_init.struct_atom)
        if sub_node_copy_init.pack.checkMapping():
            result.append(sub_node_copy_init)
            return True, sub_node_copy_init, result
        probs = check_prob(node, sub_node_copy_init)
        sub_node_copy: Node
        if probs == []:
            sub_node_copy = sub_node_copy_init

        for prob in probs:
            sub_node_copy = sub_node_copy_init.__deepcopy__()
            for i in range(len(prob)):
                b = np.where(prob[i] == True)[0][0]
                prob[i][b] = check(sub_node_copy.connect[i], node.connect[b], comp_mat)
            x = [sum(x) for x in prob]
            if 0 in x:
                continue
            else:
                for i in range(len(prob)):
                    if sub_node.connect[i].assigned_node:
                        continue
                    else:
                        b = np.where(prob[i] == True)[0][0]
                        z, node_r, mem_r = recurse(node.connect[b], sub_node_copy.connect[i], memory=memory, result=result)
                        if z:
                            for n in node_r.pack.nodes:
                                n: Node
                                if n.assigned_node is sub_node_copy.assigned_node:
                                    sub_node_copy = n
                                    break
                            if sub_node_copy.pack.checkMapping():
                                result.append(sub_node_copy)
                        else:
                            sub_node_copy = sub_node_copy_init.__deepcopy__()
                            break
        if all([x.assigned_node for x in sub_node_copy.connect]):
            return True, sub_node_copy, result
        else:
            return False, sub_node_copy, result

    comp_mat = [[Node.typeSubseteq(node, sub_node) for node in pack.nodes] for sub_node in sub_pack.nodes]
    comp_mat = np.array(comp_mat)
    ret = []

    if len(sub_pack.nodes) > len(pack.nodes):
        return []
    elif len(sub_pack.nodes) == 0:
        return [pack]

    for i, first_node in enumerate(pack.nodes):
        if comp_mat[0][i]:
            ret.append(recurse(pack.nodes[i], sub_pack.nodes[0]))

    return ret


def save_graph(system: MoleculeSystem):
    import json

    mol: Molecule = system.children[0]
    atoms = []
    for atom in mol.children:
        atom: Atom
        point = atom.point()
        if point.pick == 1.0:
            atoms.append(atom)
    struct = {i: {'type': 0, 'bond': []} for i in range(len(atoms))}
    for i, atom in enumerate(atoms[:-1]):
        struct[i]['type'] = atom.atom_type
        for bond in atom.bonds():
            if bond.parents()[bond.parents().index(atom)-1] in atoms[i+1:]:
                struct[atoms.index(bond.parents()[bond.parents().index(atom) - 1])]['type'] = bond.parents()[bond.parents().index(atom)-1].atom_type
                struct[atoms.index(bond.parents()[bond.parents().index(atom) - 1])]['bond'].append(i)
                struct[i]['bond'].append(atoms.index(bond.parents()[bond.parents().index(atom) - 1]))
    file_json = json.dumps(struct, separators=(',', ':'))
    path = '\\'.join(__file__.split('\\')[:-1])
    out = open(f'{path}\\test.json', 'w')
    out.write(file_json)
    out.close()


def test(system: MoleculeSystem, file_name='test'):
    import json

    path = '\\'.join(__file__.split('\\')[:-1])
    file = json.load(open(f'{path}\\{file_name}.json', 'r'))
    sub_pack = Pack()
    nodes = [Node(sub_pack, str_atom=Atom(np.array([0,0,0]), atom_type=file[i]['type'])) for i in file]
    for i, node in enumerate(nodes):
        for cnode_i in file[str(i)]['bond']:
            cnode = nodes[cnode_i]
            node.addConnect(cnode)

    pack = Pack()

    mol: Molecule = system.children[0]
    atoms = []
    nodes = []
    nodes_atom = {}
    for atom in mol.children:
        atom: Atom
        atoms.append(atom)
        nodes.append(Node(pack, str_atom=Atom(np.array([0, 0, 0]), atom_type=atom.atom_type)))
        nodes_atom[nodes[-1]] = atom
    struct = {i: {'type': 0, 'bond': []} for i in range(len(atoms))}
    for i, atom in enumerate(atoms[:-1]):
        struct[i]['type'] = atom.atom_type
        for bond in atom.bonds():
            if bond.parents()[bond.parents().index(atom) - 1] in atoms[i + 1:]:
                struct[atoms.index(bond.parents()[bond.parents().index(atom) - 1])]['type'] = bond.parents()[
                    bond.parents().index(atom) - 1].atom_type
                struct[atoms.index(bond.parents()[bond.parents().index(atom) - 1])]['bond'].append(i)
                struct[i]['bond'].append(atoms.index(bond.parents()[bond.parents().index(atom) - 1]))
                nodes[i].addConnect(nodes[atoms.index(bond.parents()[bond.parents().index(atom) - 1])])
    ret = findSubGraph(pack, sub_pack)
    mols = []
    for res in ret:
        if res[2]:
            sp = None
            for gr in res[2]:
                if gr.pack is not sp:
                    mol = []
                    sp = gr.pack
                    for node in sp.nodes:
                        atom = nodes_atom[node.assigned_node]
                        atom.point().pick = 1.0
                        mol.append(atom)
                    mols.append(mol)
    a = 0