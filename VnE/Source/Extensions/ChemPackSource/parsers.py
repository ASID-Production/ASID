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


from ... import point_class
import numpy as np
from ..ChemPack import PALETTE, MOLECULE_SYSTEMS
import os
from . import MoleculeClass

import debug

class FileParser:

    def __init__(self):
        self.SUPPORTED_FORMATS = {'.xyz': self.parsXyz,
                                  '.pdb': self.parsPdb,
                                  '.cif': self.parsCif}

    def fracToDec(a, b, c, al, be, ga, coords):
        al = (al / 180) * np.pi
        be = (be / 180) * np.pi
        ga = (ga / 180) * np.pi

        sin = np.sin
        cos = np.cos
        cot = lambda x: np.tan(x) ** -1
        csc = lambda x: np.sin(x) ** -1

        mat = np.array([[a * sin(be) * np.sqrt(1 - (cot(al) * cot(be) - csc(al) * csc(be) * cos(ga)) ** 2), 0, 0],
                        [a * csc(al) * cos(ga) - a * cot(al) * cos(be), b * sin(al), 0],
                        [a * cos(be), b * cos(al), c]])
        mat = mat.transpose()
        for i in range(len(coords)):
            coords[i] = (coords[i] @ mat).astype(dtype=np.float32)
        return coords

    def parsFile(self, file_path, bond=True, root=None):
        _, ext = os.path.splitext(file_path)
        ext = ext.lower()
        if ext not in self.SUPPORTED_FORMATS:
            print('File is not supported')
            return None, None
        return self.SUPPORTED_FORMATS[ext](file_path, bond, root)

    def parsXyz(self, file_path, bond=True, root=None, *args, **kwargs):
        mol_list, atom_list, bonds_l = None, None, None
        file = open(file_path, 'r')
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol_sys.file_name = file_path
        mol = MoleculeClass.Molecule(parent=mol_sys)
        i = 0
        for line in file:
            line_p = [x for x in line[:-1].split(' ') if x]
            if len(line_p) < 4:
                continue
            if PALETTE.getName(line_p[0]) != 999 and all([True if x.replace('.', '').replace('-', '').isdigit() else False for x in line_p[1:4]]):
                i += 1
                coord = np.array([float(x) for x in line_p[1:]], dtype=np.float32)
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(line_p[0]), parent=mol, name=f'{line_p[0]}{i}')
        return self.parsMolSys(mol_sys, bond, root)

    def parsPdb(self, file_path, bond=True, root=None, *args, **kwargs):
        mol_list, atom_list, bonds_l = None, None, None
        file = open(file_path, 'r')
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol_sys.file_name = file_path
        mol = MoleculeClass.Molecule(parent=mol_sys)

        def pointCreation(atom_list, atom):
            coord = atom.coord.copy()
            point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                      color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                      name=atom.name,
                                      atom_type=atom.atom_type,
                                      pdb_flag=atom.pdb_flag,
                                      pdb_atom_seq=atom.pdb_atom_seq,
                                      pdb_name=atom.pdb_name,
                                      pdb_alt_loc=atom.pdb_alt_loc,
                                      pdb_res_name=atom.pdb_res_name,
                                      pdb_chain=atom.pdb_chain,
                                      pdb_res_seq=atom.pdb_res_seq,
                                      pdb_iCode=atom.pdb_Achar,
                                      pdb_occupancy=atom.pdb_occupancy,
                                      pdb_tempFactor=atom.tempFactor,
                                      pdb_charge=atom.pdb_charge,
                                      label=atom.pdb_name)
            return point

        for line in file:
            flag = line[0:6].replace(' ', '')

            if flag == 'ATOM' or flag == 'HETATM':
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=np.float32)
                atom_type = line[76:78].replace(' ', '').capitalize()
                try:
                    atom_seq = int(line[6:11].replace(" ", ""))
                except ValueError:
                    atom_seq = -1
                atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(atom_type), parent=mol,
                                          name=f'{int(line[22:26]):>4}-{line[17:20]}{line[6:11].replace(" ", ""): >5}-{line[12:16].replace(" ", "")}',
                                          pdb_flag=flag,
                                          pdb_atom_seq=atom_seq,
                                          pdb_name=line[12:16].replace(" ", ""),
                                          pdb_alt_loc=line[16],
                                          pdb_res_name=line[17:20],
                                          pdb_chain=line[21],
                                          pdb_res_seq=int(line[22:26]),
                                          pdb_iCode=line[26],
                                          pdb_occupancy=float(line[54:60]),
                                          pdb_tempFactor=float(line[60:66]),
                                          pdb_charge=line[78:80])
                mol.addChild(atom)
        return self.parsMolSys(mol_sys, bond, root, pointCreation)

    def parsCpProp(self, file_path, bond=False, root=None, *args, **kwargs):
        mol_list, atom_list, bonds_l = None, None, None
        file = open(file_path, 'r')
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol_sys.file_name = file_path
        mol = MoleculeClass.Molecule(parent=mol_sys)

        types = {'3,-3': 6,
                 '3,-1': 7,
                 '3,+1': 8,
                 '3,+3': 9}
        attrs = [('Position (Bohr)', 'coord'),
                 ('Density of all electrons', 'density_of_all_electrons'),
                 ('Density of Alpha electrons', 'density_of_alpha_electrons'),
                 ('Density of Beta electrons', 'density_of_beta_electrons'),
                 ('Spin density of electrons', 'spin_density_of_electrons'),
                 ('Lagrangian kinetic energy G(r)', 'lagrangian_kinetic_energy_G_r'),
                 ('G(r) in X,Y,Z', 'G(r)_xyz'),
                 ('Hamiltonian kinetic energy K(r)', 'hamiltonian_kinetic_energy_K_r'),
                 ('Potential energy density V(r)', 'potential_energy_density_V_r'),
                 ('Energy density E(r) or H(r)', 'energy_density_E_r_H_r'),
                 ('Laplacian of electron density', 'laplacian_of_electron_density'),
                 ('Electron localization function (ELF)', 'electron_localization_function_ELF'),
                 ('Localized orbital locator (LOL)', 'localized_orbital_locator_LOL'),
                 ('Local information entropy', 'local_information_entropy'),
                 ('Reduced density gradient (RDG)', 'reduced_density_gradient_RDG'),
                 ('Reduced density gradient with promolecular approximation', 'reduced_density_gradient_with_promolecular_approximation'),
                 ('Sign(lambda2)*rho', 'sign_lambda2_rho'),
                 ('Sign(lambda2)*rho with promolecular approximation', 'sign_lambda2_rho_with_promolecular_approximation'),
                 ('Corr. hole for alpha, ref.', 'corr_hole_for_alpha_ref'),
                 ('Source function, ref.', 'source_function_ref'),
                 ('Wavefunction value for orbital         1 ', 'wavefunction_value_for_orbital_1'),
                 ('Average local ionization energy (ALIE)', 'average_local_ionization_energy_ALIE'),
                 ('Delta_g (under promolecular approximation)', 'delta_g_under_promolecular_approximation'),
                 ('Delta_g (under Hirshfeld partition)', 'delta_g_under_hirshfeld_partition'),
                 ('User-defined real space function', 'user_defined_real_space_function'),
                 ('ESP from nuclear charges', 'ESP_from_nuclear_charges')]
        spec_attrs = [('Components of gradient in x/y/z are', 'electron_density_gradient_xyz'),
                      ('Norm of gradient is', 'norm_of_gradient'),
                      ('Components of Laplacian in x/y/z are', 'laplacian_xyz'),
                      ('Total', 'norm_of_laplacian'),
                      ('Hessian matrix', 'hessian_matrix'),
                      ('Eigenvalues of Hessian', 'eigenvalues_of_hessian'),
                      ('Eigenvectors(columns) of Hessian', 'eigenvectors_of_hessian'),
                      ('Determinant of Hessian', 'determinant_of_Hessian'),
                      ('Ellipticity of electron density', 'ellipticity_of_electron_density'),
                      ('eta index', 'eta_index')]

        def parsCp(file, num, type):

            def pars_attr_line(file, line: str, i):
                nonlocal cp_attrs
                nonlocal coord
                if i == 0:
                    coord = np.array([float(x) for x in line[line.find(':')+1:-1].split() if x], dtype=np.float32)/1.8897
                elif i == 6:
                    value = np.array([float(x) for x in line[line.find(':')+1:-1].split() if x], dtype=np.float32)
                    cp_attrs[attrs[i][1]] = value
                elif i == 18 or i == 19:
                    loc = line.find(':')
                    loc2 = line.find(':', loc+1)
                    value = np.array([float(x) for x in line[loc+1:loc2].split() if x], dtype=np.float32)
                    value2 = [x for x in line[loc2+1:-1].split() if x][0]
                    if value2 == 'NaN':
                        value2 = None
                    else:
                        value2 = float(value2)
                    value = [value, value2]
                    cp_attrs[attrs[i][1]] = value
                else:
                    value = [float(x) for x in line[line.find(':')+1:-1].split() if x][0]
                    cp_attrs[attrs[i][1]] = value

            def pars_spec_attr_line(file, line: str, i):
                nonlocal cp_spec_attrs
                if i == 0 or i == 2:
                    line = file.__next__()
                    value = np.array([float(x) for x in line[:-1].split() if x], dtype=np.float32)
                    cp_spec_attrs[spec_attrs[i][1]] = value
                elif i == 4 or i == 6:
                    matr = []
                    for _ in range(3):
                        line = file.__next__()
                        value = [float(x) for x in line[:-1].split() if x]
                        matr.append(value)
                    value = np.array(matr, dtype=np.float32)
                    cp_spec_attrs[spec_attrs[i][1]] = value
                elif i == 5:
                    value = np.array([float(x) for x in line[line.find(':')+1:-1].split() if x], dtype=np.float32)
                    cp_spec_attrs[spec_attrs[i][1]] = value
                else:
                    value = [float(x) for x in line[line.find(':')+1:-1].split() if x][0]
                    cp_spec_attrs[spec_attrs[i][1]] = value

            cp_attrs = {x[1]: None for x in attrs[1:]}
            cp_attrs['type'] = type
            attrs_len = len(attrs)
            cp_spec_attrs = {x[1]: None for x in spec_attrs}
            spec_attrs_len = len(spec_attrs)
            pars = pars_attr_line
            coord = None
            atom_type = types[type]
            i = 0
            attrs_comp = attrs
            for line in file:
                if attrs_comp[i][0] in line:
                    pars(file, line, i)
                    i += 1
                    if i == attrs_len and pars is pars_attr_line:
                        pars = pars_spec_attr_line
                        attrs_comp = spec_attrs
                        i = 0
                    elif i == spec_attrs_len and pars is pars_spec_attr_line:
                        if file.__next__() == ' \n':
                            break
                        else:
                            raise ValueError
            cp_attrs.update(cp_spec_attrs)
            point = MoleculeClass.Atom(coord, atom_type, parent=mol, **cp_attrs)
            return point

        for line in file:
            if '----------------' in line:
                line_ed = line.split()
                line_ed = [x for x in line_ed if x]
                num = int(line_ed[2][:-1])
                type = line_ed[4][1:-1]
                parsCp(file, num, type)
        return mol_sys, (mol_list, atom_list, bonds_l)

    def parsCif(self, file_path, bond=True, root=None, *args, **kwargs):
        from gemmi import cif

        mol_list, atom_list, bonds_l = None, None, None
        block = cif.read_file(file_path).sole_block()
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.file_name = file_path
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol = MoleculeClass.Molecule(parent=mol_sys)

        sym_codes = block.find(['_symmetry_equiv_pos_as_xyz'])
        if not sym_codes:
            sym_codes = block.find(['_space_group_symop_operation_xyz'])
            sym_codes = [[i, x[0][1:-1]] for i, x in enumerate(sym_codes)]
        else:
            sym_codes = [[i, x[0]] for i, x in enumerate(sym_codes)]
        num = sym_codes[-1][0] + 1
        added = [
                 [num, 'x+1,y,z'],
                 [num+1, 'x-1,y,z'],
                 [num+2, 'x,y+1,z'],
                 [num+3, 'x,y-1,z'],
                 [num+4, 'x,y,z+1'],
                 [num+5, 'x,y,z-1'],
                ]
        sym_codes += added

        cell = block.find(['_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'])[0]
        cell = [float(cell[x][:cell[x].find('(')]) if cell[x].find('(') != -1 else float(cell[x]) for x in range(len(cell))]

        atoms = block.find(['_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z'])
        atoms = [[x[i] if i < 2 else float(x[i]) if x[i].find('(') == -1 else float(x[i][:x[i].find('(')]) for i in range(len(x))] for x in atoms]
        coords = [x[2:] for x in atoms]
        args = []
        args += cell
        args.append(coords)
        dec_coords = self.fracToDec(*args)
        for i, atom in enumerate(atoms):
            cif_data = {'cif_sym_codes': sym_codes,
                        'cif_cell_a': cell[0],
                        'cif_cell_b': cell[1],
                        'cif_cell_c': cell[2],
                        'cif_cell_al': cell[3],
                        'cif_cell_be': cell[4],
                        'cif_cell_ga': cell[5],
                        'cif_frac_coords': np.array(atom[2:], dtype=np.float32)}
            coord = np.array(dec_coords[i], dtype=np.float32)
            atom = MoleculeClass.Atom(coord.copy(), PALETTE.getName(atom[1]), parent=mol, name=atom[0], **cif_data)
        cell_coords = [[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1],
                       [0, 1, 1],
                       [1, 1, 0],
                       [1, 0, 1],
                       [1, 1, 1]]
        args = [*cell, cell_coords]
        cell_dec_coords = self.fracToDec(*args)

        mol_sys, list_tuple = self.parsMolSys(mol_sys, bond, root)
        mol_list = list_tuple[0]
        cell_list = point_class.PointsList(parent=mol_list, name='Cell', color=[0, 0, 0, 1])
        o = point_class.Point(parent=cell_list, name='o', coord=cell_dec_coords[0], color=[0, 0, 0, 1])
        a = point_class.Point(parent=cell_list, name='a', coord=cell_dec_coords[1], color=[1, 0, 0, 1])
        b = point_class.Point(parent=cell_list, name='b', coord=cell_dec_coords[2], color=[0, 1, 0, 1])
        c = point_class.Point(parent=cell_list, name='c', coord=cell_dec_coords[3], color=[0, 0, 1, 1])
        cell_rend = point_class.PointsList(parent=cell_list, name='Cell render', color=cell_list)

        oa = point_class.PointsList(parent=cell_rend, name='o-a', color=a)
        point1, point2 = (point_class.Point(parent=oa, coord=o, color=oa), point_class.Point(parent=oa, coord=a, color=oa))

        ob = point_class.PointsList(parent=cell_rend, name='o-b', color=b)
        point1, point2 = (point_class.Point(parent=ob, coord=o, color=ob), point_class.Point(parent=ob, coord=b, color=ob))

        oc = point_class.PointsList(parent=cell_rend, name='o-c', color=c)
        point1, point2 = (point_class.Point(parent=oc, coord=o, color=oc), point_class.Point(parent=oc, coord=c, color=oc))

        ab = point_class.PointsList(parent=cell_rend, name='a-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ab, coord=a, color=ab), point_class.Point(parent=ab, coord=cell_dec_coords[5], color=ab))
        ac = point_class.PointsList(parent=cell_rend, name='a-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ac, coord=a, color=ac), point_class.Point(parent=ac, coord=cell_dec_coords[6], color=ac))

        bc = point_class.PointsList(parent=cell_rend, name='b-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=bc, coord=b, color=bc), point_class.Point(parent=bc, coord=cell_dec_coords[4], color=bc))
        ba = point_class.PointsList(parent=cell_rend, name='b-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ba, coord=b, color=ba), point_class.Point(parent=ba, coord=cell_dec_coords[5], color=ba))

        cb = point_class.PointsList(parent=cell_rend, name='c-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=cb, coord=c, color=cb), point_class.Point(parent=cb, coord=cell_dec_coords[4], color=cb))
        ca = point_class.PointsList(parent=cell_rend, name='c-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ca, coord=c, color=ca), point_class.Point(parent=ca, coord=cell_dec_coords[6], color=ca))

        abc = point_class.PointsList(parent=cell_rend, name='ab-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=abc, coord=cell_dec_coords[7], color=abc), point_class.Point(parent=abc, coord=cell_dec_coords[5], color=abc))
        acb = point_class.PointsList(parent=cell_rend, name='ac-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=acb, coord=cell_dec_coords[7], color=acb),
                          point_class.Point(parent=acb, coord=cell_dec_coords[6], color=acb))
        bca = point_class.PointsList(parent=cell_rend, name='bc-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=bca, coord=cell_dec_coords[7], color=bca),
                          point_class.Point(parent=bca, coord=cell_dec_coords[4], color=bca))

        return mol_sys, list_tuple

    def parsMolSys(self, mol_sys, bond=True, root=None, point_func=None, *args, **kwargs):
        if point_func is None:
            point_func = self._pointCreation
        mol_list, atom_list, bonds_l = None, None, None
        try:
            mol = mol_sys.children[0]
        except IndexError:
            return mol_sys, (None, None, None)
        if bond:
            mol_sys.genBonds()
        if root is not None:
            mol_list = point_class.PointsList(parent=root, name=mol_sys.name)
            mol.assignPoint(mol_list)
            MOLECULE_SYSTEMS[mol_list] = mol_sys
            atom_list = point_class.PointsList(parent=mol_list, rad=0.25, name='Atoms')
            i = 1
            for atom in mol.children:
                coord = atom.coord.copy()
                point = point_func(atom_list, atom)
                atom.assignPoint(point)
                i += 1
            if bond:
                bonds_l = point_class.PointsList(parent=mol_list, rad=0.1, name='Bonds')
                bonds = []
                for atom in mol:
                    for bond in atom.bonds():
                        if bond not in bonds:
                            bond_l = point_class.PointsList(parent=bonds_l, name=f'{bond.parents()[0].point().name}_{bond.parents()[1].point().name}', rad=bonds_l)
                            b1 = point_class.Point(coord=bond.parents()[0].point(), color=bond.parents()[0].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            b2 = point_class.Point(coord=bond.parents()[1].point(), color=bond.parents()[1].point(),
                                                   rad=bond_l,
                                                   parent=bond_l)
                            bonds.append(bond)
        return mol_sys, (mol_list, atom_list, bonds_l)

    def parsWinxproOut(self, file_path, bond=True, root=None, *args, **kwarg):

        def parseCell(line, file):
            nonlocal flag
            if 'alpha          beta           gamma' in line:
                flag = True
                return
            if flag:
                par = [float(x) for x in line[:-1].split(' ') if x]
                flag = False
                methods.pop(0)
            return

        def parseAtm(line, file):
            nonlocal flag
            if 'Symop  Z       xf      yf      zf         x      y      z' in line:
                flag = True
                return
            if flag:
                if line == ' \n':
                    flag = False
                    methods.pop(0)
                    return
                data = [x for x in line[:-1].split(' ') if x]
                name = data[0]
                atom_type = int(data[1])
                coords = [float(x) for x in data[5:8]]
                atom = MoleculeClass.Atom(coords, atom_type, None, name=name)
                atoms.append(atom)

        def parseCPs(line, file):
            types = {'3,-3': 6,
                     '3,-1': 7,
                     '3,+1': 8,
                     '3,+3': 9}
            nonlocal flag
            if ' Summary of CPs: ' in line:
                flag = True
                file.__next__()
                file.__next__()
                return
            if flag:
                def checkForEnd():
                    nonlocal flag
                    if line == '\n':
                        flag = False
                        methods.pop(0)
                        return
                def parseCP():
                    nonlocal line
                    data = [x for x in line[:-1].split(' ') if x]
                    type = data[1][1:-1]
                    coords =[float(x) for x in data[8:11]]
                    cp = MoleculeClass.Atom(coords, types.get(type, 999))
                    cps.append(cp)
                    file.__next__()
                    file.__next__()
                    file.__next__()
                    file.__next__()
                    file.__next__()
                    file.__next__()
                    line = file.__next__()
                    checkForEnd()
                    return
                checkForEnd()
                while flag:
                    parseCP()
                return


        mol_list, atom_list, bonds_l = None, None, None
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.file_name = file_path
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol = MoleculeClass.Molecule(parent=mol_sys)

        cps_sys = MoleculeClass.MoleculeSystem()
        cps_sys.file_name = file_path
        cps_sys.name = 'CPs'
        cps = MoleculeClass.Molecule(parent=cps_sys)

        flag = False
        cell = None
        atoms = []
        cps = []
        methods = [parseCell, parseAtm, parseCPs]

        with open(file_path, 'r') as file:
            for line in file:
                if methods:
                    methods[0](line, file)
                else:
                    break
        mol: MoleculeClass.Molecule
        mol = mol_sys.children[0]
        for atom in atoms:
            mol.addChild(atom)

        mol = cps_sys.children[0]
        for cp in cps:
            mol.addChild(cp)

        mol_sys, list_tuple = self.parsMolSys(mol_sys, bond, root)
        return mol_sys, list_tuple

    @staticmethod
    def _pointCreation(atom_list, atom):
        coord = atom.coord.copy()
        point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                  color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                  atom_type=atom.atom_type,
                                  name=atom.name,
                                  label=atom.name)
        return point


PARSER = FileParser()