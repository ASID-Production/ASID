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
                                  '.cif': self.parsCif,
                                  'POSCAR*': self.parsPOSCAR_CONTCAR,
                                  'CONTCAR*': self.parsPOSCAR_CONTCAR}

    @staticmethod
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
        basename, ext = os.path.splitext(file_path)
        basename = os.path.basename(basename)
        ext = ext.lower()
        if ext in self.SUPPORTED_FORMATS:
            return self.SUPPORTED_FORMATS[ext](file_path, bond, root)
        else:
            e = [x for x in self.SUPPORTED_FORMATS if x.replace('*', '') in basename]
            if e:
                return self.SUPPORTED_FORMATS[e[0]](file_path, bond, root)
            else:
                return (None, None)

    def parsXyz(self, file_path, bond=True, root=None, *args, **kwargs):
        from . import Db_viewer

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

        ret = self.parsMolSys(mol_sys, bond, root)

        mol_list = ret[1][0]
        if mol_list.additional_context_actions is None:
            mol_list.addProperty('additional_context_actions',
                                 [('Export 2d diagram', lambda: Db_viewer.exportGif(file=mol_sys.file_name))],
                                 observed=False)
        else:
            mol_list.additional_context_actions.append(
                ('Export 2d diagram', lambda: Db_viewer.exportGif(file=mol_sys.file_name)))

        return ret

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
        from . import Db_viewer

        mol_list, atom_list, bonds_l = None, None, None
        blocks = cif.read_file(file_path)
        for block in blocks:
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

            if mol_list.additional_context_actions is None:
                mol_list.addProperty('additional_context_actions', [('Export 2d diagram', lambda: Db_viewer.exportGif(file=mol_sys.file_name))], observed=False)
            else:
                mol_list.additional_context_actions.append(('Export 2d diagram', lambda: Db_viewer.exportGif(file=mol_sys.file_name)))

            cell_list = self.createCellList(mol_list, cell_dec_coords)

            yield mol_sys, list_tuple

    def createCellList(self, parent, cell_dec_coords):
        cell_list = point_class.PointsList(parent=parent, name='Cell', color=[0, 0, 0, 1])
        o = point_class.Point(parent=cell_list, name='o', coord=cell_dec_coords[0], color=[0, 0, 0, 1], label='o')
        a = point_class.Point(parent=cell_list, name='a', coord=cell_dec_coords[1], color=[1, 0, 0, 1], label='a')
        b = point_class.Point(parent=cell_list, name='b', coord=cell_dec_coords[2], color=[0, 1, 0, 1], label='b')
        c = point_class.Point(parent=cell_list, name='c', coord=cell_dec_coords[3], color=[0, 0, 1, 1], label='c')
        cell_rend = point_class.PointsList(parent=cell_list, name='Cell render', color=cell_list)

        oa = point_class.PointsList(parent=cell_rend, name='o-a', color=a)
        point1, point2 = (
        point_class.Point(parent=oa, coord=o, color=oa), point_class.Point(parent=oa, coord=a, color=oa))

        ob = point_class.PointsList(parent=cell_rend, name='o-b', color=b)
        point1, point2 = (
        point_class.Point(parent=ob, coord=o, color=ob), point_class.Point(parent=ob, coord=b, color=ob))

        oc = point_class.PointsList(parent=cell_rend, name='o-c', color=c)
        point1, point2 = (
        point_class.Point(parent=oc, coord=o, color=oc), point_class.Point(parent=oc, coord=c, color=oc))

        ab = point_class.PointsList(parent=cell_rend, name='a-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ab, coord=a, color=ab),
                          point_class.Point(parent=ab, coord=cell_dec_coords[5], color=ab))
        ac = point_class.PointsList(parent=cell_rend, name='a-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ac, coord=a, color=ac),
                          point_class.Point(parent=ac, coord=cell_dec_coords[6], color=ac))

        bc = point_class.PointsList(parent=cell_rend, name='b-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=bc, coord=b, color=bc),
                          point_class.Point(parent=bc, coord=cell_dec_coords[4], color=bc))
        ba = point_class.PointsList(parent=cell_rend, name='b-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ba, coord=b, color=ba),
                          point_class.Point(parent=ba, coord=cell_dec_coords[5], color=ba))

        cb = point_class.PointsList(parent=cell_rend, name='c-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=cb, coord=c, color=cb),
                          point_class.Point(parent=cb, coord=cell_dec_coords[4], color=cb))
        ca = point_class.PointsList(parent=cell_rend, name='c-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=ca, coord=c, color=ca),
                          point_class.Point(parent=ca, coord=cell_dec_coords[6], color=ca))

        abc = point_class.PointsList(parent=cell_rend, name='ab-c', color=cell_rend)
        point1, point2 = (point_class.Point(parent=abc, coord=cell_dec_coords[7], color=abc),
                          point_class.Point(parent=abc, coord=cell_dec_coords[5], color=abc))
        acb = point_class.PointsList(parent=cell_rend, name='ac-b', color=cell_rend)
        point1, point2 = (point_class.Point(parent=acb, coord=cell_dec_coords[7], color=acb),
                          point_class.Point(parent=acb, coord=cell_dec_coords[6], color=acb))
        bca = point_class.PointsList(parent=cell_rend, name='bc-a', color=cell_rend)
        point1, point2 = (point_class.Point(parent=bca, coord=cell_dec_coords[7], color=bca),
                          point_class.Point(parent=bca, coord=cell_dec_coords[4], color=bca))
        return cell_list

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

    def parsDict(self, dict, root=None, point_func=None, name='None', *args, **kwargs):
        if point_func is None:
            point_func = self._pointCreation
        mol_list = None
        ret = None

        def rec_p(val, parent):
            if type(val) is list:
                for item in val:
                    point = point_func(parent, item)
                    item.assignPoint(point)
                return
            else:
                for key in val:
                    par = point_class.PointsList(parent=parent, rad=parent, name=key)
                    rec_p(val[key], par)

        if root is not None:
            ret = point_class.PointsList(parent=root, rad=0.25, name=name)
            rec_p(dict, ret)
        return dict, (ret, None, None)

    def parsWinxproOut(self, file_path, bond=True, root=None, *args, **kwargs):

        def parseCell(line, file):
            nonlocal flag
            if 'alpha          beta           gamma' in line:
                flag = True
                return
            if flag:
                par = [float(x) for x in line[:-1].split(' ') if x]
                cif_dict['cif_cell_a'] = par[0]
                cif_dict['cif_cell_b'] = par[1]
                cif_dict['cif_cell_c'] = par[2]
                cif_dict['cif_cell_al'] = par[3]
                cif_dict['cif_cell_be'] = par[4]
                cif_dict['cif_cell_ga'] = par[5]
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
                frac_coords = np.array([float(x) for x in data[2:5]], dtype=np.float32)
                atom = MoleculeClass.Atom(coords, atom_type, None, name=name, cif_frac_coords=frac_coords)
                atoms.append(atom)

        def parseSymm(line, file):
            nonlocal flag
            if 'CENTROSYMMETRIE =' in line:
                flag = True
                return
            if flag:
                if line == ' \n':
                    flag = False
                    methods.pop(0)
                    return
                data = [x for x in line[:-1].split(' ') if x]
                if data[0] == 'New':
                    return
                else:
                    i = 0
                    symop = ['', 'x', '', 'y', '', 'z']
                    if 'X' not in data[i]:
                        symop[0] = data[i]
                        i += 1
                    op = data[i]
                    symop[1] = op[:op.find('X')] + 'x' + op[op.find('X')+1:]
                    i += 1

                    if 'Y' not in data[i]:
                        symop[2] = data[i]
                        i += 1
                    op = data[i]
                    symop[3] = op[:op.find('Y')] + 'y' + op[op.find('Y') + 1:]
                    i += 1

                    if 'Z' not in data[i]:
                        symop[4] = data[i]
                        i += 1
                    op = data[i]
                    symop[5] = op[:op.find('Z')] + 'z' + op[op.find('Z') + 1:]
                    i += 1
                    code = f'{symop[0]}{symop[1]},{symop[2]}{symop[3]},{symop[4]}{symop[5]}'
                    symmcodes.append(code)

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
                    name = data[0]
                    type = data[1][1:-1]
                    coords =[float(x) for x in data[8:11]]
                    prop_d = {}
                    prop_d['name'] = name
                    conn = f'{data[2]}-{data[5]}'
                    prop_d['conn'] = conn
                    line = file.__next__()
                    data = [x for x in line[:-1].split(' ') if x]
                    prop_d['rho'] = float(data[2])
                    prop_d['lambda_1'] = float(data[6])
                    line = file.__next__()
                    data = [x for x in line[:-1].split(' ') if x]
                    prop_d['d2rho'] = float(data[2])
                    prop_d['lambda_2'] = float(data[6])
                    line = file.__next__()
                    data = [x for x in line[:-1].split(' ') if x]
                    prop_d['g'] = float(data[2])
                    prop_d['lambda_3'] = float(data[6])
                    line = file.__next__()
                    data = [x for x in line[:-1].split(' ') if x]
                    prop_d['v'] = float(data[2])
                    prop_d['ellipticity'] = float(data[6])
                    line = file.__next__()
                    data = [x for x in line[:-1].split(' ') if x]
                    prop_d['h'] = float(data[2])
                    cp = MoleculeClass.Atom(coords, types.get(type, 999), **prop_d)
                    cps.append(cp)
                    line = file.__next__()
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
        mol_sys.open_func = lambda: self.parsWinxproOut(file_path, bond=bond, root=root, *args, **kwargs)
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
        symmcodes = []
        methods = [parseCell, parseAtm, parseSymm, parseCPs]

        cif_dict = {'cif_sym_codes': [[0, 'x,y,z']],
                    'cif_cell_a': 1,
                    'cif_cell_b': 1,
                    'cif_cell_c': 1,
                    'cif_cell_al': 90,
                    'cif_cell_be': 90,
                    'cif_cell_ga': 90,}

        with open(file_path, 'r') as file:
            for line in file:
                if methods:
                    methods[0](line, file)
                else:
                    break

        tr_codes = ['x+1,y,z', 'x-1,y,z', 'x,y+1,z', 'x,y-1,z', 'x,y,z+1', 'x,y,z-1']
        symmcodes += tr_codes
        codes = [[i, symmcodes[i]] for i in range(len(symmcodes))]
        cif_dict['cif_sym_codes'] = codes

        mol: MoleculeClass.Molecule
        mol = mol_sys.children[0]
        for atom in atoms:
            atom.updateProps(cif_dict)
            mol.addChild(atom)

        mol = cps_sys.children[0]
        for cp in cps:
            mol.addChild(cp)

        def apCreation(atom_list, atom):
            point = self._pointCreation(atom_list, atom)
            cif = {key: atom.__getattribute__(key) for key in cif_dict.keys()}
            for key in cif:
                point.__setattr__(key, cif[key])
            return point

        mol_sys, list_tuple = self.parsMolSys(mol_sys, bond, root, point_func=apCreation)
        yield mol_sys, list_tuple

        def cppCreation(atom_list, atom):
            coord = atom.coord.copy()
            point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                      color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                      atom_type=atom.atom_type,
                                      name=atom.name,
                                      connectivity=atom.conn,
                                      label=atom.name,
                                      rho=atom.rho,
                                      d2rho=atom.d2rho,
                                      g=atom.g,
                                      v=atom.v,
                                      h=atom.h,
                                      lambda_1=atom.lambda_1,
                                      lambda_2=atom.lambda_2,
                                      lambda_3=atom.lambda_3,
                                      ellipticity=atom.ellipticity)
            return point

        cps_sys, list_tuple = self.parsMolSys(cps_sys, False, root, point_func=cppCreation)
        if list_tuple[1]:
            list_tuple[1].rad = 0.15
        yield cps_sys, list_tuple

    def parsWinxproPaths(self, file_path, bond=False, root=None, *args, **kwargs):
        mol_list, atom_list, bonds_l = None, None, None
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.file_name = file_path
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol = MoleculeClass.Molecule(parent=mol_sys)

        with open(file_path, 'r') as file:
            for line in file:
                data = [x for x in line[:-1].split() if x]
                try:
                    data[0], data[1] = float(data[0]), float(data[1])
                    data[2] = float(data[2][:-1])
                    atom = MoleculeClass.Atom(data[:3], 6, parent=mol)
                except ValueError or IndexError:
                    pass
        mol_sys, list_tuple = self.parsMolSys(mol_sys, bond, root)
        if list_tuple[1]:
            list_tuple[1].rad = 0.02
        yield mol_sys, list_tuple

    def parsAIMALLsumviz(self, file_path, bond=False, root=None, *args, **kwargs):

        def parseAtm(line, file):
            nonlocal flag
            if 'Atom      Charge                X                  Y                  Z' in line:
                flag = True
                file.__next__()
                return
            if flag:
                if line == ' \n':
                    flag = False
                    methods.pop(0)
                    return
                data = [x for x in line[:-1].split(' ') if x]
                name = data[0]
                atom_type = ''.join([x for x in name if not x.isdigit()])
                coords = [float(x) for x in data[2:5]]
                frac_coords = np.array([float(x) for x in data[2:5]], dtype=np.float32)
                atom = MoleculeClass.Atom(coords, atom_type, None, name=name, cif_frac_coords=frac_coords)
                atom.coord = atom.coord/1.88973
                atoms.append(atom)

        def parseCP(line, file):
            nonlocal flag
            nonlocal last_cp

            types = {'(3,-3)': 6,
                     '(3,-1)': 7,
                     '(3,+1)': 8,
                     '(3,+3)': 9}

            spec_props = ['GradRho', 'HessRho_EigVals', 'HessRho_EigVec1', 'HessRho_EigVec2', 'HessRho_EigVec3', 'Stress_EigVals', 'Stress_EigVec1', 'Stress_EigVec2', 'Stress_EigVec3']

            if 'CP#' in line and 'Coords' in line:
                def create_atom():
                    atom_type = types[data['Type']]
                    atom = MoleculeClass.Atom(coord, atom_type, None, name=cp_name, sup_data_dict=data)
                    atom.coord = atom.coord/1.88973
                    cps.append(atom)

                flag = True
                line = [x for x in line.split(' ') if x]
                coord = line[4:7]
                cp_name = line[1]
                last_cp = cp_name
                coord = [float(x) for x in coord]
                data = {}

                for line in file:
                    if line == '\n':
                        flag = False
                        create_atom()
                        return
                    if 'sample points along' in line:

                        f = False
                        for key in cps_dict:
                            if key in line:
                                f = True
                                dest = cps_dict[key].get(last_cp, None)
                                if dest is None:
                                    dest = []
                                    cps_dict[key][last_cp] = dest
                        if not f:
                            dest = cps_dict['other'].get(last_cp, None)
                            if dest is None:
                                dest = []
                                cps_dict['other'][last_cp] = dest

                        create_atom()
                        parsePath(line, file, dest)
                        flag = False
                        return
                    if 'Bond Ellipticity' in line:
                        line = [x for x in line.split(' ') if x]
                        name = 'Bond_Ellipticity'
                        value = line[3]
                        if 'NA' in value:
                            value = None
                        else:
                            value = float(line[3])
                        data[name] = value
                    elif '-DivStress' in line:
                        line = [x for x in line.split(' ') if x]
                        name = 'm_DivStress'
                        values = [float(x) if 'NA' not in x else None for x in line[2:5]]
                        data[name] = values
                    elif 'Type' in line:
                        line = [x for x in line.split(' ') if x]
                        name = line[0]
                        value = line[2]
                        if 'NA' in value:
                            value = None
                        else:
                            value = line[2]
                        data[name] = value
                    else:
                        check = [x in line for x in spec_props]
                        if any(check):
                            line = [x for x in line.split(' ') if x]
                            name = line[0]
                            values = [float(x) if 'NA' not in x else None for x in line[2:5]]
                            data[name] = values
                        else:
                            line = [x for x in line.split(' ') if x]
                            name = line[0]
                            value = line[2]
                            if 'NA' in value:
                                value = None
                            else:
                                value = float(line[2])
                            data[name] = value

            if flag:
                if line == ' \n':
                    flag = False
                    methods.pop(0)
                    return
                data = [x for x in line[:-1].split(' ') if x]
                name = data[0]
                atom_type = ''.join([x for x in name if not x.isdigit()])
                coords = [float(x) for x in data[2:5]]
                frac_coords = np.array([float(x) for x in data[2:5]], dtype=np.float32)
                atom = MoleculeClass.Atom(coords, atom_type, None, name=name, cif_frac_coords=frac_coords)
                atoms.append(atom)

        def parsePath(line, file, dest=None):
            nonlocal flag
            if dest is None:
                dest = []
            for line in file:
                if 'sample points along' in line:
                    if len(dest) % 2 != 0:
                        atom = MoleculeClass.Atom(dest[-1].coord, atom_type=6)
                        dest.append(atom)
                    f = False
                    for key in cps_dict:
                        if key in line:
                            f = True
                            dest = cps_dict[key].get(last_cp, None)
                            if dest is None:
                                dest = []
                                cps_dict[key][last_cp] = dest
                    if not f:
                        dest = cps_dict['other'].get(last_cp, None)
                        if dest is None:
                            dest = []
                            cps_dict['other'][last_cp] = dest
                elif line == '\n':
                    if len(dest) % 2 != 0:
                        atom = MoleculeClass.Atom(dest[-1].coord, atom_type=6)
                        dest.append(atom)
                    flag = False
                    return
                else:
                    line = line.split(' ')
                    coord = [float(x) for x in line[:-1] if x][0:3]
                    atom = MoleculeClass.Atom(coord, atom_type=6)
                    atom.coord = atom.coord/1.88973
                    dest.append(atom)
            return

        mol_list, atom_list, bonds_l = None, None, None
        mol_sys = MoleculeClass.MoleculeSystem()
        mol_sys.file_name = file_path
        mol_sys.open_func = lambda: self.parsWinxproOut(file_path, bond=bond, root=root, *args, **kwargs)
        mol_sys.name = os.path.basename(file_path).split('.')[0]
        mol = MoleculeClass.Molecule(parent=mol_sys)

        cps_sys = MoleculeClass.MoleculeSystem()
        cps_sys.file_name = file_path
        cps_sys.name = 'CPs'
        cp_mol = MoleculeClass.Molecule(parent=cps_sys)

        paths_sys = MoleculeClass.MoleculeSystem()
        paths_sys.file_name = file_path
        paths_sys.name = 'Paths'
        path_mol = MoleculeClass.Molecule(parent=paths_sys)

        flag = False
        methods = [parseAtm, parseCP]
        atoms = []
        cps = []
        cps_dict = {'first RCP attractor path': {},
                    'second RCP attractor path': {},
                    'path from RCP to BCP': {},
                    'path from BCP to atom': {},
                    'IAS +EV1 path from BCP': {},
                    'IAS -EV1 path from BCP': {},
                    'IAS +EV2 path from BCP': {},
                    'IAS -EV2 path from BCP': {},
                    'other': {}}
        last_cp = None
        paths = []

        with open(file_path, 'r') as file:
            for line in file:
                if methods:
                    methods[0](line, file)
                else:
                    break

        for atom in atoms:
            mol.addChild(atom)

        for cp in cps:
            cp_mol.addChild(cp)

        def createCpp(atom_list, atom):
            coord = atom.coord.copy()
            point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                      color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                      atom_type=atom.atom_type,
                                      name=atom.name,
                                      connectivity=None,
                                      label=atom.name,
                                      **atom.sup_data_dict)
            return point

        yield self.parsMolSys(mol_sys, bond, root)
        list_tuple = self.parsMolSys(cps_sys, False, root, createCpp)
        if list_tuple[1]:
            list_tuple[1][0].rad = 0.15
        yield list_tuple

        for path in paths:
            path_mol.addChild(path)
        list_tuple = self.parsDict(cps_dict, root, name='Paths')
        if list_tuple[1]:
            list_tuple[1][0].rad = 0.02
        yield list_tuple

    def parsPOSCAR_CONTCAR(self, file_path, bond=False, root=None, *args, **kwargs):
        from pymatgen.core import Structure

        def pymatgenStructureToMolSys(struct):
            mol_sys = MoleculeClass.MoleculeSystem()
            mol_sys.name = os.path.basename(file_path).split('.')[0]
            mol_sys.file_name = file_path
            mol = MoleculeClass.Molecule(parent=mol_sys)

            types = struct.atomic_numbers
            coords = struct.cart_coords
            names = struct.labels

            for i in range(len(types)):
                type = types[i]
                coord = coords[i]
                add_data = {'cif_sym_codes': [[0, 'x,y,z'],
                                              [1, 'x+1,y,z'],
                                              [2, 'x-1,y,z'],
                                              [3, 'x,y+1,z'],
                                              [4, 'x,y-1,z'],
                                              [5, 'x,y,z+1'],
                                              [6, 'x,y,z-1']],
                            'cif_cell_a': struct.lattice.a,
                            'cif_cell_b': struct.lattice.b,
                            'cif_cell_c': struct.lattice.c,
                            'cif_cell_al': struct.lattice.alpha,
                            'cif_cell_be': struct.lattice.beta,
                            'cif_cell_ga': struct.lattice.gamma,
                            'cif_frac_coords': np.array(struct.frac_coords[i], dtype=np.float32)}
                atom = MoleculeClass.Atom(coord, type, parent=mol, name=names[i], **add_data)
            return mol_sys

        def savePOSCAR_ToCif(molsys):
            from pymatgen.io.cif import CifWriter
            from PySide6 import QtWidgets
            file = molsys.file_name
            filename, _ = QtWidgets.QFileDialog.getSaveFileName(None, 'POSCAR to cif', filter='*.cif, *')
            w = CifWriter(Structure.from_file(file))
            w.write_file(filename)


        struct = Structure.from_file(file_path)
        mol_sys = pymatgenStructureToMolSys(struct)
        ret = self.parsMolSys(mol_sys, bond, root)

        cell_coords = [[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1],
                       [0, 1, 1],
                       [1, 1, 0],
                       [1, 0, 1],
                       [1, 1, 1]]
        args = (*struct.lattice.abc, *struct.lattice.angles, cell_coords)
        cell_dec_coords = self.fracToDec(*args)
        self.createCellList(ret[1][0], cell_dec_coords)
        a = 0

        mol_list = ret[1][0]
        if mol_list.additional_context_actions is None:
            mol_list.addProperty('additional_context_actions', [('Save to cif', lambda: savePOSCAR_ToCif(mol_sys))], observed=False)
        else:
            mol_list.additional_context_action.append([('Save to cif', lambda: savePOSCAR_ToCif(mol_sys))])
        a = 0
        return ret

    @staticmethod
    def _pointCreation(atom_list, atom):

        coord = atom.coord.copy()
        add_data = atom.sup_data_dict
        if add_data is None:
            point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                      color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                      atom_type=atom.atom_type,
                                      name=atom.name,
                                      label=atom.name)
        else:
            point = point_class.Point(parent=atom_list, coord=coord, rad=atom_list,
                                      color=PALETTE.point_dict[PALETTE.getName(atom.atom_type)],
                                      atom_type=atom.atom_type,
                                      name=atom.name,
                                      label=atom.name,
                                      **add_data)
        return point


PARSER = FileParser()