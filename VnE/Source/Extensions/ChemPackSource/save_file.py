from . import MoleculeClass
from ..ChemPack import PALETTE


class SaveFile:

    def getFormats(self):
        return SaveFile.__formats

    def save(self, molecule, filename, format):
        meth = SaveFile.__formats.get(format, None)
        if meth is not None:
            meth(self, molecule, filename)
        pass

    def save_xyz(self, molecule, filename):
        mol = molecule
        while not isinstance(mol.children[0], MoleculeClass.Atom):
            mol = molecule.children[0]

        out = open(filename, 'w')
        for atom in mol.children:
            line = f'{atom.atom_type}{atom.coord[0]:> 19.9f}{atom.coord[1]:> 17.9f}{atom.coord[2]:> 17.9f}\n'
            out.write(line)
        out.close()
        return

    def save_pdb(self, molecule, filename):
        mol = molecule
        while not isinstance(mol.children[0], MoleculeClass.Atom):
            mol = molecule.children[0]

        out = open(filename, 'w')
        i = 1
        prev_atom = mol.children[0]
        for atom in mol.children:
            if prev_atom.pdb_chain != atom.pdb_chain or prev_atom.pdb_flag != atom.pdb_flag:
                ter_line = f'TER   {i:>5}{prev_atom.pdb_res_name:>9} {prev_atom.pdb_chain}{prev_atom.pdb_res_seq:>4}{prev_atom.pdb_iCode}\n'
                out.write(ter_line)
                i += 1
            line = f'{atom.pdb_flag:<6}{i:>5} {atom.pdb_name:<4}{atom.pdb_alt_loc}{atom.pdb_res_name} {atom.pdb_chain}{atom.pdb_res_seq:>4}{atom.pdb_iCode}{atom.coord[0]:> 11.3f}{atom.coord[1]:> 8.3f}{atom.coord[2]:> 8.3f}{atom.pdb_occupancy:> 6.2f}{atom.pdb_tempFactor:> 6.2f}{PALETTE.getName(atom.atom_type):>12}{atom.pdb_charge}\n'
            out.write(line)
            prev_atom = atom
            i += 1
        out.write('END\n')
        out.close()
        return

    __formats = {'pdb': save_pdb,
                 'xyz': save_xyz}


SAVE_FILE = SaveFile()
