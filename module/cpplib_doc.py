"""
The shell for cpplib module with documentation of available functions.
"""

import cpplib
from typing import List


def search_main(graph: str, 
				data: List[str], 
				procs: int, 
				exact: bool) -> List[int]:
    """
    Search 'graph' in 'data' with graphs.
    Variables:
        graph:  graph of a molecule
                ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
        data:   the list of molecules graphs for search in
        procs:  the number of threads for multiprocessing
        exact:  is the exact search (True) or substructure (False)
    """
    return cpplib.SearchMain(graph, data, procs, exact)


def compare_graph(graph_1: str, 
				  graph_2: str, 
				  exact: bool) -> bool:
    """
    Compare 'graph_1' and 'graph_2'.
    Variables:
        graph_1:  graph of first molecule
                ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
        graph_2:  graph of second molecule
                ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
        exact:  is the exact search (True) or substructure (False)
    """
    return cpplib.CompareGraph(graph_1, graph_2, exact)


def find_molecules_in_cell(cell_params: List[float], 
						   symops: List[str], 
						   atom_types: List[int], 
						   xyz: List[float]) -> str:
    """
    Find all different graphs of molecules in a crystal.
    Variables:
        cell_params:  list of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
			(a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
		symops:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
		atom_types: list of atomic types
		xyz: list of atomic fractal coordinates - three to one atom
    """
    return cpplib.FindMoleculesInCell(cell_params, symops, atom_types, xyz)


def find_molecules_without_cell(atom_types: List[int], 
                                xyz: List[float]):
    """
    Find all different graphs of molecules in a cluster.
    Variables:
        atom_types: list of atomic types
		xyz: list of atomic cartesian (decart) coordinates - three to one atom
    """
    return cpplib.FindMoleculesWithoutCell(atom_types, xyz)


def gen_bonds(atoms: List[List[int | float]]):
    """
    Find all bonds in of molecules in a cluster.
    Variables:
        atoms: list of [list of atomic type and three real-space coordinate]
			   type is integer, coordinates are floating point numbers.
			   example [ [ 1 0.0 0.0 0.0 ], [9 0.5 0.5 0.5], ... ]
    """
    return cpplib.GenBonds(atoms)


def gen_symm(atoms: List[List[int | float]], symop: str):
    """
    Find all bonds in of molecules in a cluster.
    Variables:
        atoms: list of [list of atomic type and three real-space coordinate]
			   type is integer, coordinates are floating point numbers.
			   example [ [ 1 0.0 0.0 0.0 ], [9 0.5 0.5 0.5], ... ]
		symops:  SYMM-code for structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.	
    """
    return cpplib.GenSymm(atoms, symop)
