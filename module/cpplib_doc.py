"""
The shell for cpplib module with documentation of available functions.
"""

import cpplib
from typing import List


def search_main(graph: str, data: List[str], procs: int, exact: bool) -> List[int]:
    """
    Search 'graph' in 'data' with graphs.
    Variables:
        graph:  graph of a molecule
                ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
        data:   the list of molecules graphs for search in
        procs:  the number of physical processors for multiprocessing
        exact:  is the exact search (True) or substructure (False)
    """
    return cpplib.SearchMain(graph, data, procs, exact)


def compare_graph(graph_1: str, graph_2: str, exact: bool):
    return cpplib.CompareGraph(graph_1, graph_2, exact)


def find_molecules_in_cell(cell_params: List[float], symops: List[str], atom_types: List[int], xyz: List[float]):
    return cpplib.FindMoleculesInCell(cell_params, symops, atom_types, xyz)


def find_molecules_without_cell(atom_types: List[int], xyz: List[float]):
    return cpplib.FindMoleculesWithoutCell(atom_types, xyz)


def gen_bonds(atoms: List[List[int | float]]):
    return cpplib.GenBonds(atoms)


def gen_symm(atoms: List[List[int | float]], symop: str):
    return cpplib.GenSymm(atoms, symop)
