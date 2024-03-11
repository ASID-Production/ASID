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
      graph: graph of a molecule
        ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
      data: the list of molecules graphs for search in
      procs: the number of threads for multiprocessing
      exact: is the exact search (True) or substructure (False)
    """
    return cpplib.SearchMain(graph, data, procs, exact)


def compare_graph(graph_1: str, 
				  graph_2: str, 
				  exact: bool) -> bool:
    """
    Compare 'graph_1' and 'graph_2'.
    Variables:
      graph_1: graph of first molecule
        ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
      graph_2: graph of second molecule
        ('idx num_atoms num_bonds atom h_num atom h_num ... atom1_bond1 atom2_bond1 ... multitypes')
      exact: is the exact search (True) or substructure (False)
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
        example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
    """
    return cpplib.GenBonds(atoms)

def gen_bonds_ex(atoms: List[List[int | float]]):
    """
    Find all bonds with lengths in of molecules in a cluster.
    Variables:
      atoms: list of [list of atomic type and three real-space coordinate]
      type is integer, coordinates are floating point numbers.
        example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
    """
    return cpplib.GenBonds(atoms)


def find_distances_cell(
        cell_params: List[float], 
        symops: List[str], 
        atom_types: List[int],
        xyz: List[float],
        params: List[int | float]) -> str:
    """
    Find all distances with current params in a cluster.
    Variables:
      cell_params:  list of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
        (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
      symops:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
      atom_types: list of atomic types
      xyz: list of atomic fractal coordinates - three to one atom
      params: int type1, int type2, float min12, float max12
        min/max values could be zero - in this case bond's values will be used 
    """
    return cpplib.FindDistanceIC(cell_params, symops, atom_types, xyz, params)


def find_distances_xyz(    atom_types: List[int], 
						   xyz: List[float],
                           params: List[int | float]) -> str:
    """
    Find all distances with current params in a crystal.
    Variables:
	  atom_types: list of atomic types
	  xyz: list of atomic cartesian coordinates - three to one atom
      params: int type1, int type2, float min12, float max12
        min/max values could be zero - in this case bond's values will be used 
    """
    return cpplib.FindDistanceWC(atom_types, xyz, params)


def find_angle_cell(
        cell_params: List[float], 
        symops: List[str], 
        atom_types: List[int],
        xyz: List[float],
        params: List[int | float]) -> str:
    """
    Find all distances with current params in a cluster.
    Variables:
      cell_params:  list of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
        (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
      symops:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
      atom_types: list of atomic types
      xyz: list of atomic fractal coordinates - three to one atom
      params: int type1, int type2, int type3,                     - types
              float min12, float max12, float min23, float max23,  - distances min/max
              float min123, float max123                           - angle min/max (in deg.)
        totally 9 parameters should be defined
        angles are calculated in the range (-180, +180] degrees.
        min/max of distances could be zero - in this case bond's values will be used 
    """
    return cpplib.FindAngleIC(cell_params, symops, atom_types, xyz, params)


def find_angle_xyz( 	   atom_types: List[int], 
						   xyz: List[float],
                           params: List[int | float]) -> str:
    """
    Find all distances with current params in a crystal.
    Variables:
	  atom_types: list of atomic types
	  xyz: list of atomic cartesian coordinates - three to one atom
      params: int type1, int type2, int type3,                     - types
              float min12, float max12, float min23, float max23,  - distances min/max
              float min123, float max123                           - angle min/max (in deg.)
        totally 9 parameters should be defined
        angles are calculated in the range (-180, +180] degrees.
        min/max of distances could be zero - in this case bond's values will be used 
    """
    return cpplib.FindAngleWC(atom_types, xyz, params)


def find_torsion_cell(
        cell_params: List[float], 
        symops: List[str], 
        atom_types: List[int],
        xyz: List[float],
        params: List[int | float]) -> str:
    """
    Find all distances with current params in a cluster.
    Variables:
      cell_params:  list of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
        (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
      symops:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
      atom_types: list of atomic types
      xyz: list of atomic fractal coordinates - three to one atom
      params: int type1, int type2, int type3, int type4 -- types
              float min12, float max12, float min23, float max23, float min34, float max34, -- distances min/max
              float min123, float max123, float min234, float max234, -- angles min/max (in deg.)
              float min1234, float max1234 -- torsion min/max (in deg.)
        totally 16 parameters should be defined
        angles and torsions are calculated in the range (-180, +180] degrees.
        min/max of distances could be zero - in this case bond's values will be used 
    """
    return cpplib.FindTorsionIC(cell_params, symops, atom_types, xyz, params)


def find_torsion_xyz( 	   atom_types: List[int], 
						   xyz: List[float],
                           params: List[int | float]) -> str:
    """
    Find all distances with current params in a crystal.
    Variables:
	  atom_types: list of atomic types
	  xyz: list of atomic cartesian coordinates - three to one atom
      params: int type1, int type2, int type3, int type4 -- types
              float min12, float max12, float min23, float max23, float min34, float max34, -- distances min/max
              float min123, float max123, float min234, float max234, -- angles min/max (in deg.)
              float min1234, float max1234 -- torsion min/max (in deg.)
        totally 16 parameters should be defined
        angles and torsions are calculated in the range (-180, +180] degrees.
        min/max of distances could be zero - in this case bond's values will be used 
    """
    return cpplib.FindTorsionWC(atom_types, xyz, params)