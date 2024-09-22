from typing import List, Dict, Tuple


def GenBonds(atoms: List[Tuple[int, float, float, float]]) -> Dict["bonds" : List[(int,int)]]:
    """
        Finds all bonds in a cluster (real space).
        Variables:
          atoms: tuple of [list of atomic type and three real-space coordinate].
            AtomType is integer, coordinates are cartesian floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns: List of tuples ( AtomIndex1, AtomIndex2 ), where
          AtomIndex1 and AtomIndex2: indexes in parameter List 'atoms' of atoms,
            which form a bond.
    """
    ...
def GenBondsEx(atoms: List[Tuple[int, float, float, float]]) -> Dict["bonds" : List[(int,int,float)]]:
    """
        Finds all bonds and it's length in a cluster (real space).
        Variables:
          atoms: tuple of [list of atomic type and three real-space coordinate].
            AtomType is integer, coordinates are cartesian floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns: List of Tuples ( AtomIndex1, AtomIndex2, Length ), where
          AtomIndex1 and AtomIndex2: indexes in parameter List 'atoms' of atoms,
            which form a bond.
          Length: distance between atoms with indexes AtomIndex1 and AtomIndex2.
    """
    ...
def SearchMain(graph: str, data: List[str], nprocs: int, exact: bool) -> List[int]:
    """
        Search 'graph' in 'data' with graphs.
        Variables:
          graph: String representation of request molecular graph.
          data: List of string representation of molecular graphs.
          nprocs: the number of threads for multiprocessing.
          exact: boolean flag for exact search (True) or substructure search (False).
        Returns: List of successful IDs.
    """
    ...
def CompareGraph(graph_1: str, graph_2: str, exact: bool) -> bool:
    """
        Compare 'graph_1' and 'graph_2'.
        Variables:
          graph_1: String representation of request molecular graph.
          graph_2: String representation of target molecular graph.
          exact: boolean flag for exact search (True) or substructure search (False)
    """
    ...
def FindMoleculesInCell(cell_params: List[float,float,float,float,float,float], symms: List[str], atoms: List[Tuple[int, float, float, float]]) -> Dict[
    "graph_str": str,
    "error_str": str, 
    "xyz_block": Dict["count": int,
                      "atoms": Dict["x": float,
                                    "y": float,
                                    "z": float,
                                    "init_idx": int],
                      "bonds": List[(int,int)]]]:
    """
        Find all different graphs of molecules in a crystal.
        Variables:
          cell_params: List of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
            (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
          symms: SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
          atoms: Tuple of [list of atomic type and three real-space coordinates].
            AtomType is integer, coordinates are cartesian floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns:
          Dictionary with keys ["graph_str", "error_str", "xyz_block"], where
            "graph_str": String representation of united molecular graph.
            "error_str": String with error message. Empty if no error happened.
            "xyz_block": Dictionary with keys ["count", "atoms", "bonds"], where
              "count": Integer, equal number of equivalent molecules.
              "atoms": Dictionary with keys ["x", "y", "z", "init_idx"], where
                "x", "y", "z": Floating point coordinates of atom.
                "init_idx": Integer, equal index of atom in parameter List 'atoms'.
              "bonds": List of bond Tuples (int,int).
    """
    ...
def FindMoleculesWithoutCell(atoms: List[Tuple[int, float, float, float]]) -> Dict[
    "graph_str": str,
    "error_str": str, 
    "xyz_block": Dict["count": int,
                      "atoms": Dict["x": float,
                                    "y": float,
                                    "z": float,
                                    "init_idx": int],
                      "bonds": List[(int,int)]]]:
    """
        Find all different graphs of molecules in a cluster.
        Variables:
          atoms: tuple of [list of atomic type and three real-space coordinates].
            AtomType is integer, coordinates are cartesian floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns:
          Dictionary with keys ["graph_str", "error_str", "xyz_block"], where
            "graph_str": String representation of united molecular graph.
            "error_str": String with error message. Empty if no error happened.
            "xyz_block": Dictionary with keys ["count", "atoms", "bonds"], where
              "count": Integer, equal number of equivalent molecules.
              "atoms": Dictionary with keys ["x", "y", "z", "init_idx"], where
                "x", "y", "z": Floating point coordinates of atom.
                "init_idx": Integer, equal index of atom in parameter List 'atoms'.
              "bonds": List of bond Tuples (int,int).
    """
    ...
def GenSymm(atoms: List[Tuple[int, float, float, float]], flags: int, symms: List[str]) -> List[Tuple[int, float, float, float]]:
    """
        Add symmetry-generated atoms to parameter 'atoms' (internal coordinates).
        Variables:
          atoms: Tuple of [list of atomic type and three real-space coordinates].
            AtomType is integer, coordinates are cartesian floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          flags: Integer of two independent boolean flags [0,3]:
            0x01: Move to cell all generated atoms. ATTENTION! Did not change coordinates
              already existed in 'atoms' parameter.
            0x02: First of all move center of mass of existed atoms.
          symms: SYMM-codes for generation.
        Returns: Modified parameter 'atoms'.
    """
    ...
def FindDistanceIC(cell_params: List[float,float,float,float,float,float], symms: List[str], atoms: List[Tuple[int, float, float, float]], params: List[int | float]) -> Dict["distances": List[Tuple[int,int,float]]]:
    """
        Find all distances with current params in unit cell.
        Variables:
          cell_params: List of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
            (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
          symms:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
          atoms: List of Tuples of atomic type and three internal coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List with exactly 4 parameters:
            type1, type2  -- Atomic types : Integers
            min12, max12  -- Minimal and maximal values of distances : Floating point numbers
            min/max values could be zero - in this case bond's values will be used
        Returns:
          Dictionary with key ["distances"], where
            "distances": List of Tuples (Index1, Index2, Value), where
              Index1 and Index2: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: distance between atoms with indexes Index1 and Index2.
    """
    ...
def FindDistanceWC(atoms: List[Tuple[int, float, float, float]], params: List[int | float]) -> Dict["distances": List[Tuple[int,int,float]]]:
    """
        Find all distances with current params in cartesian cluster.
        Variables:
          atoms: List of Tuples of atomic type and three cartesian coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List with exactly 4 parameters:
            type1, type2  -- Atomic types : Integers
            min12, max12  -- Minimal and maximal values of distances : Floating point numbers
            min/max values could be zero - in this case bond's values will be used
        Returns:
          Dictionary with key ["distances"], where
            "distances": List of Tuples (Index1, Index2, Value), where
              Index1 and Index2: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: distance between atoms with indexes Index1 and Index2.
    """
    ...
def FindAngleIC(cell_params: List[float,float,float,float,float,float], symms: List[str], atoms: List[Tuple[int, float, float, float]], params: List[int | float]) -> Dict["angles": List[Tuple[int,int,int,float]]]:
    """
        Find all angles with current params in unit cell.
        Variables:
          cell_params: List of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
            (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
          symms:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
          atoms: List of Tuples of atomic type and three internal coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List of exactly 9 parameters:
            type1, type2, type3         -- Atomic types : Integers
            min12, max12, min23, max23  -- Minimal and maximal values of distances : Floating point numbers
            min123, max123              -- Minimal and maximal values of angles (in deg.) : Floating point numbers
            Min/Max of distances could be zero - in this case bond's values will be used.
            Angles will be compared in the range (-180, +180] degrees.
        Returns:
          Dictionary with key ["angles"], where
            "angles": List of Tuples (Index1, Index2, Index3, Value), where
              Index1, Index2 and Index3: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: angle in degrees of [Index1, Index2, Index3].
    """
    ...
def FindAngleWC(atoms: List[Tuple[int, float, float, float]], params: List[int | float]) -> Dict["angles": List[Tuple[int,int,int,float]]]:
    """
        Find all angles with current params in cartesian cluster.
        Variables:
          atoms: List of Tuples of atomic type and three cartesian coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List of exactly 9 parameters:
            type1, type2, type3         -- Atomic types : Integers
            min12, max12, min23, max23  -- Minimal and maximal values of distances : Floating point numbers
            min123, max123              -- Minimal and maximal values of angles (in deg.) : Floating point numbers
            Min/Max of distances could be zero - in this case bond's values will be used.
            Angles will be compared in the range (-180, +180] degrees.
        Returns:
          Dictionary with key ["angles"], where
            "angles": List of Tuples (Index1, Index2, Index3, Value), where
              Index1, Index2 and Index3: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: angle in degrees of [Index1, Index2, Index3].
    """
    ...
def FindTorsionIC(cell_params: List[float,float,float,float,float,float], symms: List[str], atoms: List[(int, float, float, float)], params: List[int | float]) -> Dict["tors": List[Tuple[int,int,int,int,float]]]:
    """
        Find all torsion angles with current params in unit cell.
        Variables:
          cell_params: List of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
            (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
          symms:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
          atoms: List of Tuples of atomic type and three internal coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List of exactly 16 parameters:
            type1, type2, type3, type4                -- Atomic types : Integers
            min12, max12, min23, max23, min34, max34  -- Minimal and maximal values of distances : Floating point numbers
            min123, max123, min234, max234            -- Minimal and maximal values of angles (in deg.) : Floating point numbers
            min1234, max1234                          -- Minimal and maximal values of torsion angles (in deg.) : Floating point numbers
            Min/Max of distances could be zero - in this case bond's values will be used.
            Angles and torsion angles will be compared in the range (-180, +180] degrees.
        Returns:
          Dictionary with key ["tors"], where
            "tors": List of Tuples (Index1, Index2, Index3, Index4, Value), where
              Index1, Index2, Index3 and Index4: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: torsion angle in degrees of [Index1, Index2, Index3, Index4].
    """
    ...
def FindTorsionWC(atoms: List[Tuple[int, float, float, float]], params: List[int | float]) -> Dict["tors": List[Tuple[int,int,int,int,float]]]:
    """
        Find all angles with current params in cartesian cluster.
        Variables:
          atoms: List of Tuples of atomic type and three cartesian coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          params: List of exactly 16 parameters:
            type1, type2, type3, type4                -- Atomic types : Integers
            min12, max12, min23, max23, min34, max34  -- Minimal and maximal values of distances : Floating point numbers
            min123, max123, min234, max234            -- Minimal and maximal values of angles (in deg.) : Floating point numbers
            min1234, max1234                          -- Minimal and maximal values of torsion angles (in deg.) : Floating point numbers
            Min/Max of distances could be zero - in this case bond's values will be used.
            Angles and torsion angles will be compared in the range (-180, +180] degrees.
        Returns:
          Dictionary with key ["tors"], where
            "tors": List of Tuples (Index1, Index2, Index3, Index4, Value), where
              Index1, Index2, Index3 and Index4: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: torsion angle in degrees of [Index1, Index2, Index3, Index4].
    """
    ...
def FindDAT_IC(cell_params: List[float,float,float,float,float,float], symms: List[str], atoms: List[(int, float, float, float)]) -> Dict[
    "bonds":  List[Tuple[int,int,float]],
    "angles": List[Tuple[int,int,int,float]],
    "tors":   List[Tuple[int,int,int,int,float]]]:
    """
        Create dictionary with distances, angles and torsions in cell.
        Variables:
          cell_params: List of exactly 6 cell parameters in strict order: [a, b, c, alpha, beta, gamma].
            (a, b, c - are translation vectors (in Angstroms) and alpha, beta, gamma - are angles (in degrees)
          symms:  SYMM-codes of structure. Should contain 'x,y,z' (equivalent) as first ([0]) symmetry - it is ignored.
          atoms: List of Tuples of atomic type and three internal coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns:
          Dictionary with keys ["bonds","angles","tors"], where
            "bonds": List of Tuples (Index1, Index2, Value), where
              Index1 and Index2: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: distance between atoms with indexes Index1 and Index2.
            "angles": List of Tuples (Index1, Index2, Index3, Value), where
              Index1, Index2 and Index3: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: angle in degrees of [Index1, Index2, Index3].
            "tors": List of Tuples (Index1, Index2, Index3, Index4, Value), where
              Index1, Index2, Index3 and Index4: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: torsion angle in degrees of [Index1, Index2, Index3, Index4].
    """
    ...
def FindDAT_WC(atoms: List[(int, float, float, float)]) -> Dict[
    "bonds":  List[Tuple[int,int,float]],
    "angles": List[Tuple[int,int,int,float]],
    "tors":   List[Tuple[int,int,int,int,float]]]:
    """
        Create dictionary with distances, angles and torsions in xyz.
        Variables:
          atoms: List of Tuples of atomic type and three cartesian coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
        Returns:
          Dictionary with keys ["bonds","angles","tors"], where
            "bonds": List of Tuples (Index1, Index2, Value), where
              Index1 and Index2: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: distance between atoms with indexes Index1 and Index2.
            "angles": List of Tuples (Index1, Index2, Index3, Value), where
              Index1, Index2 and Index3: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: angle in degrees of [Index1, Index2, Index3].
            "tors": List of Tuples (Index1, Index2, Index3, Index4, Value), where
              Index1, Index2, Index3 and Index4: Integer indexes in parameter List 'atoms' of atoms,
                which are on requested distance.
              Value: torsion angle in degrees of [Index1, Index2, Index3, Index4].
    """
    ...
def himp(atoms: List[Tuple[int, float, float, float]], value: float | List[float]) -> List[Tuple[int, float, float, float]]:
    """
        Moves hydrogens to the nearest atom.
        Variables:
          atoms: List of Tuples of atomic type and three cartesian coordinates.
            AtomType is integer, coordinates are floating point numbers.
            example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
          value: Could be a single value or list of values. 
            single value: move all hydrogens on the 'value' distance to the nearest atom
            List of values: move all hydrogens on the 'value[type]' distance to the nearest 
              atom, depending on atom type. 
            Attention! value[0] - ignored, but should not be empty! So, value[6] - is a C-H distance.
        Returns:
          Dictionary with keys ["atoms","error_str"], where
            "atoms": List of Tuples of atomic type and three cartesian coordinates.
              AtomType is integer, coordinates are floating point numbers.
              example [ [ 1, 0.0, 0.0, 0.0 ], [9, 0.5, 0.5, 0.5], ... ]
            "error_str": optional key. String with error message. Exists only if something gone wrong.
    """
    ...