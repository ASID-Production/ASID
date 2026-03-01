from typing import List, Dict, Tuple


def GenBonds(atoms: List[Tuple[int, float, float, float]]) -> Dict[str, List[Tuple[int, int]]]:
    """
    Identify all bonded atom index pairs in a cluster given Cartesian atom coordinates.
    
    Parameters:
        atoms (List[Tuple[int, float, float, float]]): List of atoms where each entry is
            (AtomType, x, y, z). AtomType is an integer; x, y, z are Cartesian coordinates.
    
    Returns:
        Dict[str, List[Tuple[int, int]]]: Dictionary mapping a string key to a list of
        bonded atom index pairs (AtomIndex1, AtomIndex2). Indices refer to positions in
        the input `atoms` list.
    """
    ...
def GenBondsEx(atoms: List[Tuple[int, float, float, float]]) -> Dict[str, List[(int,int,float)]]:
    """
    Find all bonds in a cluster and report their lengths.
    
    Parameters:
        atoms (List[Tuple[int, float, float, float]]): List of atoms where each atom is a tuple
            (AtomType, x, y, z). AtomType is an integer; x, y, z are Cartesian coordinates as floats.
    
    Returns:
        Dict[str, List[Tuple[int, int, float]]]: A dictionary mapping a string key to a list of
            tuples (AtomIndex1, AtomIndex2, Length). AtomIndex1 and AtomIndex2 are indices into the
            input `atoms` list identifying bonded atom pairs; Length is the distance between them.
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
def Cluster(cell_params: List[float], 
            symms: List[str], 
            atoms: List[Tuple[int, float, float, float]], 
            anchors: List[Tuple[float,float,float,float]]) -> Dict[
    "points":  List[Dict[
        "index":      int,
        "point_frac": Tuple[float,float,float],
        "shift":      Tuple[int,int,int],
        "symmref":    int,
        "type":       int],
    "hasPolymer": bool]]:
    """
            Build a cluster of symmetry-generated atoms around the given anchor positions.
            
            Parameters:
                cell_params (List[float]): Six cell parameters [a, b, c, alpha, beta, gamma] where a,b,c are cell lengths (Å)
                    and alpha,beta,gamma are angles (degrees).
                symms (List[str]): SYMM-codes for symmetry generation; the first entry (equivalent to "x,y,z") is ignored.
                atoms (List[Tuple[int, float, float, float]]): List of atoms as (AtomType, x, y, z) where AtomType is an int
                    and x,y,z are internal (fractional) coordinates.
                anchors (List[Tuple[float, float, float, float]]): List of anchors as (x, y, z, radius) where x,y,z are internal
                    coordinates and radius is generation radius in angstroms.
            
            Returns:
                dict: A dictionary with two keys:
                    "points": List[dict] — each dict describes a generated atom with keys:
                        "index" (int): index of the atom in the original atoms list,
                        "point_frac" (Tuple[float, float, float]): fractional coordinates (fx, fy, fz),
                        "shift" (Tuple[int, int, int]): integer translation shift (dx, dy, dz) relative to the symmetry reference,
                        "symmref" (int): index of the SYMM-code used as reference,
                        "type" (int): atomic type.
                    "hasPolymer" (bool): True if the cell contains a polymer (MOF) structure, False otherwise.
            """
    ...
def VoronoiCalculation(cell_params: List[float], 
                       symms: List[str], 
                       atoms: List[Tuple[int, float, float, float]], 
                       bools: List[bool], 
                       cutoff: float) -> Dict[
    "vertices":  List[Tuple[float,float,float]],
    "edges":  List[List[int,int]],
    "polygons":  List[Dict[
        "vertexes": List[int],
        "edges":  List[int],
        "atoms":  List[int]]],
    "polyhedra": List[Dict[
        "vertexes": List[int],
        "edges":  List[int],
        "polygons":  List[int]]]]:
    """
                       Compute the Voronoi tessellation for a set of atoms given a unit cell and symmetry information.
                       
                       Parameters:
                           cell_params (List[float]): Six unit-cell parameters in order [a, b, c, alpha, beta, gamma]
                               where a, b, c are lengths in angstroms and alpha, beta, gamma are angles in degrees.
                           symms (List[str]): List of SYMM-codes for the structure; the first entry (equivalent to "x,y,z") is ignored.
                           atoms (List[Tuple[int, float, float, float]]): List of atoms as (AtomType, x, y, z) using internal (fractional) coordinates.
                           bools (List[bool]): Flags indicating which atoms should be covered by Voronoi cells; order corresponds to `atoms`.
                           cutoff (float): Distance cutoff (in angstroms) used when building the tessellation.
                       
                       Returns:
                           dict: A dictionary with the following keys:
                               "vertices": List[Tuple[float, float, float]] — coordinates of all Voronoi vertices and cell centers (as 3-tuples).
                               "edges": List[List[int, int]] — list of edges given as pairs of vertex indices.
                               "polygons": List[Dict] — each polygon is a dict with:
                                   "vertexes" (List[int]) — indices of vertices forming the polygon (ordered around the polygon),
                                   "edges" (List[int]) — indices of edges that bound the polygon,
                                   "atoms" (List[int]) — indices of atoms adjacent to the polygon.
                               "polyhedra": List[Dict] — each polyhedron is a dict with:
                                   "vertexes" (List[int]) — indices of vertices belonging to the polyhedron,
                                   "edges" (List[int]) — indices of edges belonging to the polyhedron,
                                   "polygons" (List[int]) — indices of polygons that form the polyhedron.
                       
                       Notes:
                           - The sequence of polyhedra corresponds to the atoms flagged `True` in `bools`, in the same order.
                           - All indices in "edges", "polygons", and "polyhedra" refer to positions in the returned lists.
                       """
    ...
