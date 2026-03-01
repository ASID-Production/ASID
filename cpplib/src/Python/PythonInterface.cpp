// Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
// This file is part of ASID - Atomistic Simulation Instruments and Database
// For more information see <https://github.com/ASID-Production/ASID>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// ******************************************************************************************
//  Author:      Alexander A. Korlyukov (head)
//  ORCID:       0000-0002-5600-9886
//  Author:      Alexander D. Volodin (author of cpplib)
//  ORCID:       0000-0002-3522-9193
//  Author:      Petr A. Buikin (author of api_database)
//  ORCID:       0000-0001-9243-9915
//  Author:      Alexander R. Romanenko (author of VnE)
//  ORCID:       0009-0003-5298-6836
//
// ******************************************************************************************
#include <array>
#include <cmath>
#include <list>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../Functions/Functions.h"
#include "../BaseHeaders/BaseTypes.h"
#include "../BaseHeaders/DebugMes.h"
#include "../Classes/Distances.h"
#include "../Classes/FindMolecules.h"
#include "../Classes/Geometry.h"
#include "../Classes/Voronoi.h"
#include "../Python/PythonInterface.h"

extern const cpplib::Distances* p_distances;

using namespace cpplib;
using namespace cpplib::currents;
using namespace cpplib::basic_types;
using PointType = cpplib::FAM_Cell::PointType;

enum class ErrorState {
	OK = 0,
	UnknownError,
};
struct Prepare_WC {
	std::vector<cpplib::basic_types::AtomTypeBase> types;
	std::vector<cpplib::geometry::Point<FloatingPointType>> points;
	explicit Prepare_WC(PyObject* otuples) {
		Py_ssize_t s = PyList_Size(otuples);
		types.reserve(static_cast<size_t>(s));
		points.reserve(static_cast<size_t>(s));

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* o_tuple = PyList_GetItem(otuples, i);
			types.push_back(static_cast<AtomTypeBase>(PyLong_AsLong(PyTuple_GetItem(o_tuple, 0))));
			points.emplace_back(static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 1))),
								static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 2))),
								static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 3))));
		}
	}
};
struct Prepare_IC : public Prepare_WC {
	std::array<cpplib::basic_types::FloatingPointType, 6> cell;
	std::vector<const char*> symm;
	/**
	 * @brief Construct a Prepare_IC by reading unit-cell parameters and symmetry labels from Python objects.
	 *
	 * Initializes the six-element cell array from a Python sequence of six numeric values and fills the symmetry
	 * label vector from a Python list of Unicode strings. The atomic tuples are parsed by the base Prepare_WC constructor.
	 *
	 * @param ocell Python sequence (length 6) of numeric cell parameters in the order expected by the code.
	 * @param osymm Python list of symmetry label strings; each entry is converted to a C string stored in `symm`.
	 * @param otuples Python iterable of atom tuples passed to the Prepare_WC base constructor for type/coordinate parsing.
	 */
	Prepare_IC(PyObject* ocell, PyObject* osymm, PyObject* otuples) : Prepare_WC(otuples) {
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}
		const Py_ssize_t s = PyList_Size(osymm);
		symm.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			Py_ssize_t us;
			symm[i] = PyUnicode_AsUTF8AndSize(PyList_GetItem(osymm, i), &us);
		}
	}
};

extern "C" {
	inline static void useDistances(PyObject* self);
}

template <char times>
/**
 * @brief Parse a sequence of distance parameter pairs from a Python list and resolve defaults.
 *
 * Reads 2*`times` numeric entries from the Python list `oparams` beginning at index `d` and returns them as an array of (min, max) distance pairs for atom-type pairs defined by `type`. If an extracted value equals 0, the function initializes (if necessary) the global distance data and replaces that value with the corresponding minimum or maximum distance for the atom-type pair.
 *
 * @tparam times Number of distance pairs to parse.
 * @param self Python module/object used to initialize distance data when needed.
 * @param oparams Python list containing numeric parameters; must contain at least `d + 2*times` entries.
 * @param type Array of atom type indices of size `times + 1`; pair i is (type[i], type[i+1]).
 * @param d Index in `oparams` where parsing starts; advanced by 2 for each parsed pair (overall increment by 2*`times`).
 * @return std::array<std::pair<FloatingPointType, FloatingPointType>, times> Array of (min, max) distance pairs (values in project FloatingPointType units).
 *
 * @note The function may call `useDistances(self)` and access the global `p_distances` to obtain default distances when an input value is 0.
 */
static std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, times> FindDParamsParse(PyObject* self, PyObject* oparams, const std::array<int, times + 1> type, char& d) {
	std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, times> value;
	for (char i = 0; i < times; i++, d += 2) {
		value[i].first = static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(oparams, d)));
		value[i].second = static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(oparams, d + 1)));

		if (value[i].first == 0) {
			useDistances(self);
			value[i].first = p_distances->minDistance(static_cast<cpplib::basic_types::AtomTypeBase>(type[i]), static_cast<cpplib::basic_types::AtomTypeBase>(type[i + 1]));
		}

		if (value[i].second == 0) {
			useDistances(self);
			value[i].second = p_distances->maxDistance(static_cast<cpplib::basic_types::AtomTypeBase>(type[i]), static_cast<cpplib::basic_types::AtomTypeBase>(type[i + 1]));
		}
	}
	return value;
}
template <char times>
static std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, times> FindATParamsParse(PyObject* oparams, char& d) {
	std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, times> value;
	for (char i = 0; i < times; i++, d += 2) 
{
		value[i].first = static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(oparams, d)));
		value[i].second = static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(oparams, d + 1)));
	}
	return value;
}
extern "C" {

	inline static ErrorState pyListToVectorFloat(PyObject* plist, std::vector<cpplib::basic_types::FloatingPointType>* pret) {
		const Py_ssize_t s = PyList_Size(plist);
		std::vector<cpplib::basic_types::FloatingPointType>& ret = *pret;
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			ret[i] = static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i)));
		}
		return ErrorState::OK;
	}
	inline static ErrorState pyListToVectorCurPoint(PyObject* plist, std::vector<PointType>* retp) {
		std::vector<PointType>& ret = *retp;
		const Py_ssize_t s = PyList_Size(plist);
		ret.reserve(s);
		for (Py_ssize_t i = 0; i < s; i += 3) {
			ret.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i))),
							 static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i + 1))),
							 static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i + 2))));
		}
		return ErrorState::OK;
	}
	inline static ErrorState pyListToVectorInt(PyObject* plist, std::vector<int>* pret) {
		const Py_ssize_t s = PyList_Size(plist);
		std::vector<int>& ret = *pret;
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			ret[i] = static_cast<int>(PyLong_AsLong(PyList_GetItem(plist, i)));
		}
		return ErrorState::OK;
	}
	/**
	 * @brief Convert a Python list of strings to a C-style string pointer vector.
	 *
	 * Populates the provided vector with UTF-8 C-string pointers obtained from
	 * each element of the Python list and resizes the vector to the list length.
	 *
	 * @param plist Python list whose elements are expected to be Unicode strings.
	 * @param pret Pointer to the vector that will be filled with `const char*`
	 *             pointers to UTF-8 encoded string data.
	 *             The pointers refer to memory owned by the corresponding Python
	 *             string objects and remain valid only while those Python objects
	 *             are alive and unchanged.
	 * @return ErrorState::OK on success.
	 */
	inline static ErrorState pyListToVectorCharP(PyObject* plist, std::vector<const char*>* pret) {

		std::vector<const char*>& ret = *pret;
		const Py_ssize_t s = PyList_Size(plist);
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			Py_ssize_t us;
			ret[i] = PyUnicode_AsUTF8AndSize(PyList_GetItem(plist, i), &us);
		}
		return ErrorState::OK;
	}
	/**
	 * @brief Parses a Python list of (type, x, y, z) tuples into parallel vectors of atom types and points.
	 *
	 * Clears and fills `types` and `points` so that each entry corresponds to one tuple from `o_list`.
	 *
	 * @param o_list Python list where each element is a 4-tuple: (integer atom type, float x, float y, float z).
	 * @param types Output vector populated with atom types cast to AtomTypeBase.
	 * @param points Output vector populated with 3D points (PointType) from the tuple coordinates.
	 * @return ErrorState `OK` on successful parse.
	 */
	inline static ErrorState pyTXYZparse(PyObject* o_list, std::vector<AtomTypeBase>* types, std::vector<PointType>* points) {
		Py_ssize_t s = PyList_Size(o_list);
		types->clear();
		types->reserve(static_cast<int>(s));
		points->clear();
		points->reserve(static_cast<int>(s));

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* o_tuple = PyList_GetItem(o_list, i);
			types->push_back(static_cast<AtomTypeBase>(PyLong_AsLong(PyTuple_GetItem(o_tuple, 0))));
			points->emplace_back(static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 1))),
								 static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 2))),
								 static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 3))));
		}
		return ErrorState::OK;
	}
	/**
	 * @brief Ensure the global distances instance is initialized.
	 *
	 * Initializes a single static Distances object from the "BondLength.ini" file
	 * located in the same directory as the provided module object and stores its
	 * address in the global pointer `p_distances`. If `p_distances` is already
	 * non-null, the function returns without performing any work.
	 *
	 * @param self Python module or module-like object whose `__file__` attribute is used
	 *             to locate the directory containing `BondLength.ini`.
	 */
	inline static void useDistances(PyObject* self) {
		if (p_distances != nullptr) {
			return;
		}
		Py_ssize_t us; 
		PyObject* file_obj = PyObject_GetAttrString(self, "__file__");
		std::string full(PyUnicode_AsUTF8AndSize(file_obj, &us));
		Py_DECREF(file_obj);

		auto found = full.find_last_of("\\/");
		auto bond_filename = full.substr(0, found + 1) + "BondLength.ini";

		static Distances dist(bond_filename);
		p_distances = &dist;
	}

	/**
	 * @brief Generate a list of bonds for the provided atom tuples in Cartesian coordinates.
	 *
	 * Parses the Python sequence in @p arg as a list of (type, x, y, z) tuples, computes bonds using
	 * the library's distance rules, and returns a Python dictionary containing the bond list.
	 *
	 * @param arg Python sequence of atom tuples (type, x, y, z).
	 * @return PyObject* Python dictionary { "bonds": <list> } where <list> is a sequence of bond records. */
	static PyObject* cpplib_GenBonds(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *p_distances;
		Prepare_WC all(arg);

		FAM_Struct famstr(std::move(all.types), std::move(all.points));
		std::string errM;
		auto&& bonds = famstr.findBonds(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		PyObject* lst = py_util::convert(bonds);
		return Py_BuildValue("{s:N}",
							 "bonds", lst);
	}
	/**
	 * @brief Generate bonds (including bond lengths) for a set of atoms provided as tuples.
	 *
	 * Parses the provided Python sequence of atom tuples and computes the bonding list with distances.
	 *
	 * @param arg Python iterable of 4-element tuples (atom_type, x, y, z) describing atoms.
	 * @return PyObject* A Python dict with a single key "bonds" whose value is a list of tuples (i, j, length)
	 * where `i` and `j` are atom indices and `length` is the bond length as a floating-point value.
	 */
	static PyObject* cpplib_GenBondsEx(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *p_distances;
		Prepare_WC all(arg);
		FAM_Struct famstr(std::move(all.types), std::move(all.points));
		std::string errM;
		// TODO: Don't use famstr anymore
		auto&& bonds = famstr.findBondsEx(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		PyObject* lst = py_util::convert(bonds);
		return Py_BuildValue("{s:N}",
							 "bonds", lst);
	}

	/**
	 * @brief Parse Python arguments and perform a substructure search via SearchMain.
	 *
	 * Parses a Python argument tuple expected as (search: str, data: list[str], np: int, exact: bool),
	 * converts the data list to a vector of C strings, calls SearchMain(search, data, np, exact),
	 * and returns the converted Python representation of SearchMain's result.
	 *
	 * @param self Unused module/self pointer.
	 * @param args Python argument tuple: (search, data, np, exact).
	 * @return PyObject* The Python-converted result of SearchMain. Returns Py_None if argument parsing fails.
	 */
	static PyObject* cpplib_SearchMain(PyObject* self, PyObject* args) {
		const char* search = NULL;
		PyObject* o = NULL;
		int np = 0;
		int exact = 0;
		if (!PyArg_ParseTuple(args, "sOip", &search, &o, &np, &exact)) {
			Py_RETURN_NONE;
		}
		std::vector< const char*> data; pyListToVectorCharP(o, &data);
		const auto ret = SearchMain(search, std::move(data), np, (exact != 0));
		return py_util::convert(ret);
	}
	/**
	 * @brief Compare two graph representations using an exactness flag.
	 *
	 * Parses a Python tuple (s1, s2, exact) where s1 and s2 are graph strings
	 * and `exact` is a boolean; returns whether the two graphs match.
	 *
	 * If argument parsing fails, the function returns None.
	 *
	 * @param self Unused Python module or object pointer.
	 * @param args A Python tuple (s1, s2, exact) to be parsed.
	 * @returns PyObject* `True` if the graphs match under the given exactness flag, `False` otherwise; or `None` if parsing fails.
	 */
	static PyObject* cpplib_CompareGraph(PyObject* self, PyObject* args) {
		const char* s1 = NULL;
		const char* s2 = NULL;
		int b = 0;
		if (!PyArg_ParseTuple(args, "ssp", &s1, &s2, &b)) {
			Py_RETURN_NONE;
		}
		if (CompareGraph(s1, s2, (b != 0))) {
			Py_RETURN_TRUE;
		}
		else {
			Py_RETURN_FALSE;
		}
	}

	/**
	 * @brief Find molecules inside a unit cell and return their graph and XYZ representation.
	 *
	 * Parses a unit cell, symmetry, and atom tuples from Python arguments, searches for molecules
	 * within the cell, and returns a Python dictionary describing the result.
	 *
	 * @return PyObject* A Python dictionary with keys:
	 *   - "graph_str": C string with the found molecular graph representation.
	 *   - "error_str": C string with any error message (empty on success).
	 *   - "xyz_block": list of molecule entries; each entry is a dict with:
	 *       - "count": integer molecule identifier.
	 *       - "atoms": list of atom dicts with keys "x", "y", "z" (float coordinates) and "init_idx" (long).
	 *       - "bonds": list of 2-tuples (int, int) describing bonds between atom indices.
	 *
	 * Returns Python None if argument parsing fails.
	 */
	static PyObject* cpplib_FindMoleculesInCell(PyObject* self, PyObject* args) {
		// TODO: Repair memory leaks!
		useDistances(self);
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuple = NULL;
		if (!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuple)) {
			Py_RETURN_NONE;
		}

		Prepare_IC all(ocell, osymm, otuple);

		auto [graph, error, ret] = FindMoleculesInCell(all.cell, all.symm, all.types, all.points);
		PyObject* o_xyz_block = PyList_New(0);

		for (const auto& [atoms, id, bonds] : ret) {
			PyObject* o_molecule = PyList_New(0);
			for (auto& [point, type] : atoms) {
				PyObject* o_atom = Py_BuildValue("{s:f,s:f,s:f,s:l}",
												 "x", float(point[0]),
												 "y", float(point[1]),
												 "z", float(point[2]),
												 "init_idx", long(id));
				PyList_Append(o_molecule, o_atom);
			}
			PyObject* o_bonds = PyList_New(0);
			for (const auto& bond : bonds) {
				PyObject* o_bond1 = Py_BuildValue("(ii)", int(bond.first), int(bond.second));
				PyList_Append(o_bonds, o_bond1);
			}

			PyList_Append(o_xyz_block, Py_BuildValue("{s:l,s:O,s:O}",
													 "count", long(id),
													 "atoms", o_molecule,
													 "bonds", o_bonds));
		}
		return Py_BuildValue("{s:s,s:s,s:O}",
							 "graph_str", graph.c_str(),
							 "error_str", error.c_str(),
							 "xyz_block", o_xyz_block);
	}
	/**
	 * @brief Find molecules from a list of atomic tuples without unit cell information.
	 *
	 * Parses the provided Python sequence of atomic tuples, groups atoms into distinct
	 * molecules, and returns a Python dictionary containing a graph representation,
	 * any error message, and a detailed block describing found molecules.
	 *
	 * @param otuple Python iterable of atom tuples. Each tuple must contain an atom
	 * type (first element) followed by three numeric coordinates (x, y, z).
	 *
	 * @return PyObject* A Python dict with keys:
	 * - "graph_str": C-string serialization of the molecule graph.
	 * - "error_str": C-string containing an error message or empty string on success.
	 * - "xyz_block": list of molecule dictionaries; each molecule dict contains:
	 *   - "count": long number of atoms in the molecule.
	 *   - "atoms": list of atom dictionaries, each with keys:
	 *     - "x", "y", "z": floating-point Cartesian coordinates.
	 *     - "init_idx": original atom index (long).
	 *   - "bonds": list of 2-tuples (int, int) describing bonded atom index pairs.
	 */
	static PyObject* cpplib_FindMoleculesWithoutCell(PyObject* self, PyObject* otuple) {
		// TODO: Repair memory leaks!
		useDistances(self);

		Prepare_WC all(otuple);

		auto ret = FindMoleculesWithoutCell(all.types, all.points);
		PyObject* o_xyz_block = PyList_New(0);

		for (auto& mol : std::get<2>(ret)) {
			PyObject* o_molecule = PyList_New(0);
			for (auto& atom : std::get<0>(mol)) {
				PyObject* o_atom = Py_BuildValue("{s:f,s:f,s:f,s:l}",
												 "x", cpplib::basic_types::FloatingPointType(std::get<0>(atom)[0]),
												 "y", cpplib::basic_types::FloatingPointType(std::get<0>(atom)[1]),
												 "z", cpplib::basic_types::FloatingPointType(std::get<0>(atom)[2]),
												 "init_idx", long(std::get<1>(atom)));
				PyList_Append(o_molecule, o_atom);
			}
			PyObject* o_bonds = PyList_New(0);
			for (auto& bond : std::get<2>(mol)) {
				PyObject* o_bond1 = Py_BuildValue("(ii)", int(bond.first), int(bond.second));
				PyList_Append(o_bonds, o_bond1);
			}

			PyList_Append(o_xyz_block, Py_BuildValue("{s:l,s:O,s:O}",
													 "count", long(std::get<1>(mol)),
													 "atoms", o_molecule,
													 "bonds", o_bonds));
		}
		return Py_BuildValue("{s:s,s:s,s:O}",
							 "graph_str", std::get<0>(ret).c_str(),
							 "error_str", std::get<1>(ret).c_str(),
							 "xyz_block", o_xyz_block);
	}

	/**
	 * @brief Generate symmetry-equivalent atom positions and append them to a Python list.
	 *
	 * Parses a Python argument tuple (otuples, flags, osymm), generates atoms produced by
	 * applying the provided symmetry operations, and appends new atom tuples to `otuples`.
	 *
	 * @param args Python argument tuple expected to contain:
	 *   - otuples: Python list of atom tuples (type, x, y, z) to which new atoms will be appended.
	 *   - flags: byte-valued flags where bit 0 enables moving generated atoms into the unit cell
	 *     and bit 1 enables centering symmetry by moving the center of mass into the cell.
	 *   - osymm: Python list of symmetry operation definitions (strings) to apply.
	 * @return The same Python list object passed as `otuples` with generated symmetry atom tuples appended;
	 *         `Py_None` if argument parsing fails.
	 */
	static PyObject* cpplib_GenSymm(PyObject* self, PyObject* args) {
		// TODO: Repair memory leaks!
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		std::byte flags;

		if (!PyArg_ParseTuple(args, "OBO", &otuples, &flags, &osymm)) {
			Py_RETURN_NONE;
		}
		bool movetocell = (flags & std::byte(1)) != std::byte(0);
		bool movemasstocell = (flags & std::byte(2)) != std::byte(0);

		Prepare_WC all(otuples);
		std::vector<const char*> nsymm; pyListToVectorCharP(osymm, &nsymm);
		std::vector<geometry::Symm<FloatingPointType>> symm;
		const size_t ss = nsymm.size();
		for (size_t i = 0; i < ss; i++) {
			symm.emplace_back(nsymm[i]);
		}

		// Move center of mass to cell
		const auto s_points = all.points.size();
		const auto sf_points = static_cast<FloatingPointType>(s_points);
		if (movemasstocell) {
			PointType centerofmass(0, 0, 0);
			for (size_t i = 0; i < s_points; i++) {
				centerofmass += all.points[i];
			}
			centerofmass /= sf_points;
			PointType ceilmass(std::ceil(centerofmass[0]),
							   std::ceil(centerofmass[1]),
							   std::ceil(centerofmass[2]));

			for (size_t i = 0; i < ss; i++) {
				PointType movedcenter = symm[i].GenSymm(centerofmass);
				PointType ceilmoved(std::ceil(movedcenter[0]),
									std::ceil(movedcenter[1]),
									std::ceil(movedcenter[2]));
				symm[i].point += ceilmass - ceilmoved;
			}
		}

		FAM_Struct famstr(std::move(all.types), std::move(all.points));
		FAM_Cell fcell(FAM_Cell::base(32, 32, 32, 90, 90, 90, true));
		WITH_LOG_M(fcell, GenerateSymm, famstr, symm, movetocell, true);

		for (Py_ssize_t i = s_points; i < famstr.sizePoints; i++) {
			PyList_Append(otuples, Py_BuildValue("(Ifff)",
												 static_cast<unsigned int>(famstr.types[std::get<0>(famstr.parseIndex[i])]),
												 static_cast<cpplib::basic_types::FloatingPointType>(famstr.points[i][0]),
												 static_cast<cpplib::basic_types::FloatingPointType>(famstr.points[i][1]),
												 static_cast<cpplib::basic_types::FloatingPointType>(famstr.points[i][2])));
		}
		Py_INCREF(otuples);
		return otuples;
	}

	/**
	 * @brief Compute distances between atom pairs of a specified type pair within the given cell and symmetry.
	 *
	 * Expects a Python tuple of four objects: (ocell, osymm, otuples, oparams).
	 *
	 * @param self Unused Python module/self pointer.
	 * @param args A 4-item tuple:
	 *   - ocell: sequence of six floating-point cell parameters.
	 *   - osymm: sequence of symmetry descriptor strings.
	 *   - otuples: sequence of atom tuples (type, x, y, z) in Cartesian coordinates.
	 *   - oparams: sequence where the first two elements are integer atom-type indices specifying the pair to analyze;
	 *              additional elements provide optional distance parameter controls.
	 *
	 * @return A Python dictionary with key "distances" whose value is a list of 3-tuples (i, j, d):
	 *         - i: first atom index (unsigned int)
	 *         - j: second atom index (unsigned int)
	 *         - d: distance between atoms i and j as a float
	 *
	 * @note Returns None immediately if argument parsing fails.
	 */
	static PyObject* cpplib_FindDistanceIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;

		if (!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			Py_RETURN_NONE;
		}
		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 2> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		char d = 2;
		auto value = FindDParamsParse<1>(self, oparams, type, d);

		auto res = WITH_LOG(FindDistanceIC, all.cell, all.symm, all.types, all.points, type, value[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<float>(std::get<2>(res[i]))));
		}
		return Py_BuildValue("{s:O}", "distances", lst);
	}
	/**
	 * @brief Extracts distances between atom pairs of specified types from world-coordinate tuples.
	 *
	 * Parses Python arguments (otuples, oparams), interprets the first two entries of oparams as
	 * atom-type indices to search for, computes matching pairwise distances from otuples, and
	 * returns them as a Python dictionary under the "distances" key.
	 *
	 * @param otuples Python list of tuples representing atoms; each tuple is expected to contain an
	 *                atom type followed by three coordinates (x, y, z).
	 * @param oparams Python list where the first two elements are integer atom-type indices. If
	 *                distance bounds in additional parameters are zero, default bond-length limits
	 *                may be used to resolve them.
	 * @return PyObject* A Python dict { "distances": list } where each item in the list is a
	 *         tuple (I, I, F): the two atom indices (unsigned int) and the distance between them
	 *         as a float.
	 */
	static PyObject* cpplib_FindDistanceWC(PyObject* self, PyObject* args) {
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			Py_RETURN_NONE;
		}
		Prepare_WC all(otuples);

		const std::array<int, 2> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		char d = 2;
		auto value = FindDParamsParse<1>(self, oparams, type, d);


		auto res = FindDistanceWC(all.types, all.points, type, value[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<float>(std::get<2>(res[i]))));
		}
		return Py_BuildValue("{s:O}", "distances", lst);
	}

	/**
	 * @brief Compute angles for a specified triplet of atom types using cell and symmetry information.
	 *
	 * Parses (ocell, osymm, otuples, oparams) from the Python argument tuple, interprets the first
	 * three entries of `oparams` as the atom-type indices that define the angle (central atom is the second),
	 * applies distance and angular parameter filters, and returns all matching angles.
	 *
	 * @param ocell Python sequence representing the unit cell (six floating-point values).
	 * @param osymm Python sequence of symmetry labels/operations corresponding to `otuples`.
	 * @param otuples Python sequence of atomic tuples (type and coordinates) for the input structure.
	 * @param oparams Python list of parameters where the first three elements are integer type indices
	 *                (i, j, k) and the remaining values supply distance/angle parameter pairs.
	 * @return PyObject* A Python dict with key "angles" whose value is a list of tuples `(I, I, I, F)`:
	 *         the three integer atom indices followed by the angle in degrees as a float.
	 */
	static PyObject* cpplib_FindAngleIC(PyObject* self, PyObject* args) {
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			Py_RETURN_NONE;
		}
		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 3> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		auto value_d = FindDParamsParse<2>(self, oparams, type, d);
		auto value_a = FindATParamsParse<1>(oparams, d);


		auto res = FindAngleIC(all.cell, all.symm, all.types, all.points, type, value_d, value_a[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(res[i])))));
		}
		return Py_BuildValue("{s:O}", "angles", lst);
	}
	/**
	 * @brief Finds angles (three-atom motifs) in a world-coordinate atom list matching a specified type pattern.
	 *
	 * Expects Python arguments (otuples, oparams) where:
	 * - otuples is a sequence of atom tuples (type, x, y, z) in world coordinates.
	 * - oparams is a list whose first three entries are integer atom-type indices (A, B, C) followed by distance/angle parameter values;
	 *   zero values for distance bounds trigger lookup from the global distances table.
	 *
	 * @return A Python dictionary with key "angles" whose value is a list of 4-tuples (I, I, I, F):
	 *         the three integer atom indices (A, B, C) followed by the measured angle in degrees as a float.
	 */
	static PyObject* cpplib_FindAngleWC(PyObject* self, PyObject* args) {
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			Py_RETURN_NONE;
		}

		Prepare_WC all(otuples);

		const std::array<int, 3> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
									   static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
									   static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		auto value_d = FindDParamsParse<2>(self, oparams, type, d);
		auto value_a = FindATParamsParse<1>(oparams, d);


		auto res = FindAngleWC(all.types, all.points, type, value_d, value_a[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(res[i])))));
		}
		return Py_BuildValue("{s:O}", "angles", lst);
	}
	/**
	 * @brief Find torsion angles for specified atom type indices within an input cell.
	 *
	 * Parses four type indices and parameter tuples from the provided Python arguments,
	 * runs the torsion search on the prepared unit cell and atomic coordinates, and
	 * returns the found torsions as Python structures.
	 *
	 * @param self Unused Python module/state pointer.
	 * @param args A Python tuple containing (ocell, osymm, otuples, oparams):
	 *             - ocell: cell parameters,
	 *             - osymm: symmetry information,
	 *             - otuples: list of atom tuples (types and coordinates),
	 *             - oparams: list whose first four entries are the integer atom-type indices
	 *               used to select torsions and subsequent entries supply distance/angle parameters.
	 *
	 * @return A Python dict with key "tors" mapped to a list of 5-tuples (I, I, I, I, F),
	 *         where the first four entries are atom indices (unsigned ints) identifying the torsion
	 *         and the fifth entry is the torsion angle in degrees.
	 */
	static PyObject* cpplib_FindTorsionIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			Py_RETURN_NONE;
		}

		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 4> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};

		char d = type.size();
		auto value_d = FindDParamsParse<3>(self, oparams, type, d);
		auto value_a = FindATParamsParse<2>(oparams, d);
		auto value_t = FindATParamsParse<2>(oparams, d);

		auto res = WITH_LOG(FindTorsionIC, all.cell, all.symm, all.types, all.points, type, value_d, value_a, value_t[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<unsigned int>(std::get<3>(res[i])),
											 static_cast<cpplib::basic_types::FloatingPointType>(cpplib::geometry::RadtoGrad(std::get<4>(res[i])))));
		}
		return Py_BuildValue("{s:O}", "tors", lst);
	}
	/**
	 * @brief Compute torsion relationships from a set of world-coordinate atom tuples and parameters.
	 *
	 * Parses Python arguments (otuples, oparams), interprets the first four entries of oparams as atom
	 * type indices, extracts distance/angle/torsion search parameters, and finds matching torsions in
	 * the provided points.
	 *
	 * @param self Unused Python module object pointer.
	 * @param args A Python tuple (otuples, oparams):
	 *             - otuples: sequence of atom tuples (types and 3D coordinates) in world coordinates.
	 *             - oparams: list that begins with four integer atom-type indices followed by distance,
	 *               angle and torsion parameter values required by the search helpers.
	 * @return PyObject* A Python dictionary with a single key "tors" whose value is a list of 5-tuples
	 *         (I, I, I, I, F): the four atom indices (unsigned integers) that form each torsion and the
	 *         torsion angle in degrees as a float.
	 */
	static PyObject* cpplib_FindTorsionWC(PyObject* self, PyObject* args) {

		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			Py_RETURN_NONE;
		}

		Prepare_WC all(otuples);

		const std::array<int, 4> type{static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};

		char d = type.size();
		auto value_d = FindDParamsParse<3>(self, oparams, type, d);
		auto value_a = FindATParamsParse<2>(oparams, d);
		auto value_t = FindATParamsParse<2>(oparams, d);

		auto res = FindTorsionWC(all.types, all.points, type, value_d, value_a, value_t[0]);

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++) {
			PyList_Append(lst, Py_BuildValue("(IIIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<unsigned int>(std::get<3>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<4>(res[i])))));
		}
		return Py_BuildValue("{s:O}", "tors", lst);
	}
	/**
	 * @brief Compute bonds, angles, and torsions for a structure given cell and symmetry.
	 *
	 * Parses Python arguments (ocell, osymm, otuples), constructs an internal representation
	 * of the unit cell and atoms, and returns a Python dictionary containing lists of
	 * detected bonds, angles, and torsions.
	 *
	 * @param self Unused Python module/self pointer.
	 * @param args A Python tuple (ocell, osymm, otuples):
	 *        - ocell: sequence of 6 floating-point cell parameters.
	 *        - osymm: sequence of symmetry specification entries.
	 *        - otuples: sequence of atomic tuples (type and coordinates).
	 * @return PyObject* A Python dict with keys:
	 *         - "bonds": list of tuples (I, I, F) where the two integers are atom indices and
	 *           the float is the bond distance.
	 *         - "angles": list of tuples (I, I, I, F) where the three integers are atom indices
	 *           and the float is the angle in degrees.
	 *         - "tors": list of tuples (I, I, I, I, F) where the four integers are atom indices
	 *           and the float is the torsion angle in degrees.
	 */
	static PyObject* cpplib_FindDAT_IC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		if (!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuples)) {
			Py_RETURN_NONE;
		}
		Prepare_IC all(ocell, osymm, otuples);

		auto dat = FindDAT_IC(all.cell, all.symm, all.types, all.points);
		auto& datdist = std::get<0>(dat);
		auto& datang = std::get<1>(dat);
		auto& dattor = std::get<2>(dat);
		PyObject* ret = PyDict_New();

		PyObject* list_d = PyList_New(0);
		for (size_t i = 0; i < datdist.size(); i++) {
			PyList_Append(list_d, Py_BuildValue("(IIf)",
												static_cast<unsigned int>(std::get<0>(datdist[i])),
												static_cast<unsigned int>(std::get<1>(datdist[i])),
												static_cast<float>(std::get<2>(datdist[i]))));
		}
		PyDict_SetItemString(ret, "bonds", list_d);

		PyObject* list_a = PyList_New(0);
		for (size_t i = 0; i < datang.size(); i++) {
			PyList_Append(list_a, Py_BuildValue("(IIIf)",
												static_cast<unsigned int>(std::get<0>(datang[i])),
												static_cast<unsigned int>(std::get<1>(datang[i])),
												static_cast<unsigned int>(std::get<2>(datang[i])),
												static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(datang[i])))));
		}
		PyDict_SetItemString(ret, "angles", list_a);

		PyObject* list_t = PyList_New(0);
		for (size_t i = 0; i < dattor.size(); i++) {
			PyList_Append(list_t, Py_BuildValue("(IIIIf)",
												static_cast<unsigned int>(std::get<0>(dattor[i])),
												static_cast<unsigned int>(std::get<1>(dattor[i])),
												static_cast<unsigned int>(std::get<2>(dattor[i])),
												static_cast<unsigned int>(std::get<3>(dattor[i])),
												static_cast<float>(cpplib::geometry::RadtoGrad(std::get<4>(dattor[i])))));
		}
		PyDict_SetItemString(ret, "tors", list_t);

		return ret;
	}
	/**
	 * @brief Compute bond, angle, and torsion data from a list of world-coordinate atoms.
	 *
	 * Parses the provided Python sequence of atom tuples and computes detected
	 * pairwise distances (bonds), triplet angles, and quadruplet torsions, returning
	 * them packaged as Python lists inside a dictionary.
	 *
	 * @param otuples Python iterable of atom tuples; each item is expected to contain
	 *                an atom type followed by three Cartesian coordinates (x, y, z).
	 *
	 * @return PyObject* A Python dictionary with three keys:
	 *         - "bonds": list of tuples (i, j, d) where i and j are atom indices and d is
	 *                    the distance as a float.
	 *         - "angles": list of tuples (i, j, k, a) where i,j,k are atom indices and
	 *                     a is the angle in degrees as a float.
	 *         - "tors": list of tuples (i, j, k, l, t) where i,j,k,l are atom indices and
	 *                   t is the torsion angle in degrees as a float.
	 */
	static PyObject* cpplib_FindDAT_WC(PyObject* self, PyObject* otuples) {

		Prepare_WC all(otuples);

		auto dat = FindDAT_WC(all.types, all.points);
		auto& datdist = std::get<0>(dat);
		auto& datang = std::get<1>(dat);
		auto& dattor = std::get<2>(dat);
		PyObject* ret = PyDict_New();

		PyObject* list_d = PyList_New(0);
		for (size_t i = 0; i < datdist.size(); i++) {
			PyList_Append(list_d, Py_BuildValue("(IIf)",
												static_cast<unsigned int>(std::get<0>(datdist[i])),
												static_cast<unsigned int>(std::get<1>(datdist[i])),
												static_cast<float>(std::get<2>(datdist[i]))));
		}
		PyDict_SetItemString(ret, "bonds", list_d);

		PyObject* list_a = PyList_New(0);
		for (size_t i = 0; i < datang.size(); i++) {
			PyList_Append(list_a, Py_BuildValue("(IIIf)",
												static_cast<unsigned int>(std::get<0>(datang[i])),
												static_cast<unsigned int>(std::get<1>(datang[i])),
												static_cast<unsigned int>(std::get<2>(datang[i])),
												static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(datang[i])))));
		}
		PyDict_SetItemString(ret, "angles", list_a);

		PyObject* list_t = PyList_New(0);
		for (size_t i = 0; i < dattor.size(); i++) {
			PyList_Append(list_t, Py_BuildValue("(IIIIf)",
												static_cast<unsigned int>(std::get<0>(dattor[i])),
												static_cast<unsigned int>(std::get<1>(dattor[i])),
												static_cast<unsigned int>(std::get<2>(dattor[i])),
												static_cast<unsigned int>(std::get<3>(dattor[i])),
												static_cast<float>(cpplib::geometry::RadtoGrad(std::get<4>(dattor[i])))));
		}
		PyDict_SetItemString(ret, "tors", list_t);

		return ret;
	}

	/**
	 * @brief Repositions hydrogen atoms toward their nearest neighbor according to a himp mapping.
	 *
	 * For each atom with type equal to 1, finds its nearest neighboring atom and moves the hydrogen
	 * along the vector from that neighbor toward the hydrogen by a fraction specified in the himp map.
	 * The himp argument may be a single float (applied to all neighbor types) or a list of floats
	 * indexed by atom type. If the himp list does not contain an entry for a required neighbor type,
	 * the function returns an error dictionary describing the missing entry.
	 *
	 * @returns PyObject* A Python dictionary. On success: {"atoms": <list of atom coordinate tuples>} with updated positions.
	 * On error: {"error_str": "<message>"}.
	 */
	static PyObject* cpplib_himp(PyObject* self, PyObject* args) {
		// args [2] = [(type,x,y,z), ... ], 
		//            length : cpplib::basic_types::FloatingPointType or [cpplib::basic_types::FloatingPointType, ...]

		PyObject* o_tuple = NULL;
		PyObject* o_himp = NULL;
		std::vector<FloatingPointType> himp;

		if (!PyArg_ParseTuple(args, "OO", &o_tuple, &o_himp)) {
			Py_RETURN_NONE;
		}
		Prepare_WC all(o_tuple);

		bool simplehimp = PyFloat_CheckExact(o_himp);
		if (simplehimp) {
			himp.resize(cpplib::constants::mend_size, static_cast<FloatingPointType>(PyFloat_AsDouble(o_himp)));
		}
		else {
			pyListToVectorFloat(o_himp, &himp);
		}
		const auto himp_s = himp.size();
		// Code section
		const auto s = all.types.size();
		for (int i = 0; i < s; i++) {
			if (all.types[i] != 1) continue;
			auto dist = static_cast<FloatingPointType>(INFINITY);
			int best = i;
			for (int j = 0; j < s; j++) {
				if (j == i) continue;
				auto temp = PointType::distance(all.points[i], all.points[j]);
				if (temp < dist) {
					dist = temp;
					best = j;
				}
			}
			if (himp_s <= all.types[best]) {
				std::string err = std::string("Too short himp list: type ") + std::to_string(all.types[best]) + " is not exist.";
				return Py_BuildValue("{s:s}", "error_str", err.c_str());
			}
			all.points[i] = (all.points[best] + ((all.points[i] - all.points[best]) * (himp[all.types[best]] / dist)));
		}

		// returns List[Tuple(atom1, atom2), ...] 
		return Py_BuildValue("{s:N}",
							 "atoms", py_util::convert(all.points));
	}

	/**
	 * @brief Perform a substructure search between two molecule strings.
	 *
	 * Parses two C strings from the Python call, interprets the first as the query
	 * molecule and the second as the database/target molecule, runs a full search,
	 * and returns whether the query was found inside the target.
	 *
	 * @param self Unused Python module/object pointer.
	 * @param args Python argument tuple containing two strings: the query and target.
	 * @return PyObject* Python True if the query matches the target, Python False if
	 * the query does not match, or Python None if argument parsing failed.
	 */
	static PyObject* cpplib_SubSearch(PyObject* self, PyObject* args) {

		const char* s1 = NULL;
		const char* s2 = NULL;
		if (!PyArg_ParseTuple(args, "ss", &s1, &s2)) {
			Py_RETURN_NONE;
		}
		bool res = false;

		cpplib::SearchGraph graph;
		auto&& inputpair = cpplib::MoleculeParser<AtomTypeRequest>::Read(s1);
		graph.setupInput(std::move(inputpair.first));
		auto datg = cpplib::MoleculeParser<AtomTypeRequest>::Read(s2).first.makeCopyEx<AtomTypeData>();
		graph.setupData(std::move(datg));
		graph.prepareToSearch();
		if (graph.startFullSearch(false)) {
			Py_RETURN_TRUE;
		}
		else {
			Py_RETURN_FALSE;
		}
	}

	/**
	 * @brief Compare and compact a structure using the Compaq routine and return errors and XYZ block.
	 *
	 * Parses a Python tuple of (ocell, osymm, otuple), initializes distance data, builds an internal
	 * Prepare_IC from the inputs, runs Compaq (with logging) and returns the resulting errors list
	 * and XYZ block.
	 *
	 * @param self Unused Python module/object pointer.
	 * @param args Python tuple containing (ocell, osymm, otuple).
	 * @return PyObject* A Python dictionary with keys:
	 *         - "errors": list of error descriptions produced by Compaq.
	 *         - "xyz_block": XYZ representation of the (possibly modified) structure.
	 */
	static PyObject* cpplib_compaq(PyObject* self, PyObject* args) {
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuple = NULL;

		if (!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuple)) {
			Py_RETURN_NONE;
		}
		useDistances(self);

		Prepare_IC all(ocell, osymm, otuple);

		auto [ret_v, ret_l] = WITH_LOG(Compaq, all.cell, all.symm, all.types, all.points);


		return Py_BuildValue("{s:N,s:N}",
							 "errors", py_util::convert(ret_l),
							 "xyz_block", py_util::convert(ret_v));
	}
	/**
	 * @brief Reorders and normalizes a serialized AtomTypeData database string.
	 *
	 * Takes a Python Unicode string representing an AtomTypeData database, sorts/restructures
	 * its entries into a canonical form, and returns the resulting Unicode string.
	 *
	 * @param arg Python Unicode string containing the database to sort.
	 * @return PyObject* Python Unicode string with the sorted and normalized database.
	 */
	static PyObject* cpplib_SortDatabase(PyObject* self, PyObject* arg) {
		Py_ssize_t us;
		const auto ret = cpplib::MoleculeParser<cpplib::currents::AtomTypeData>::ResortString(PyUnicode_AsUTF8AndSize(arg, &us));
		return PyUnicode_FromString(ret.c_str());
	}

	/**
	 * @brief Create clusters from atomic coordinates and anchor points within a unit cell.
	 *
	 * Parses Python arguments (cell, symmetry, atom tuples, anchor list, radius), builds an internal
	 * representation, runs the clustering routine, and returns cluster points plus a polymer flag.
	 *
	 * @param self Unused Python module/self pointer.
	 * @param args A Python tuple: (cell, symm, tuples, anchors, radius)
	 *   - cell: unit cell parameters (sequence of six floats).
	 *   - symm: symmetry labels/list corresponding to the cell.
	 *   - tuples: list of atom tuples (type and coordinates).
	 *   - anchors: list of 4-tuples (x, y, z, weight) specifying anchor points and their weights.
	 *   - radius: floating-point over-radius used for clustering.
	 *
	 * @return A Python dictionary with keys:
	 *   - "points": converted cluster points (implementation-specific structure).
	 *   - "hasPolymer": Python boolean set to `True` if a polymeric component was detected, `False` otherwise.
	 *
	 * If argument parsing fails, the function returns Python None.
	 */
	static PyObject* cpplib_ClusterCreate(PyObject* self, PyObject* args) {
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* ocoords = NULL;
		cpplib::basic_types::FloatingPointType over_radius = 0;
		if (!PyArg_ParseTuple(args, "OOOOf", &ocell, &osymm, &otuples, &ocoords, &over_radius)) {
			Py_RETURN_NONE;
		}

		LOG_INTERFACE_GUARD("cpplib_ClusterCreate");
		Prepare_IC all(ocell, osymm, otuples);

		Py_ssize_t s = PyList_Size(ocoords);
		std::vector<Cluster::AnchorType> anchors;
		anchors.reserve(static_cast<size_t>(s));

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* o_tuple = PyList_GetItem(ocoords, i);
			anchors.emplace_back(cpplib::geometry::Point<FloatingPointType>(
				static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 0))),
				static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 1))),
				static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 2)))),
				static_cast<cpplib::basic_types::FloatingPointType>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 3))));
		}
		bool hasPolymer = false;
		auto ret = WITH_LOG(ClusterCreate, all.cell, all.symm, all.types, all.points, anchors, over_radius, hasPolymer);

		return Py_BuildValue("{s:N,s:N}",
							 "points", py_util::convert(ret),
							 "hasPolymer", PyBool_FromLong(hasPolymer?1:0));
	}

	/**
	 * @brief Compute Voronoi cells for an atomic configuration and return cell and unit-cell data.
	 *
	 * Parses a Python argument tuple of the form (cell, symm, tuples, bools, cutoff) and builds a symmetry-aware
	 * unit cell, computes spatial bonds, constructs a Voronoi diagram, and packages the result for Python.
	 *
	 * The expected Python arguments are:
	 *  - cell: sequence of six floats describing the unit cell (a, b, c, alpha, beta, gamma)
	 *  - symm: list of symmetry descriptors (strings or symmetry objects accepted by Prepare_IC)
	 *  - tuples: list of atomic tuples (type and Cartesian coordinates) used to build the asymmetric unit
	 *  - bools: list of integers (0 or 1) with one entry per asymmetric-unit site indicating whether the
	 *           corresponding site should be treated as included when computing Voronoi cells
	 *  - cutoff: float radius used to build the spatial grid for bond detection
	 *
	 * @param self Unused Python module/object pointer (standard C API convention).
	 * @param args Python tuple described above.
	 * @return PyObject* A new reference to a Python dictionary with keys:
	 *   - "voronoi_cells": an object produced by py_util::convert(vf) representing fused Voronoi polyhedra
	 *   - "unit_cell": an object produced by py_util::convert(buildresult.atoms) representing the unit-cell atoms
	 *   Returns Python None if argument parsing fails.
	 */
	static PyObject* cpplib_Voronoi(PyObject* self, PyObject* args) {
		using Diagram = cpplib::voronoi::VoronoiDiagram;

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* obools = NULL;
		float cutoff = 6.0;
		if (!PyArg_ParseTuple(args, "OOOOf", &ocell, &osymm, &otuples, &obools, &cutoff)) {
			Py_RETURN_NONE;
		}

		LOG_INTERFACE_GUARD("cpplib_Voronoi");

		Prepare_IC all(ocell, osymm, otuples);
		auto ps = all.points.size();
		std::vector<int> intbools;
		intbools.reserve(ps);
		pyListToVectorInt(obools, &intbools);
		std::vector<bool> bools(ps);

		for (int i = 0; i < ps; i++) {
			bools[i] = intbools[i] != 0;
		}

		geometry::Cell cell(all.cell);

		bools.resize(all.points.size(), false);

		std::vector<geometry::Symm<FloatingPointType>> symmvec;
		symmvec.reserve(all.symm.size());
		for (int i = 0; i < all.symm.size(); i++)
		{
			symmvec.emplace_back(all.symm[i]);
		}

		cluster_detail::UnitCellBuilder ucb(symmvec);
		auto buildresult = ucb.build(all.points, all.types);


		geometry::SpatialGrid<FloatingPointType> space;
		space.build(buildresult.atoms.points, cell, cutoff);
		auto bonds = WITH_LOG_M(space, get_bonds, false);

		// Flags intentionally correspond only to the asymmetric-unit inputs; 
		// VoronoiDiagram resizes the flag vector and treats symmetry-expanded sites as false.
		Diagram diag(buildresult.atoms.points,bonds, cell, bools);

		auto ce = diag.extractCells();
				
		cpplib::voronoi::VoronoiFused vf(ce, cell.fracToCart());
		vf.polyhedra.resize(all.points.size());

		// Build return value
		return Py_BuildValue("{s:N,s:N}",
							 "voronoi_cells", py_util::convert(vf),
							 "unit_cell", py_util::convert(buildresult.atoms));
	}

	static struct PyMethodDef methods[] = {
		{ "GenBonds", cpplib_GenBonds, METH_O, "Generate bond list"},
		{ "GenBondsEx", cpplib_GenBondsEx, METH_O, "Generate bond list with length"},
		{ "SearchMain", cpplib_SearchMain, METH_VARARGS, "Compare graph with data"},
		{ "CompareGraph", cpplib_CompareGraph, METH_VARARGS, "Compare two graphs"},
		{ "FindMoleculesInCell", cpplib_FindMoleculesInCell, METH_VARARGS, "Create graph from cell"},
		{ "FindMoleculesWithoutCell", cpplib_FindMoleculesWithoutCell, METH_O, "Create graph from xyz"},
		{ "GenSymm", cpplib_GenSymm, METH_VARARGS, "Generates symmetry by symm code"},
		{ "FindDistanceIC", cpplib_FindDistanceIC, METH_VARARGS, "Find distances with current parameters in cell"},
		{ "FindDistanceWC", cpplib_FindDistanceWC, METH_VARARGS, "Find distances with current parameters in xyz"},
		{ "FindAngleIC", cpplib_FindAngleIC, METH_VARARGS, "Find angles with current parameters in cell"},
		{ "FindAngleWC", cpplib_FindAngleWC, METH_VARARGS, "Find angles with current parameters in xyz"},
		{ "FindTorsionIC", cpplib_FindTorsionIC, METH_VARARGS, "Find torsions with current parameters in cell"},
		{ "FindTorsionWC", cpplib_FindTorsionWC, METH_VARARGS, "Find torsions with current parameters in xyz"},
		{ "FindDAT_IC", cpplib_FindDAT_IC, METH_VARARGS, "Create dictionary with distances, angles and torsions in cell"},
		{ "FindDAT_WC", cpplib_FindDAT_WC, METH_O, "Create dictionary with distances, angles and torsions in xyz"},
		{ "himp", cpplib_himp, METH_VARARGS, "Moves hydrogens to the nearest atom"},
		{ "SubSearch", cpplib_SubSearch, METH_VARARGS, "Compare two graphs"},
		{ "compaq", cpplib_compaq, METH_VARARGS, "Do the same as Olex2 'compaq' function"},
		{ "SortDatabase", cpplib_SortDatabase, METH_O, "Sort graph"},
		{ "Cluster", cpplib_ClusterCreate, METH_VARARGS, "Create cluster"},
		{ "VoronoiCalculation", cpplib_Voronoi, METH_VARARGS, "Calculate Voronoi cells"},


		{ NULL, NULL, 0, NULL }
	};

	static PyModuleDef cpplib_module = {
		PyModuleDef_HEAD_INIT, "cpplib", NULL, -1, methods,
		NULL, NULL, NULL, NULL
	};

	PyMODINIT_FUNC PyInit_cpplib(void) {
		return PyModule_Create(&cpplib_module);
	}
}