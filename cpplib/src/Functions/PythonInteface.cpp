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
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "AllInOneAndCurrent.h"
#include "Functions.h"
#include "../Classes/Interfaces.h"
#include "../BaseHeaders/DebugMes.h"
#include <list>

using namespace cpplib::currents;

// Utilities section

enum class ErrorState {
	OK = 0,
	UnknownError,
};
struct Prepare_WC {
	std::vector<cpplib::currents::AtomTypeData> types;
	std::vector<cpplib::currents::PointType> points;
	Prepare_WC(PyObject* otuples) {
		Py_ssize_t s = PyList_Size(otuples);
		types.reserve(static_cast<size_t>(s));
		points.reserve(static_cast<size_t>(s));

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* o_tuple = PyList_GetItem(otuples, i);
			types.push_back(static_cast<AtomTypeData>(PyLong_AsLong(PyTuple_GetItem(o_tuple, 0))));
			points.emplace_back(static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 1))),
								static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 2))),
								static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 3))));
		}
	}
};
struct Prepare_IC : public Prepare_WC {
	std::array<float, 6> cell;
	std::vector<const char*> symm;
	Prepare_IC(PyObject* ocell, PyObject* osymm, PyObject* otuples) : Prepare_WC(otuples) {
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}
		const Py_ssize_t s = PyList_Size(osymm);
		deb_write("pyListToVectorCharP s =", s);
		symm.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			symm[i] = PyUnicode_AsUTF8(PyList_GetItem(osymm, i));
		}
	}
};

extern "C" {
inline static void useDistances(PyObject * self);
}

template <char times>
static std::array<std::pair<float, float>, times> FindDParamsParse(PyObject* self, PyObject* oparams, const std::array<int, times+1> type, char& d) {
	std::array<std::pair<float, float>, times> value;
	for (char i = 0; i < times; i++, d += 2)
	{
		value[i].first = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oparams, d)));
		value[i].second = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oparams, d + 1)));

		if (value[i].first == 0) {
			useDistances(self);
			value[i].first = p_distances->minDistance(type[i], type[i + 1]);
		}

		if (value[i].second == 0) {
			useDistances(self);
			value[i].second = p_distances->maxDistance(type[i], type[i + 1]);
		}
	}
	return value;
}
template <char times>
static std::array<std::pair<float, float>, times> FindATParamsParse(PyObject* oparams, char& d) {
	std::array<std::pair<float, float>, times> value;
	for (char i = 0; i < times; i++, d += 2)
	{
		value[i].first = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oparams, d)));
		value[i].second = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oparams, d + 1)));
	}
	return value;
}
extern "C" {

	inline static ErrorState pyListToVectorFloat(PyObject* plist, std::vector<float>* pret) {
		const Py_ssize_t s = PyList_Size(plist);
		std::vector<float>& ret = *pret;
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			ret[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(plist, i)));
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
	inline static ErrorState pyListToVectorCharP(PyObject* plist, std::vector<const char*>* pret) {

		std::vector<const char*>& ret = *pret;
		deb_write("pyListToVectorCharP called");
		const Py_ssize_t s = PyList_Size(plist);
		deb_write("pyListToVectorCharP s =", s);
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			deb_write("pyListToVectorCharP start i = ", i);

			ret[i] = PyUnicode_AsUTF8(PyList_GetItem(plist, i));
			deb_write("pyListToVectorCharP end   i = ", i);
		}
		deb_write("pyListToVectorCharP return");
		return ErrorState::OK;
	}
	inline static ErrorState pyTXYZparse(PyObject* o_list, std::vector<AtomTypeData>* types, std::vector<PointType>* points) {
		Py_ssize_t s = PyList_Size(o_list);
		types->clear();
		types->reserve(static_cast<int>(s));
		points->clear();
		points->reserve(static_cast<int>(s));
		
		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* o_tuple = PyList_GetItem(o_list, i);
			//if (o_tuple == NULL) return ErrorState::UnknownError;
			types->push_back(static_cast<AtomTypeData>(PyLong_AsLong(PyTuple_GetItem(o_tuple, 0))));
			points->emplace_back(static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 1))),
								 static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 2))),
								 static_cast<float>(PyFloat_AsDouble(PyTuple_GetItem(o_tuple, 3))));
		}
		return ErrorState::OK;
	}
	inline static void useDistances(PyObject* self) {
		deb_write("useDistances check p_dist");
		if (p_distances != nullptr) {
			return;
		}
		deb_write("useDistances parse __file__");
		std::string full(PyUnicode_AsUTF8(PyObject_GetAttrString(self, "__file__")));

		auto found = full.find_last_of("\\/");
		auto bond_filename = full.substr(0, found + 1) + "BondLength.ini";

		deb_write("Create dist");
		static DistancesType dist(bond_filename);
		p_distances = &dist;
	}

	// Python section
	static PyObject* cpplib_GenBonds(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *(p_distances);
		Prepare_WC all(arg);
		
		FAMStructType famstr(std::move(all.types), std::move(all.points));
		std::string errM;
		auto&& bonds = famstr.findBonds(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < bonds.size(); i++)
		{
			PyList_Append(lst, Py_BuildValue("(ll)",
											 static_cast<long>(bonds[i].first),
											 static_cast<long>(bonds[i].second)));
		}
		return Py_BuildValue("{s:O}",
							 "bonds", lst);
	}
	static PyObject* cpplib_GenBondsEx(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *(p_distances);
		Prepare_WC all(arg);
		FAMStructType famstr(std::move(all.types), std::move(all.points));
		std::string errM;
		auto&& bonds = famstr.findBondsEx(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < bonds.size(); i++)
		{
			PyList_Append(lst, Py_BuildValue("(llf)",
											 static_cast<long>(bonds[i].first),
											 static_cast<long>(bonds[i].second),
											 static_cast<float>(bonds[i].length)));
		}
		return Py_BuildValue("{s:O}",
							 "bonds", lst);
	}

	static PyObject* cpplib_SearchMain(PyObject* self, PyObject* args) {
		const char* search = NULL;
		PyObject* o = NULL;
		int np = 0;
		int exact = 0;
		if (!PyArg_ParseTuple(args, "sOip", &search, &o, &np, &exact)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		const Py_ssize_t s = PyList_Size(o);
		const int ds = static_cast<int>(s);

		deb_write("search = ", search);
		deb_write("np = ", np);
		deb_write("exact = ", exact);
		std::vector< const char*> data; pyListToVectorCharP(o, &data);
		deb_write("data.size = ", data.size());
		deb_write("py_SearchMain invoke SearchMain");
		const auto ret = SearchMain(search, std::move(data), np, (exact != 0));
		deb_write("py_SearchMain closes SearchMain");
		const auto ret_s = ret.size();
		PyObject* ret_o = PyList_New(0);
		deb_write("py_SearchMain create return list");

		for (Py_ssize_t i = 0; i < ret_s; i++)
		{
			PyList_Append(ret_o, PyLong_FromLong(ret[i]));
		}
		deb_write("py_SearchMain return");
		return ret_o;
	}
	static PyObject* cpplib_CompareGraph(PyObject* self, PyObject* args)
	{
		deb_write("cpplib_CompareGraph: start");
		const char* s1 = NULL;
		const char* s2 = NULL;
		int b = 0;
		deb_write("cpplib_CompareGraph: arg parse start");
		if(!PyArg_ParseTuple(args, "ssp", &s1, &s2, &b)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		deb_write("s1 = ", s1);
		deb_write("s2 = ", s2);
		deb_write("exact = ", b);
		deb_write("cpplib_CompareGraph: invoke CompareGraph");
		if (CompareGraph(s1, s2, (b != 0))) {
			Py_RETURN_TRUE;
		}
		else {
			Py_RETURN_FALSE;
		}
	}

	static PyObject* cpplib_FindMoleculesInCell(PyObject* self, PyObject* args) {
		useDistances(self);
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuple = NULL;
		if(!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuple)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}

		Prepare_IC all(ocell, osymm, otuple);

		auto ret = FindMoleculesInCell(all.cell, all.symm, all.types, all.points);
		PyObject* o_xyz_block = PyList_New(0);

		for (auto & mol : std::get<2>(ret))
		{
			PyObject* o_molecule = PyList_New(0);
			for (auto& atom : std::get<0>(mol)) {
				PyObject* o_atom = Py_BuildValue("{s:f,s:f,s:f,s:l}",
												 "x", float(std::get<0>(atom).get(0)),
												 "y", float(std::get<0>(atom).get(1)),
												 "z", float(std::get<0>(atom).get(2)),
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
		// List[Tuple(atom1, atom2), ...] под ключом 'bonds' 
		return Py_BuildValue("{s:s,s:s,s:O}",
							 "graph_str", std::get<0>(ret).c_str(),
							 "error_str", std::get<1>(ret).c_str(),
							 "xyz_block", o_xyz_block);
	}
	static PyObject* cpplib_FindMoleculesWithoutCell(PyObject* self, PyObject* otuple) {
		useDistances(self);

		Prepare_WC all(otuple);

		deb_write("cpplib_FindMoleculesWithoutCell invoke FindMoleculesInCell");
		auto ret = FindMoleculesWithoutCell(all.types, all.points);
		deb_write("cpplib_FindMoleculesWithoutCell returned from FindMoleculesInCell");
		PyObject* o_xyz_block = PyList_New(0);

		for (auto& mol : std::get<2>(ret))
		{
			PyObject* o_molecule = PyList_New(0);
			for (auto& atom : std::get<0>(mol)) {
				PyObject* o_atom = Py_BuildValue("{s:f,s:f,s:f,s:l}",
												 "x", float(std::get<0>(atom).get(0)),
												 "y", float(std::get<0>(atom).get(1)),
												 "z", float(std::get<0>(atom).get(2)),
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
		// List[Tuple(atom1, atom2), ...] под ключом 'bonds' 
		return Py_BuildValue("{s:s,s:s,s:O}",
							 "graph_str", std::get<0>(ret).c_str(),
							 "error_str", std::get<1>(ret).c_str(),
							 "xyz_block", o_xyz_block);
	}

	static PyObject* cpplib_GenSymm(PyObject* self, PyObject* args) {
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		unsigned char flags;

		deb_write("cpplib_GenSymm: Parsing start");
		if(!PyArg_ParseTuple(args, "OBO", &otuples, &flags, &osymm)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		bool movetocell = flags & 1;
		bool movemasstocell = flags & 2;

		Prepare_WC all(otuples);
		deb_write("cpplib_GenSymm: atom parsing ended");

		deb_write("cpplib_GenSymm: pyListToVectorCharP invoking");
		std::vector<const char*> nsymm; pyListToVectorCharP(osymm, &nsymm);
		deb_write("cpplib_GenSymm: pyListToVectorCharP returned");

		deb_write("cpplib_GenSymm: symm parsing started");
		std::vector<SymmType> symm;
		const size_t ss = nsymm.size();
		for (size_t i = 0; i < ss; i++)
		{
			symm.emplace_back(nsymm[i]);
		}
		deb_write("cpplib_GenSymm: Parsing ended");

		// Move center of mass to cell
		const auto s_points = all.points.size();
		if (movemasstocell) {
			deb_write("cpplib_GenSymm: Move center of mass started");
			PointType centerofmass(0, 0, 0);
			for (size_t i = 0; i < s_points; i++)
			{
				centerofmass += all.points[i];
			}
			centerofmass /= s_points;
			PointType ceilmass(std::ceil(centerofmass.get(0)), 
							   std::ceil(centerofmass.get(1)), 
							   std::ceil(centerofmass.get(2)));

			for (size_t i = 0; i < ss; i++)
			{
				PointType movedcenter = symm[i].GenSymm(centerofmass);
				PointType ceilmoved(std::ceil(movedcenter.get(0)), 
										std::ceil(movedcenter.get(1)), 
										std::ceil(movedcenter.get(2)));
				symm[i].point += ceilmass - ceilmoved;
			}

			deb_write("cpplib_GenSymm: Move center of mass ended");
		}

		FAMStructType famstr(std::move(all.types), std::move(all.points));
		FAMCellType fcell(FAMCellType::base(32, 32, 32, 90, 90, 90, true));
		fcell.GenerateSymm(famstr, symm, movetocell);

		deb_write("cpplib_GenSymm: famstr.types.size() = ", famstr.types.size());
		deb_write("cpplib_GenSymm: famstr.points.size() = ", famstr.points.size());
		deb_write("cpplib_GenSymm: famstr.parseIndex.size() = ", famstr.parseIndex.size());
		deb_write("cpplib_GenSymm: famstr.sizePoints = ", famstr.sizePoints);

		for (Py_ssize_t i = s_points; i < famstr.sizePoints; i++)
		{
			PyList_Append(otuples, Py_BuildValue("(Ifff)",
												 static_cast<unsigned int>(famstr.types[famstr.parseIndex[i]]),
												 static_cast<float>(famstr.points[i].get(0)),
												 static_cast<float>(famstr.points[i].get(1)),
												 static_cast<float>(famstr.points[i].get(2))));
		}
		Py_INCREF(otuples);
		return otuples;
	}

	static PyObject* cpplib_FindDistanceIC(PyObject* self, PyObject* args) {
		deb_write("cpplib_FindDistanceIC: invoked");

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;

		if(!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 2> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		char d = 2;
		auto value = FindDParamsParse<1>(self, oparams, type, d);

		deb_write("cpplib_FindDistanceIC: all.types.size: ", all.types.size());
		deb_write("cpplib_FindDistanceIC: all.points.size: ", all.points.size());
		deb_write("cpplib_FindDistanceIC: param type1: ", type[0]);
		deb_write("cpplib_FindDistanceIC: param type2: ", type[1]);
		deb_write("cpplib_FindDistanceIC: param min: ", value[0].first);
		deb_write("cpplib_FindDistanceIC: param max: ", value[0].second);
		deb_write("cpplib_FindDistanceIC: invoke FindDistanceIC");
		auto res =  FindDistanceIC(all.cell, all.symm, all.types, all.points, type, value[0]);
		deb_write("cpplib_FindDistanceIC: FindDistanceIC returned");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<float>(std::get<2>(res[i]))));
		}
		deb_write("cpplib_FindDistanceIC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "distances", lst);
	}
	static PyObject* cpplib_FindDistanceWC(PyObject* self, PyObject* args) {
		deb_write("cpplib_FindDistanceWC: invoked");
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if(!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		Prepare_WC all(otuples);

		const std::array<int, 2> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		char d = 2;
		auto value = FindDParamsParse<1>(self, oparams, type, d);


		deb_write("cpplib_FindDistanceWC: all.types.size: ", all.types.size());
		deb_write("cpplib_FindDistanceWC: all.points.size: ", all.points.size());
		deb_write("cpplib_FindDistanceWC: param type1: ", type[0]);
		deb_write("cpplib_FindDistanceWC: param type2: ", type[1]);
		deb_write("cpplib_FindDistanceWC: param min: ", value[0].first);
		deb_write("cpplib_FindDistanceWC: param max: ", value[0].second);
		deb_write("cpplib_FindDistanceWC: invoke FindDistanceWC");
		auto res = FindDistanceWC(all.types, all.points, type, value[0]);
		deb_write("cpplib_FindDistanceWC: FindDistanceWC return");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<float>(std::get<2>(res[i]))));
		}
		deb_write("cpplib_FindDistanceWC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "distances", lst);
	}

	static PyObject* cpplib_FindAngleIC(PyObject* self, PyObject* args) {
		deb_write("cpplib_FindAngleIC: invoked");
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 3> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		auto value_d = FindDParamsParse<2>(self, oparams, type, d);
		auto value_a = FindATParamsParse<1>(oparams, d);


		deb_write("cpplib_FindAngleIC: invoke FindAngleIC");
		auto res = FindAngleIC(all.cell, all.symm, all.types, all.points, type, value_d, value_a[0]);
		deb_write("cpplib_FindAngleIC: FindAngleIC return");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(res[i])))));
		}
		deb_write("cpplib_FindAngleIC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "angles", lst);
	}
	static PyObject* cpplib_FindAngleWC(PyObject* self, PyObject* args) {
		deb_write("cpplib_FindAngleWC: invoked");
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}

		Prepare_WC all(otuples);

		const std::array<int, 3> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
			                           static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
									   static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		auto value_d = FindDParamsParse<2>(self, oparams, type, d);
		auto value_a = FindATParamsParse<1>(oparams, d);


		deb_write("cpplib_FindAngleWC: invoke FindAngleWC");
		auto res = FindAngleWC(all.types, all.points, type, value_d, value_a[0]);
		deb_write("cpplib_FindAngleWC: FindAngleWC return");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<3>(res[i])))));
		}
		deb_write("cpplib_FindAngleWC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "angles", lst);
	}
	static PyObject* cpplib_FindTorsionIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}

		Prepare_IC all(ocell, osymm, otuples);

		const std::array<int, 4> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};

		char d = type.size();
		auto value_d = FindDParamsParse<3>(self, oparams, type, d);
		auto value_a = FindATParamsParse<2>(oparams, d);
		auto value_t = FindATParamsParse<2>(oparams, d);

		deb_write("cpplib_FindTorsionIC: invoke FindTorsionIC");
		auto res = FindTorsionIC(all.cell, all.symm, all.types, all.points, type, value_d, value_a, value_t[0]);
		deb_write("cpplib_FindTorsionIC: FindTorsionIC return");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<unsigned int>(std::get<3>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<4>(res[i])))));
		}
		deb_write("cpplib_FindTorsionIC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "tors", lst);
	}
	static PyObject* cpplib_FindTorsionWC(PyObject* self, PyObject* args) {

		PyObject* otuples = NULL;
		PyObject* oparams = NULL;
		if (!PyArg_ParseTuple(args, "OO", &otuples, &oparams)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}

		Prepare_WC all(otuples);

		const std::array<int, 4> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};

		char d = type.size();
		auto value_d = FindDParamsParse<3>(self, oparams, type, d);
		auto value_a = FindATParamsParse<2>(oparams, d);
		auto value_t = FindATParamsParse<2>(oparams, d);

		deb_write("cpplib_FindTorsionWC: invoke FindTorsionWC");
		auto res = FindTorsionWC(all.types, all.points, type, value_d, value_a, value_t[0]);
		deb_write("cpplib_FindTorsionWC: FindTorsionWC return");

		auto res_s = res.size();
		PyObject* lst = PyList_New(0);
		for (size_t i = 0; i < res_s; i++)
		{
			PyList_Append(lst, Py_BuildValue("(IIIIf)",
											 static_cast<unsigned int>(std::get<0>(res[i])),
											 static_cast<unsigned int>(std::get<1>(res[i])),
											 static_cast<unsigned int>(std::get<2>(res[i])),
											 static_cast<unsigned int>(std::get<3>(res[i])),
											 static_cast<float>(cpplib::geometry::RadtoGrad(std::get<4>(res[i])))));
		}
		deb_write("cpplib_FindTorsionWC: List[...] completed. size = ", res_s);
		return Py_BuildValue("{s:O}", "tors", lst);
	}
	static PyObject* cpplib_FindDAT_IC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuples = NULL;
		if (!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuples)) {
			deb_write("! Critic Error: Parse Error - return None");
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
	static PyObject* cpplib_FindDAT_WC(PyObject* self, PyObject* otuples) {
		deb_write("cpplib_FindDAT_WC started");

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

	static PyObject* cpplib_himp(PyObject* self, PyObject* args) {
		// args [2] = [(type,x,y,z), ... ], 
		//            length : float or [float, ...]

		PyObject* o_tuple = NULL;
		PyObject* o_himp = NULL;
		std::vector<FloatingPointType> himp;

		if (!PyArg_ParseTuple(args, "OO", &o_tuple, &o_himp)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		Prepare_WC all(o_tuple);

		bool simplehimp = PyFloat_CheckExact(o_himp);
		if (simplehimp) {
			himp.resize(cpplib::mend_size, static_cast<FloatingPointType>(PyFloat_AsDouble(o_himp)));
		}
		else {
			pyListToVectorFloat(o_himp, &himp);
		}
		const int himp_s = himp.size();
		// Code section
		const int s = all.types.size();
		for (int i = 0; i < s; i++)
		{
			if (all.types[i] != 1) continue;
			FloatingPointType dist = INFINITY;
			int best = i;
			for (int j = 0; j < s; j++)
			{
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
			all.points[i] = (all.points[best] + ((all.points[i] - all.points[best])* (himp[all.types[best]] / dist)));
		}
		
		// Return section
		PyObject* o_xyz_block = PyList_New(0);
		PyObject* o_ret = PyDict_New();
		for (int i = 0; i < s; i++)
		{
			PyObject* o_atom = Py_BuildValue("(fff)",
											 static_cast<float>(all.points[i].get(0)),
											 static_cast<float>(all.points[i].get(1)),
											 static_cast<float>(all.points[i].get(2)));
			PyList_Append(o_xyz_block, o_atom);
		}
		// List[Tuple(atom1, atom2), ...] под ключом 'bonds' 
		return Py_BuildValue("{s:O}",
							 "atoms", o_xyz_block);
	}

	static PyObject* cpplib_SubSearch(PyObject* self, PyObject* args) {
		deb_write("cpplib_SubSearch started");

		const char* s1 = NULL;
		const char* s2 = NULL;
		if (!PyArg_ParseTuple(args, "ss", &s1, &s2)) {
			deb_write("! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}
		deb_write("s1 = ", s1);
		deb_write("s2 = ", s2);
		deb_write("cpplib_SubSearch: invoke CompareGraph");
		bool res = false;

		cpplib::SearchGraph graph;

		auto&& inputpair = SearchGraphType::RequestGraphType::ReadInput(s1);
		graph.setupInput(std::move(inputpair.first));
		deb_write("cpplib_SubSearch start ReadData");
		auto datg = SearchGraphType::RequestGraphType::ReadInput(s2).first.makeCopyEx<AtomTypeData>();
		graph.setupData(std::move(datg));
		deb_write("cpplib_SubSearch start prepareSearch");
		graph.prepareToSearch();
		deb_write("cpplib_SubSearch start FullSearch");
		if( graph.startFullSearch(false)){
			Py_RETURN_TRUE;
		}
		else {
			Py_RETURN_FALSE;
		}
	}

	static PyObject* cpplib_compaq(PyObject* self, PyObject* args) {
		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otuple = NULL;

		deb_write("cpplib_compaq argument parsing start");
		if (!PyArg_ParseTuple(args, "OOO", &ocell, &osymm, &otuple)) {
			deb_write("cpplib_compaq! Critic Error: Parse Error - return None");
			Py_RETURN_NONE;
		}

		deb_write("cpplib_compaq call useDistances");
		useDistances(self);

		deb_write("cpplib_compaq call Prepare_IC");
		Prepare_IC all(ocell, osymm, otuple);

		deb_write("cpplib_compaq call Compaq");
		auto ret = Compaq(all.cell, all.symm, all.types, all.points);
		deb_write("cpplib_compaq Compaq successful");

		PyObject* o_xyz_block = PyList_New(0);
		PyObject* o_errors = PyList_New(0);

		deb_write("cpplib_compaq create List[txyz]");
		deb_write("cpplib_compaq std::get<0>(ret).size() = ", std::get<0>(ret).size());
		for (int i = 0; i < std::get<0>(ret).size(); i++)
		{
			PyObject* o_atom = Py_BuildValue("(fff)",
											 static_cast<float>(std::get<0>(ret)[i].get(0)),
											 static_cast<float>(std::get<0>(ret)[i].get(1)),
											 static_cast<float>(std::get<0>(ret)[i].get(2)));
			PyList_Append(o_xyz_block, o_atom);
		}

		deb_write("cpplib_compaq create List[error_str]");
		while(std::get<1>(ret).empty() == false)
		{
			PyList_Append(o_errors, Py_BuildValue("s", (std::get<1>(ret)).front().c_str()));
			(std::get<1>(ret)).pop_front();
		}

		deb_write("cpplib_compaq return Dict {errors, xyz_block}");
		return Py_BuildValue("{s:O,s:O}",
							 "errors", o_errors,
							 "xyz_block", o_xyz_block);
	}
	static PyObject* cpplib_SortDatabase(PyObject* self, PyObject* arg) {
		const auto ret = cpplib::currents::SearchGraphType::DatabaseGraphType::ResortString(PyUnicode_AsUTF8(arg));
		return PyUnicode_FromString(ret.c_str());
	}
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

	{ NULL, NULL, 0, NULL }
};

static struct PyModuleDef cpplib_module = {
	PyModuleDef_HEAD_INIT, "cpplib", NULL, -1, methods,
	NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_cpplib(void)
{
	return PyModule_Create(&cpplib_module);
}