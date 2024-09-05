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

using namespace cpplib::currents;

// Utilities section

extern "C" {
	inline static void pyListToVectorFloat(PyObject* plist, std::vector<float>* pret) {
		const Py_ssize_t s = PyList_Size(plist);
		std::vector<float>& ret = *pret;
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			ret[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(plist, i)));
		}
	}
	inline static void pyListToVectorCurPoint(PyObject* plist, std::vector<PointType>* retp) {
		std::vector<PointType>& ret = *retp;
		const Py_ssize_t s = PyList_Size(plist);
		ret.reserve(s);
		for (Py_ssize_t i = 0; i < s; i += 3) {
			ret.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i))),
							 static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i + 1))),
							 static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(plist, i + 2))));
		}
	}
	inline static void pyListToVectorInt(PyObject* plist, std::vector<int>* pret) {
		const Py_ssize_t s = PyList_Size(plist);
		std::vector<int>& ret = *pret;
		ret.resize(s);
		for (Py_ssize_t i = 0; i < s; i++) {
			ret[i] = static_cast<int>(PyLong_AsLong(PyList_GetItem(plist, i)));
		}
	}
	inline static void pyListToVectorCharP(PyObject* plist, std::vector<const char*>* pret) {

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
	}

	inline static void useDistances(PyObject* self) {
		if (p_distances != nullptr)
			return;
		std::string full(PyUnicode_AsUTF8(PyObject_GetAttrString(self, "__file__")));
		auto found = full.find_last_of("\\/");
		auto bond_filename = full.substr(0, found + 1) + "BondLength.ini";
		static DistancesType dist(bond_filename);
		p_distances = &dist;
	}
	// Python section
	static PyObject* cpplib_GenBonds(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *(p_distances);
		const Py_ssize_t s = PyList_Size(arg);
		std::vector<AtomTypeData> types;
		types.reserve(s);
		std::vector<PointType> points;
		points.reserve(s);

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* tp = PyList_GetItem(arg, i);
			types.emplace_back(static_cast<AtomTypeData>(PyLong_AsLong(PyList_GetItem(tp, 0))));
			points.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 1))),
								static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 2))),
								static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 3))));
		}
		FAMStructType famstr(std::move(types), std::move(points));
		std::string errM;
		auto&& bonds = famstr.findBonds(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		std::string line;
		for (size_t i = 0; i < bonds.size(); i++)
		{
			line += std::to_string(bonds[i].first) + ':' + std::to_string(bonds[i].second) + '\n';
		}
		return PyUnicode_FromString(line.c_str());
	}
	static PyObject* cpplib_GenBondsEx(PyObject* self, PyObject* arg) {
		useDistances(self);
		auto& distances = *(p_distances);
		const Py_ssize_t s = PyList_Size(arg);
		std::vector<AtomTypeData> types;
		types.reserve(s);
		std::vector<PointType> points;
		points.reserve(s);

		for (Py_ssize_t i = 0; i < s; i++) {
			PyObject* tp = PyList_GetItem(arg, i);
			types.emplace_back(static_cast<AtomTypeData>(PyLong_AsLong(PyList_GetItem(tp, 0))));
			points.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 1))),
								static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 2))),
								static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 3))));
		}
		FAMStructType famstr(std::move(types), std::move(points));
		std::string errM;
		auto&& bonds = famstr.findBondsEx(distances, errM, [](const PointType& p1, const PointType& p2) {return (p1 - p2).r(); }).first;

		std::string line;
		for (size_t i = 0; i < bonds.size(); i++)
		{
			line += std::to_string(bonds[i].first) + ':' + std::to_string(bonds[i].second)
				+ ':' + std::to_string(bonds[i].length) + '\n';
		}
		return PyUnicode_FromString(line.c_str());
	}

	static PyObject* cpplib_SearchMain(PyObject* self, PyObject* args) {
		const char* search = NULL;
		PyObject* o = NULL;
		int np = 0;
		int exact = 0;
		PyArg_ParseTuple(args, "sOip", &search, &o, &np, &exact);
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
		deb_write("py_CompareGraph in");
		const char* s1 = NULL;
		const char* s2 = NULL;
		int b = 0;
		deb_write("py_CompareGraph arg parse start");
		PyArg_ParseTuple(args, "ssp", &s1, &s2, &b);
		deb_write("s1 = ", s1);
		deb_write("s2 = ", s2);
		deb_write("exact = ", b);
		deb_write("py_CompareGraph invoke CompareGraph");
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
		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otypes, &oxyz);

		std::array<float, 6> cell;
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}

		std::vector<const char*> symm; pyListToVectorCharP(osymm, &symm);
		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		auto ret = FindMoleculesInCell(cell, symm, types, xyz);
		PyObject* o_xyz_block = PyList_New(0);

		PyObject* o_ret = PyDict_New();

		PyDict_SetItemString(o_ret, "graph_str", PyUnicode_FromString(ret.first.c_str()));

		for (auto & mol : ret.second)
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
		return Py_BuildValue("{s:s,s:O}", 
							 "graph_str", ret.first.c_str(),
							 "xyz_block", o_xyz_block);
	}

	static PyObject* cpplib_FindMoleculesWithoutCell(PyObject* self, PyObject* args) {
		useDistances(self);

		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;

		PyArg_ParseTuple(args, "OO", &otypes, &oxyz);


		Py_ssize_t s = PyList_Size(otypes);
		const int types_s = static_cast<int>(s);

		std::vector<int> types;
		std::vector<float> xyz;
		for (Py_ssize_t i = 0; i < s; i++) {
			types[i] = static_cast<int>(PyLong_AsLong(PyList_GetItem(otypes, i)));
			xyz[i * 3] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3)));
			xyz[i * 3 + 1] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 1)));
			xyz[i * 3 + 2] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 2)));
		}

		auto ret = FindMoleculesWithoutCell(types, xyz);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_GenSymm(PyObject* self, PyObject* args) {
		PyObject* osymm = NULL;
		PyObject* arg = NULL;
		unsigned char flags;

		deb_write("Parsing start");
		if(!PyArg_ParseTuple(args, "ObO", &arg, &flags, &osymm)) {
			deb_write("Parse Error");
			Py_RETURN_NONE;
		}
		bool movetocell = flags & 1;
		bool movemasstocell = flags & 2;

		const Py_ssize_t s_points = PyList_Size(arg);
		deb_write("arg size = ", s_points);
		std::vector<AtomTypeData> types;
		types.reserve(s_points);
		std::vector<PointType> points;
		points.reserve(s_points);

		deb_write("Parsing: atom parsing started");
		for (Py_ssize_t i = 0; i < s_points; i++) {
			PyObject* tp = PyList_GetItem(arg, i);
			types.emplace_back(PyLong_AsLong(PyList_GetItem(tp, 0)));
			points.emplace_back(PyFloat_AsDouble(PyList_GetItem(tp, 1)), PyFloat_AsDouble(PyList_GetItem(tp, 2)), PyFloat_AsDouble(PyList_GetItem(tp, 3)));
		}
		deb_write("pyListToVectorCharP invoking");
		std::vector<const char*> nsymm; pyListToVectorCharP(osymm, &nsymm);

		deb_write("Parsing: symm parsing started");
		std::vector<SymmType> symm;
		const size_t ss = nsymm.size();
		for (size_t i = 0; i < ss; i++)
		{
			symm.emplace_back(nsymm[i]);
		}
		deb_write("Parsing ended");
		if (movemasstocell) {
			deb_write("Move center of mass started");
			PointType centerofmass(0, 0, 0);
			for (size_t i = 0; i < s_points; i++)
			{
				centerofmass += points[i];
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

			deb_write("Move center of mass ended");
		}

		FAMStructType famstr(std::move(types), std::move(points));
		FAMCellType fcell(FAMCellType::base(32, 32, 32, 90, 90, 90, true));
		fcell.GenerateSymm(famstr, symm, movetocell);

		deb_write("famstr.types.size() = ", famstr.types.size());
		deb_write("famstr.points.size() = ", famstr.points.size());
		deb_write("famstr.parseIndex.size() = ", famstr.parseIndex.size());
		deb_write("famstr.sizePoints = ", famstr.sizePoints);

		for (Py_ssize_t i = s_points; i < famstr.sizePoints; i++)
		{
			PyObject* lst = PyList_New(0);
			PyList_Append(lst, PyLong_FromLong(famstr.types[famstr.parseIndex[i]]));
			PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(0)));
			PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(1)));
			PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(2)));
			PyList_Append(arg, lst);
		}
		Py_INCREF(arg);
		return arg;
	}

	static PyObject* cpplib_FindDistanceIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOOOO", &ocell, &osymm, &otypes, &oxyz, &oparams);

		std::array<float, 6> cell;
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}

		std::vector<const char*> symm; pyListToVectorCharP(osymm, &symm);
		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 2> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		std::pair<float, float> value;
		value.first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)));
		value.second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)));


		if (value.first == 0) {
			useDistances(self);
			value.second = p_distances->minDistance(type[0], type[1]);
		}
		if (value.second == 0) {
			useDistances(self);
			value.second = p_distances->maxDistance(type[0], type[1]);
		}


		std::string ret = FindDistanceIC(cell, symm, types, xyz, type, value);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindDistanceWC(PyObject* self, PyObject* args) {

		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOO", &otypes, &oxyz, &oparams);

		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 2> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1)))};
		std::pair<float, float> value;
		value.first = static_cast<float>(PyLong_AsLong(PyList_GetItem(oparams, 2)));
		value.second = static_cast<float>(PyLong_AsLong(PyList_GetItem(oparams, 3)));


		if (value.first == 0) {
			useDistances(self);
			value.second = p_distances->minDistance(type[0], type[1]);
		}
		if (value.second == 0) {
			useDistances(self);
			value.second = p_distances->maxDistance(type[0], type[1]);
		}


		std::string ret = FindDistanceWC(types, xyz, type, value);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindAngleIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOOOO", &ocell, &osymm, &otypes, &oxyz, &oparams);

		std::array<float, 6> cell;
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}

		std::vector<const char*> symm; pyListToVectorCharP(osymm, &symm);
		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 3> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		std::array<std::pair<float, float>, 2> value_d;
		for (char i = 0; i < 2; i++, d += 2)
		{
			value_d[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_d[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

			if (value_d[i].first == 0) {
				useDistances(self);
				value_d[i].second = p_distances->minDistance(type[i], type[i + 1]);
			}

			if (value_d[i].second == 0) {
				useDistances(self);
				value_d[i].second = p_distances->maxDistance(type[i], type[i + 1]);
			}
		}

		std::pair<float, float> value_a;

		value_a.first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
		value_a.second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

		std::string ret = FindAngleIC(cell, symm, types, xyz, type, value_d, value_a);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindAngleWC(PyObject* self, PyObject* args) {

		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOO", &otypes, &oxyz, &oparams);

		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 3> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2)))};
		char d = type.size();
		std::array<std::pair<float, float>, 2> value_d;
		for (char i = 0; i < 2; i++, d += 2)
		{
			value_d[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_d[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

			if (value_d[i].first == 0) {
				useDistances(self);
				value_d[i].second = p_distances->minDistance(type[i], type[i + 1]);
			}

			if (value_d[i].second == 0) {
				useDistances(self);
				value_d[i].second = p_distances->maxDistance(type[i], type[i + 1]);
			}
		}

		std::pair<float, float> value_a;

		value_a.first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
		value_a.second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

		std::string ret = FindAngleWC(types, xyz, type, value_d, value_a);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindTorsionIC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOOOO", &ocell, &osymm, &otypes, &oxyz, &oparams);

		std::array<float, 6> cell;
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}

		std::vector<const char*> symm; pyListToVectorCharP(osymm, &symm);
		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 4> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};
		char d = type.size();
		std::array<std::pair<float, float>, 3> value_d;
		for (char i = 0; i < 3; i++, d += 2)
		{
			value_d[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_d[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

			if (value_d[i].first == 0) {
				useDistances(self);
				value_d[i].second = p_distances->minDistance(type[i], type[i + 1]);
			}

			if (value_d[i].second == 0) {
				useDistances(self);
				value_d[i].second = p_distances->maxDistance(type[i], type[i + 1]);
			}
		}
		std::array<std::pair<float, float>, 2> value_a;
		for (char i = 0; i < 2; i++, d += 2)
		{
			value_a[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_a[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));
		}
		std::pair<float, float> value_t;

		value_t.first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
		value_t.second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

		std::string ret = FindTorsionIC(cell, symm, types, xyz, type, value_d, value_a, value_t);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindTorsionWC(PyObject* self, PyObject* args) {

		PyObject* oxyz = NULL;
		PyObject* otypes = NULL;
		PyObject* oparams = NULL;
		PyArg_ParseTuple(args, "OOO", &otypes, &oxyz, &oparams);

		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		const std::array<int, 4> type {static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 0))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 1))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 2))),
										static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, 3)))};
		char d = type.size();
		std::array<std::pair<float, float>, 3> value_d;
		for (char i = 0; i < 3; i++, d += 2)
		{
			value_d[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_d[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

			if (value_d[i].first == 0) {
				useDistances(self);
				value_d[i].second = p_distances->minDistance(type[i], type[i + 1]);
			}

			if (value_d[i].second == 0) {
				useDistances(self);
				value_d[i].second = p_distances->maxDistance(type[i], type[i + 1]);
			}
		}
		std::array<std::pair<float, float>, 2> value_a;
		for (char i = 0; i < 2; i++, d += 2)
		{
			value_a[i].first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
			value_a[i].second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));
		}
		std::pair<float, float> value_t;

		value_t.first = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d)));
		value_t.second = static_cast<int>(PyLong_AsLong(PyList_GetItem(oparams, d + 1)));

		std::string ret = FindTorsionWC(types, xyz, type, value_d, value_a, value_t);

		return PyUnicode_FromString(ret.c_str());
	}

	static PyObject* cpplib_FindDAT_IC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otypes = NULL;
		PyObject* oxyz = NULL;
		PyArg_ParseTuple(args, "OOOO", &ocell, &osymm, &otypes, &oxyz);

		std::array<float, 6> cell;
		for (Py_ssize_t i = 0; i < 6; i++) {
			cell[i] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(ocell, i)));
		}

		std::vector<const char*> symm; pyListToVectorCharP(osymm, &symm);
		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		auto dat = FindDAT_IC(cell, symm, types, xyz);
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
												static_cast<float>(std::get<3>(datang[i]))));
		}
		PyDict_SetItemString(ret, "angles", list_a);

		PyObject* list_t = PyList_New(0);
		for (size_t i = 0; i < dattor.size(); i++) {
			PyList_Append(list_t, Py_BuildValue("(IIIIf)", 
												static_cast<unsigned int>(std::get<0>(dattor[i])),
												static_cast<unsigned int>(std::get<1>(dattor[i])),
												static_cast<unsigned int>(std::get<2>(dattor[i])),
												static_cast<unsigned int>(std::get<3>(dattor[i])),
												static_cast<float>(std::get<4>(dattor[i]))));
		}
		PyDict_SetItemString(ret, "angles", list_t);

		return ret;
	}
	static PyObject* cpplib_FindDAT_WC(PyObject* self, PyObject* args) {

		PyObject* ocell = NULL;
		PyObject* osymm = NULL;
		PyObject* otypes = NULL;
		PyObject* oxyz = NULL;
		PyArg_ParseTuple(args, "OO", &otypes, &oxyz);

		std::vector<int> types; pyListToVectorInt(otypes, &types);
		std::vector<float> xyz; pyListToVectorFloat(oxyz, &xyz);

		auto dat = FindDAT_WC(types, xyz);
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
												static_cast<float>(std::get<3>(datang[i]))));
		}
		PyDict_SetItemString(ret, "angles", list_a);

		PyObject* list_t = PyList_New(0);
		for (size_t i = 0; i < dattor.size(); i++) {
			PyList_Append(list_t, Py_BuildValue("(IIIIf)", 
												static_cast<unsigned int>(std::get<0>(dattor[i])),
												static_cast<unsigned int>(std::get<1>(dattor[i])),
												static_cast<unsigned int>(std::get<2>(dattor[i])),
												static_cast<unsigned int>(std::get<3>(dattor[i])),
												static_cast<float>(std::get<4>(dattor[i]))));
		}
		PyDict_SetItemString(ret, "angles", list_t);

		return ret;
	}

}

// PyDict_SetItemString


static struct PyMethodDef methods[] = {
	{ "GenBonds", cpplib_GenBonds, METH_O, "Generate bond list"},
	{ "GenBondsEx", cpplib_GenBondsEx, METH_O, "Generate bond list with length"},
	{ "GenSymm", cpplib_GenSymm, METH_VARARGS, "Generates symmetry by symm code"},
	{ "SearchMain", cpplib_SearchMain, METH_VARARGS, "Compare graph with data"},
	{ "CompareGraph", cpplib_CompareGraph, METH_VARARGS, "Compare two graphs"},
	{ "FindMoleculesInCell", cpplib_FindMoleculesInCell, METH_VARARGS, "Create graph from cell"},
	{ "FindMoleculesWithoutCell", cpplib_FindMoleculesWithoutCell, METH_VARARGS, "Create graph from xyz"},
	{ "FindDistanceIC", cpplib_FindDistanceIC, METH_VARARGS, "Find distances with current parameters in cell"},
	{ "FindDistanceWC", cpplib_FindDistanceWC, METH_VARARGS, "Find distances with current parameters in xyz"},
	{ "FindAngleIC", cpplib_FindAngleIC, METH_VARARGS, "Find angles with current parameters in cell"},
	{ "FindAngleWC", cpplib_FindAngleWC, METH_VARARGS, "Find angles with current parameters in xyz"},
	{ "FindTorsionIC", cpplib_FindTorsionIC, METH_VARARGS, "Find torsions with current parameters in cell"},
	{ "FindTorsionWC", cpplib_FindTorsionWC, METH_VARARGS, "Find torsions with current parameters in xyz"},
	{ "FindDAT_IC", cpplib_FindDAT_IC, METH_VARARGS, "Create dictionary with distances, angles and torsions in cell"},
	{ "FindDAT_WC", cpplib_FindDAT_WC, METH_VARARGS, "Create dictionary with distances, angles and torsions in xyz"},

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