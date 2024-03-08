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
#include "../BaseHeaders/DebugMes.h"

inline static void useDistances(PyObject* self) {
	if (p_distances != nullptr)
		return;
	std::string full(PyUnicode_AsUTF8(PyObject_GetAttrString(self, "__file__")));
	auto found = full.find_last_of("\\/");
	auto bond_filename = full.substr(0, found + 1) + "BondLength.ini";
	static CurrentDistances dist(bond_filename);
	p_distances = &dist;
}

static PyObject* cpplib_GenBonds(PyObject* self, PyObject* arg) {
	useDistances(self);
	auto& distances = *(p_distances);
	const Py_ssize_t s = PyList_Size(arg);
	std::vector<AtomType> types;
	types.reserve(s);
	std::vector<CurrentPoint> points;
	points.reserve(s);

	for (Py_ssize_t i = 0; i < s; i++) {
		PyObject* tp = PyList_GetItem(arg, i);
		types.emplace_back(static_cast<AtomType>(PyLong_AsLong(PyList_GetItem(tp, 0))));
		points.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 1))),
			static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 2))),
			static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 3))));
	}
	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr(std::move(types), std::move(points));
	std::string errM;
	auto&& bonds = famstr.findBonds(distances, errM, [](const CurrentPoint& p1, const CurrentPoint& p2) {return (p1 - p2).r(); }).first;

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
	std::vector<AtomType> types;
	types.reserve(s);
	std::vector<CurrentPoint> points;
	points.reserve(s);

	for (Py_ssize_t i = 0; i < s; i++) {
		PyObject* tp = PyList_GetItem(arg, i);
		types.emplace_back(static_cast<AtomType>(PyLong_AsLong(PyList_GetItem(tp, 0))));
		points.emplace_back(static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 1))),
			static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 2))),
			static_cast<FloatingPointType>(PyFloat_AsDouble(PyList_GetItem(tp, 3))));
	}
	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr(std::move(types), std::move(points));
	std::string errM;
	auto&& bonds = famstr.findBondsEx(distances, errM, [](const CurrentPoint& p1, const CurrentPoint& p2) {return (p1 - p2).r(); }).first;

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
	std::vector<const char*> data(ds, nullptr);

	deb_write("search = ", search);
	deb_write("np = ", np);
	deb_write("exact = ", exact);
	deb_write("py_CompareGraph invoke CompareGraph");
	for (Py_ssize_t i = 0; i < s; i++) {
		data[i] = PyUnicode_AsUTF8(PyList_GetItem(o, i));
	}
	const auto ret = SearchMain(search, std::move(data), np, (exact != 0));
	const auto ret_s = ret.size();
	PyObject* ret_o = PyList_New(0);
	for (Py_ssize_t i = 0; i <= ret[0]; i++)
	{
		PyList_Append(ret_o, PyLong_FromLong(ret[i]));
	}
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

	Py_ssize_t s = PyList_Size(osymm);
	const int symm_s = static_cast<int>(s);

	std::vector<const char*> symm(symm_s, nullptr);
	for (Py_ssize_t i = 0; i < s; i++) {
		symm[i] = PyUnicode_AsUTF8(PyList_GetItem(osymm, i));
	}

	s = PyList_Size(otypes);
	const int types_s = static_cast<int>(s);

	std::vector<int> types;
	std::vector<float> xyz;
	for (Py_ssize_t i = 0; i < s; i++) {
		types[i] = static_cast<int>(PyLong_AsLong(PyList_GetItem(otypes, i)));
		xyz[i * 3] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3)));
		xyz[i * 3 + 1] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 1)));
		xyz[i * 3 + 2] = static_cast<float>(PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 2)));
	}

	auto ret = FindMoleculesInCell(cell, symm, types, xyz);

	return PyUnicode_FromString(ret.c_str());
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

// [[type,x,y,z],...], str
static PyObject* cpplib_GenSymm(PyObject* self, PyObject* args) {
	char* str = NULL;
	PyObject* arg = NULL;

	PyArg_ParseTuple(args, "Os", &arg, &str);

	const Py_ssize_t s = PyList_Size(arg);
	std::vector<AtomType> types;
	types.reserve(s);
	std::vector<CurrentPoint> points;
	points.reserve(s);

	for (Py_ssize_t i = 0; i < s; i++) {
		PyObject* tp = PyList_GetItem(arg, i);
		types.emplace_back(PyLong_AsLong(PyList_GetItem(tp, 0)));
		points.emplace_back(PyFloat_AsDouble(PyList_GetItem(tp, 1)), PyFloat_AsDouble(PyList_GetItem(tp, 2)), PyFloat_AsDouble(PyList_GetItem(tp, 3)));
	}

	std::vector<geometry::Symm<FloatingPointType>> symm;
	symm.emplace_back(str, false);

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr(std::move(types), std::move(points));
	FAM_Cell<FloatingPointType> fcell(CurrentCell(10, 10, 10, 90, 90, 90, true));
	fcell.GenerateSymm(famstr, symm);


	std::string line;
	for (Py_ssize_t i = s; i < famstr.sizePoints; i++)
	{
		PyObject* lst = PyList_New(0);
		PyList_Append(lst, PyLong_FromLong(famstr.types[famstr.parseIndex[i]]));
		PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(0)));
		PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(1)));
		PyList_Append(lst, PyFloat_FromDouble(famstr.points[i].get(2)));
		PyList_Append(arg, lst);
	}
	return PyUnicode_FromString(line.c_str());
}




static struct PyMethodDef methods[] = {
	{ "GenBonds", cpplib_GenBonds, METH_O, "Generate bond list"},
	{ "GenBondsEx", cpplib_GenBondsEx, METH_O, "Generate bond list with length"},
	{ "GenSymm", cpplib_GenSymm, METH_VARARGS, "Generates symmetry by symm code"},
	{ "SearchMain", cpplib_SearchMain, METH_VARARGS, "Compare graph with data"},
	{ "CompareGraph", cpplib_CompareGraph, METH_VARARGS, "Compare two graphs"},
	{ "FindMoleculesInCell", cpplib_FindMoleculesInCell, METH_VARARGS, "Create graph from cell"},
	{ "FindMoleculesWithoutCell", cpplib_FindMoleculesWithoutCell, METH_VARARGS, "Create graph from xyz"},

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