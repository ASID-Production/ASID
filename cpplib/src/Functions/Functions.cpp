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
#include "../BaseHeaders/DebugMes.h"
#include "../Classes/Interfaces.h"
#include "../Platform/Definitions.h"
#include "AllInOneAndCurrent.h"

#include <thread>
#include <vector>
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static void ChildThreadFunc(const CurrentRequestGraph& input, const AtomicIDType MaxAtom, SearchDataInterface<MolecularIDType, size_type>& dataInterface, const bool exact);

const CurrentDistances* p_distances = nullptr;

bool CompareGraph(const char* search1, const char* search2, const bool exact) {
	deb_write("CompareGraph start");
	CurrentSearchGraph graph;

	deb_write("CompareGraph CurrentSearchGraph start ReadInput");
	graph.setupInput(CurrentRequestGraph::ReadInput(search1));
	deb_write("CompareGraph CurrentSearchGraph start ReadData");
	graph.setupData(CurrentDatabaseGraph::ReadData(search2));
	deb_write("CompareGraph CurrentSearchGraph start prepareSearch");
	graph.prepareToSearch();
	deb_write("CompareGraph CurrentSearchGraph start FullSearch");
	return graph.startFullSearch(exact);
}
int* SearchMain(const char* search, const char** data, const int data_s, const int np, const bool exact) {
	static std::vector<int> result;

	CurrentRequestGraph input = CurrentRequestGraph::ReadInput(search);
	SearchDataInterface<MolecularIDType, size_type> databuf(data, data_s);
	std::vector<std::thread> threads;
	const size_t nThreads = std::min(std::min(static_cast<unsigned int>(np), std::thread::hardware_concurrency()), static_cast<unsigned int>(data_s)) - 1;
	threads.reserve(nThreads);
	auto ma = input.findStart();

	for (size_t i = 0; i < nThreads; i++) {
		threads.emplace_back(ChildThreadFunc, std::cref(input), ma, std::ref(databuf), exact);
	}
	ChildThreadFunc(input, ma, databuf, exact);

	for (size_t i = 0; i < nThreads; i++) {
		threads[i].join();
	}

	result = databuf.getAllResults();
	return &(result[0]);
}

const char* FindMoleculesInCell(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s) {

	static std::string ret;
	auto& distances = *(p_distances);
	if (p_distances->isReady() == false) {
		ret = std::string(";Error! Could not open BondLength.ini");
		return ret.c_str();
	}

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr;
	ParseData(famstr, types, xyz, types_s);

	FAM_Cell<FloatingPointType> fc(CurrentCell(unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], true));

	std::vector<geometry::Symm<FloatingPointType>> symmv;
	symmv.reserve(symm_s - 1);
	for (int i = 1; i < symm_s; i++) {
		symmv.emplace_back(symm[i]);
	}

	fc.GenerateSymm(famstr, symmv);
	famstr.sizePoints = static_cast<AtomicIDType>(famstr.points.size());
	famstr.types.reserve(famstr.sizePoints);
	for (size_type i = famstr.sizeUnique; i < famstr.sizePoints; i++) {
		famstr.types.emplace_back(famstr.types[famstr.parseIndex[i]]);
	}
	fc.CreateSupercell(famstr.points, famstr.findCutoff(distances));
	std::string errorMsg;
	auto res = famstr.findBonds(distances, errorMsg, [fc](const CurrentPoint& p1, const CurrentPoint& p2) {return fc.distanceInCell(p1, p2); });

	CurrentFindMolecules fm(std::move(famstr));

	ret = fm.findMolecules(distances, res.first, res.second, errorMsg);

	return ret.c_str();
}
const char* FindMoleculesWithoutCell(const int* types, const float* xyz, const int types_s) {
	static std::string ret;
	auto& distances = *(p_distances);

	if (distances.isReady() == false) {
		ret = std::string(";Error! Could not open BondLength.ini");
		return ret.c_str();
	}

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr;
	ParseData(famstr, types, xyz, types_s);
	std::string errorMsg;
	auto res = famstr.findBonds(distances, errorMsg, [](const CurrentPoint& p1, const CurrentPoint& p2) {return (p1 - p2).r(); });

	CurrentFindMolecules fm(std::move(famstr));

	ret = fm.findMolecules(distances, res.first, res.second, errorMsg);

	return ret.c_str();
}
const char* FindDistanceWC(const int* types, const float* xyz, const int types_s,const int type1, const int type2, const float min_value, const float max_value) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> fs;
	ParseData(fs, types, xyz, types_s);
	CurrentFindGeometry fg(fs);
	const auto raw = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min_value, max_value);
	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ";\n";
	}
	return res.c_str();
}
const char* FindDistanceIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr;
	ParseData(famstr, types, xyz, types_s);

	FAM_Cell<FloatingPointType> fc(CurrentCell(unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], true));

	std::vector<geometry::Symm<FloatingPointType>> symmv;
	symmv.reserve(symm_s - 1);
	for (int i = 1; i < symm_s; i++) {
		symmv.emplace_back(symm[i]);
	}

	fc.GenerateSymm(famstr, symmv);
	famstr.sizePoints = static_cast<AtomicIDType>(famstr.points.size());
	famstr.types.reserve(famstr.sizePoints);
	for (size_type i = famstr.sizeUnique; i < famstr.sizePoints; i++) {
		famstr.types.emplace_back(famstr.types[famstr.parseIndex[i]]);
	}
	fc.CreateSupercell(famstr.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < famstr.sizePoints; i++)
	{
		famstr.points[i] = fc.fracToCart() * famstr.points[i];
	}

	CurrentFindGeometry fg(famstr);
	const auto raw = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min_value, max_value);
	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += (std::to_string(std::get<2>(raw[i])) + ';') + '\n';
	}
	return res.c_str();
}

const char* FindAngleWC(const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> fs;
	ParseData(fs, types, xyz, types_s);
	CurrentFindGeometry fg(fs);
	const auto raw12 = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min12, max12);
	const auto raw23 = fg.findDistance(static_cast<AtomType>(type2), static_cast<AtomType>(type3), min23, max23);
	const auto raw = fg.findAngle(raw12, raw23, geometry::GradtoRad(min123), geometry::GradtoRad(max123));
	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += (std::to_string(geometry::RadtoGrad(std::get<3>(raw[i]))) + ';') + '\n';
	}
	return res.c_str();
}
const char* FindAngleIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, 
						const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr;
	ParseData(famstr, types, xyz, types_s);

	FAM_Cell<FloatingPointType> fc(CurrentCell(unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], true));

	std::vector<geometry::Symm<FloatingPointType>> symmv;
	symmv.reserve(symm_s - 1);
	for (int i = 1; i < symm_s; i++) {
		symmv.emplace_back(symm[i]);
	}

	fc.GenerateSymm(famstr, symmv);
	famstr.sizePoints = static_cast<AtomicIDType>(famstr.points.size());
	famstr.types.reserve(famstr.sizePoints);
	for (size_type i = famstr.sizeUnique; i < famstr.sizePoints; i++) {
		famstr.types.emplace_back(famstr.types[famstr.parseIndex[i]]);
	}
	fc.CreateSupercell(famstr.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < famstr.sizePoints; i++)
	{
		famstr.points[i] = fc.fracToCart() * famstr.points[i];
	}

	CurrentFindGeometry fg(famstr);

	const auto raw12 = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min12, max12);
	const auto raw23 = fg.findDistance(static_cast<AtomType>(type2), static_cast<AtomType>(type3), min23, max23);
	const auto raw = fg.findAngle(raw12, raw23, geometry::GradtoRad(min123), geometry::GradtoRad(max123));
	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += (std::to_string(geometry::RadtoGrad(std::get<3>(raw[i]))) + ';') + '\n';
	}
	return res.c_str();
}

const char* FindTorsionWC(const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
	const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> fs;
	ParseData(fs, types, xyz, types_s);
	CurrentFindGeometry fg(fs);
	const auto raw12 = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min12, max12);
	const auto raw23 = fg.findDistance(static_cast<AtomType>(type2), static_cast<AtomType>(type3), min23, max23);
	const auto raw34 = fg.findDistance(static_cast<AtomType>(type3), static_cast<AtomType>(type4), min34, max34);
	const auto raw123 = fg.findAngle(raw12, raw23, geometry::GradtoRad(min123), geometry::GradtoRad(max123));
	const auto raw234 = fg.findAngle(raw23, raw34, geometry::GradtoRad(min234), geometry::GradtoRad(max234));
	const auto raw = fg.findTorsion(raw123, raw234, geometry::GradtoRad(min1234), geometry::GradtoRad(max1234));

	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += std::to_string(std::get<3>(raw[i])) + ':';
		res += std::to_string(geometry::RadtoGrad(std::get<4>(raw[i]))) + ";\n";
	}
	return res.c_str();
}
const char* FindTorsionIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
	const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234) {
	static std::string res;

	FAM_Struct<AtomType, AtomicIDType, FloatingPointType> famstr;
	ParseData(famstr, types, xyz, types_s);

	FAM_Cell<FloatingPointType> fc(CurrentCell(unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], true));

	std::vector<geometry::Symm<FloatingPointType>> symmv;
	symmv.reserve(symm_s - 1);
	for (int i = 1; i < symm_s; i++) {
		symmv.emplace_back(symm[i]);
	}

	fc.GenerateSymm(famstr, symmv);
	famstr.sizePoints = static_cast<AtomicIDType>(famstr.points.size());
	famstr.types.reserve(famstr.sizePoints);
	for (size_type i = famstr.sizeUnique; i < famstr.sizePoints; i++) {
		famstr.types.emplace_back(famstr.types[famstr.parseIndex[i]]);
	}
	fc.CreateSupercell(famstr.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < famstr.sizePoints; i++)
	{
		famstr.points[i] = fc.fracToCart() * famstr.points[i];
	}

	CurrentFindGeometry fg(famstr);

	const auto raw12 = fg.findDistance(static_cast<AtomType>(type1), static_cast<AtomType>(type2), min12, max12);
	const auto raw23 = fg.findDistance(static_cast<AtomType>(type2), static_cast<AtomType>(type3), min23, max23);
	const auto raw34 = fg.findDistance(static_cast<AtomType>(type3), static_cast<AtomType>(type4), min34, max34);
	const auto raw123 = fg.findAngle(raw12, raw23, geometry::GradtoRad(min123), geometry::GradtoRad(max123));
	const auto raw234 = fg.findAngle(raw23, raw34, geometry::GradtoRad(min234), geometry::GradtoRad(max234));
	const auto raw = fg.findTorsion(raw123, raw234, geometry::GradtoRad(min1234), geometry::GradtoRad(max1234));

	const auto raws = raw.size();
	res.clear();
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += std::to_string(std::get<3>(raw[i])) + ':';
		res += (std::to_string(geometry::RadtoGrad(std::get<4>(raw[i]))) + ';') + '\n';
	}
	return res.c_str();
}

// Single thread function
static void ChildThreadFunc(const CurrentRequestGraph& input, const AtomicIDType MaxAtom, SearchDataInterface<MolecularIDType, size_type>& dataInterface, const bool exact) {
	CurrentSearchGraph graph;
	while (true) {
		auto next = dataInterface.getNext();
		if (next == nullptr) {
			return;
		}
		graph.setupInput(input.makeCopy());
		CurrentDatabaseGraph molData = CurrentDatabaseGraph::ReadData(next);
		auto id = molData.getID();
		graph.setupData(move(molData));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}

inline static void useDistances (PyObject* self) {
	if (p_distances != nullptr)
		return;
	std::string full(PyUnicode_AsUTF8(PyObject_GetAttrString(self, "__file__")));
	auto found = full.find_last_of("\\/");
	auto bond_filename = full.substr(0, found+1)+"BondLength.ini";
	static CurrentDistances dist(bond_filename);
	p_distances = &dist;
}
// [[type,x,y,z],[...]...] in, str out
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

// SearchMain(
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
	int* ret = SearchMain(search, data.data(), ds, np, (exact != 0));
	PyObject* ret_o = PyList_New(ret[0]);
	for (Py_ssize_t i = 1; i <= ret[0]; i++)
	{
		PyList_SetItem(ret_o, i - 1, PyLong_FromLong(ret[i]));
	}
	return ret_o;
}

static PyObject * cpplib_CompareGraph(PyObject * self, PyObject * args)
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

	float cell[6];
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

	std::vector<int> types(types_s);
	std::vector<float> xyz(types_s*3);
	for (Py_ssize_t i = 0; i < s; i++) {
		types[i] = PyLong_AsLong(PyList_GetItem(otypes, i));
		xyz[i * 3] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3));
		xyz[i * 3 + 1] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 1));
		xyz[i * 3 + 2] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 2));
	}

	const char* ret = FindMoleculesInCell(cell, &(symm[0]), symm_s, &(types[0]), &(xyz[0]), types_s);

	return PyUnicode_FromString(ret);
}

static PyObject* cpplib_FindMoleculesWithoutCell(PyObject* self, PyObject* args) {
	useDistances(self);

	PyObject* oxyz = NULL;
	PyObject* otypes = NULL;

	PyArg_ParseTuple(args, "OO", &otypes, &oxyz);


	Py_ssize_t s = PyList_Size(otypes);
	const int types_s = static_cast<int>(s);

	std::vector<int> types(types_s);
	std::vector<float> xyz(types_s);
	for (Py_ssize_t i = 0; i < s; i++) {
		types[i] = PyLong_AsLong(PyList_GetItem(otypes, i));
		xyz[i * 3] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3));
		xyz[i * 3 + 1] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 1));
		xyz[i * 3 + 2] = PyFloat_AsDouble(PyList_GetItem(oxyz, i * 3 + 2));
	}

	const char* ret = FindMoleculesWithoutCell( &(types[0]), &(xyz[0]), types_s);

	return PyUnicode_FromString(ret);
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
	symm.emplace_back(str,false);

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