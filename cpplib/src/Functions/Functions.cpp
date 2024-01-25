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
#include "../Classes/Interfaces.h"
#include "../Platform/Definitions.h"
#include "AllInOneAndCurrent.h"

#include <thread>
#include <vector>

static void ChildThreadFunc(const CurrentMoleculeGraph& input, const AtomicIDType MaxAtom, SearchDataInterface<MolecularIDType, size_type>& dataInterface, const bool exact);

static const CurrentDistances distances("./modules/c_modules/BondLength.ini");

API bool CompareGraph(const char* search1, const char* search2, const bool exact) {
	CurrentSearchGraph graph;
	graph.setupInput(CurrentMoleculeGraph(search1));
	graph.setupData(CurrentMoleculeGraph(search2));
	graph.prepareToSearch();
	return graph.startFullSearch(exact);
}

API int* SearchMain(const char* search, const char** data, const int data_s, const int np, const bool exact) {
	static std::vector<int> result;

	CurrentMoleculeGraph input(search);
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

API const char* FindMoleculesInCell(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s) {
	static std::string ret;

	if (distances.isReady() == false) {
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
	famstr.sizePoints = famstr.points.size();
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
API const char* FindMoleculesWithoutCell(const int* types, const float* xyz, const int types_s) {
	static std::string ret;

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
API const char* FindDistanceWC(const int* types, const float* xyz, const int types_s,const int type1, const int type2, const float min_value, const float max_value) {
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
API const char* FindDistanceIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value) {
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
	famstr.sizePoints = famstr.points.size();
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


API const char* FindAngleWC(const int* types, const float* xyz, const int types_s, 
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
API const char* FindAngleIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, 
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
	famstr.sizePoints = famstr.points.size();
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

API const char* FindTorsionWC(const int* types, const float* xyz, const int types_s,
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
API const char* FindTorsionIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
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
	famstr.sizePoints = famstr.points.size();
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
static void ChildThreadFunc(const CurrentMoleculeGraph& input, const AtomicIDType MaxAtom, SearchDataInterface<MolecularIDType, size_type>& dataInterface, const bool exact) {
	CurrentSearchGraph graph;
	while (true) {
		auto next = dataInterface.getNext();
		if (next == nullptr) {
			return;
		}
		graph.setupInput(input.makeCopy());
		CurrentMoleculeGraph molData(next);
		auto id = molData.getID();
		graph.setupData(move(molData));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}
