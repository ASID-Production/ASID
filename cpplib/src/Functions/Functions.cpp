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
#include "Functions.h"
#include "AllInOneAndCurrent.h"

#include <thread>
#include <vector>

using namespace cpplib::currents;

static void ChildThreadFunc(const SearchGraphType::RequestGraphType& input, const SearchGraphType::AtomIndex MaxAtom, SearchDataInterfaceType& dataInterface, const bool exact);

const DistancesType* p_distances = nullptr;

bool CompareGraph(const char* search1, const char* search2, const bool exact) {
	deb_write("CompareGraph start");
	SearchGraphType graph;

	deb_write("CompareGraph CurrentSearchGraph start ReadInput");
	auto&& inputpair = SearchGraphType::RequestGraphType::ReadInput(search1);
	graph.setupInput(std::move(inputpair.first));
	deb_write("CompareGraph CurrentSearchGraph start ReadData");
	graph.setupData(SearchGraphType::DatabaseGraphType::ReadData(search2, inputpair.second));
	deb_write("CompareGraph CurrentSearchGraph start prepareSearch");
	graph.prepareToSearch();
	deb_write("CompareGraph CurrentSearchGraph start FullSearch");
	return graph.startFullSearch(exact);
}
std::vector<int> SearchMain(const char* search, std::vector<const char*>&& data, const int np, const bool exact) {
	auto&& inputpair = SearchGraphType::RequestGraphType::ReadInput(search);
	SearchDataInterfaceType databuf(std::move(data), std::move(inputpair.second));
	std::vector<std::thread> threads;
	const size_t nThreads = std::min(std::min(static_cast<unsigned int>(np), std::thread::hardware_concurrency()),
										static_cast<unsigned int>(databuf.size())) - 1;
	threads.reserve(nThreads);
	auto ma = inputpair.first.findStart();

	for (size_t i = 0; i < nThreads; i++) {
		threads.emplace_back(ChildThreadFunc, std::cref(inputpair.first), ma, std::ref(databuf), exact);
	}
	ChildThreadFunc(inputpair.first, ma, databuf, exact);

	for (size_t i = 0; i < nThreads; i++) {
		threads[i].join();
	}

	return databuf.getAllResults();
}

std::pair<std::string, std::vector<std::pair<std::vector<std::tuple<PointType, AtomIndex>>, int>>>
	FindMoleculesInCell(const std::array<float, 6>& unit_cell, std::vector<const char*>& symm, std::vector<int>& types, std::vector<float>& xyz) {

	using return_type = std::pair<std::string, std::vector<std::pair<std::vector<std::tuple<PointType, AtomIndex>>, int>>>;

	auto& distances = *(p_distances);
	if (p_distances->isReady() == false) {
		return std::make_pair(std::string(";Error! Could not open BondLength.ini"),
								return_type::second_type());
	}
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, types, xyz);

	fc.CreateSupercell(fs.points, fs.findCutoff(distances));

	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, [fc](const PointType& p1, const PointType& p2) {return fc.distanceInCell(p1, p2); });

	for (size_t i = 0; i < fs.sizePoints; i++)
	{
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindMoleculesType fm(std::move(fs));

	return fm.findMolecules(distances, res.first, res.second, errorMsg);
}
std::string FindMoleculesWithoutCell(const std::vector<int>& types, std::vector<float>& xyz) {
	std::string ret;
	auto& distances = *(p_distances);

	if (distances.isReady() == false) {
		return std::string(";Error! Could not open BondLength.ini");
	}

	FAMStructType fs;
	ParseDataType(fs, types, xyz);
	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, PointType::distance);

	FindMoleculesType fm(std::move(fs));

	ret = fm.findMolecules(distances, res.first, res.second, errorMsg).first;

	return ret.c_str();
}

std::string FindDistanceWC(std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 2>& type,
							const std::pair<float, float>& value)
{
	FAMStructType fs;
	ParseDataType(fs, types, xyz);
	FindGeometryType fg(fs);
	const auto raw = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[0]),
										static_cast<FindGeometryType::AtomType>(type[1]),
										value);
	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ";\n";
	}
	return res;
}

std::string FindDistanceIC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 2>& type,
							const std::pair<float, float>& value) {
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, types, xyz);

	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++)
	{
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometryType fg(fs);
	const auto raw = fg.findDistance(static_cast<FindGeometryType::AtomType>(std::get<0>(type)),
										static_cast<FindGeometryType::AtomType>(std::get<1>(type)),
										value);
	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += (std::to_string(std::get<2>(raw[i])) + ';') + '\n';
	}
	return res;
}

std::string FindAngleWC(std::vector<int>& types,
						std::vector<float>& xyz,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a)
{
	FAMStructType fs;
	ParseDataType(fs, types, xyz);
	FindGeometryType fg(fs);
	const auto raw12 = fg.findDistance(static_cast<FindGeometryType::AtomType>(std::get<0>(type)),
										static_cast<FindGeometryType::AtomType>(std::get<1>(type)),
										value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometryType::AtomType>(std::get<1>(type)),
										static_cast<FindGeometryType::AtomType>(std::get<2>(type)),
										value_d[1]);
	const auto raw = fg.findAngle(raw12, raw23, std::make_pair(cpplib::geometry::GradtoRad(value_a.first), cpplib::geometry::GradtoRad(value_a.second)));
	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += (std::to_string(cpplib::geometry::RadtoGrad(std::get<3>(raw[i]))) + ';') + '\n';
	}
	return res;
}

std::string FindAngleIC(const std::array<float, 6>& unit_cell,
						std::vector<const char*>& symm,
						std::vector<int>& types,
						std::vector<float>& xyz,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a)
{
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, types, xyz);

	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++)
	{
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometryType fg(fs);
	const auto raw12 = fg.findDistance(static_cast<FindGeometryType::AtomType>(std::get<0>(type)),
										static_cast<FindGeometryType::AtomType>(std::get<1>(type)),
										value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometryType::AtomType>(std::get<1>(type)),
										static_cast<FindGeometryType::AtomType>(std::get<2>(type)),
										value_d[1]);
	const auto raw = fg.findAngle(raw12, raw23, std::make_pair(cpplib::geometry::GradtoRad(value_a.first), cpplib::geometry::GradtoRad(value_a.second)));
	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += (std::to_string(cpplib::geometry::RadtoGrad(std::get<3>(raw[i]))) + ';') + '\n';
	}
	return res;
}

std::string FindTorsionWC(std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 4>& type,
							const std::array<std::pair<float, float>, 3>& value_d,
							const std::array<std::pair<float, float>, 2>& value_a,
							const std::pair<float, float>& value_t)
{
	FAMStructType fs;
	ParseDataType(fs, types, xyz);
	FindGeometryType fg(fs);

	const auto raw12 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[0]), static_cast<FindGeometryType::AtomType>(type[1]), value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[1]), static_cast<FindGeometryType::AtomType>(type[2]), value_d[1]);
	const auto raw34 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[2]), static_cast<FindGeometryType::AtomType>(type[3]), value_d[2]);
	const auto raw123 = fg.findAngle(raw12, raw23,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[0].first), cpplib::geometry::GradtoRad(value_a[0].second)));
	const auto raw234 = fg.findAngle(raw23, raw34,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[1].first), cpplib::geometry::GradtoRad(value_a[1].second)));
	const auto raw = fg.findTorsion(raw123, raw234,
									std::make_pair(cpplib::geometry::GradtoRad(value_t.first), cpplib::geometry::GradtoRad(value_t.second)));

	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += std::to_string(std::get<3>(raw[i])) + ':';
		res += std::to_string(cpplib::geometry::RadtoGrad(std::get<4>(raw[i]))) + ";\n";
	}
	return res;
}
std::string FindTorsionIC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 4>& type,
							const std::array<std::pair<float, float>, 3>& value_d,
							const std::array<std::pair<float, float>, 2>& value_a,
							const std::pair<float, float>& value_t) {
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, types, xyz);

	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++)
	{
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometryType fg(fs);

	const auto raw12 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[0]), static_cast<FindGeometryType::AtomType>(type[1]), value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[1]), static_cast<FindGeometryType::AtomType>(type[2]), value_d[1]);
	const auto raw34 = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[2]), static_cast<FindGeometryType::AtomType>(type[3]), value_d[2]);
	const auto raw123 = fg.findAngle(raw12, raw23,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[0].first), cpplib::geometry::GradtoRad(value_a[0].second)));
	const auto raw234 = fg.findAngle(raw23, raw34,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[1].first), cpplib::geometry::GradtoRad(value_a[1].second)));
	const auto raw = fg.findTorsion(raw123, raw234,
									std::make_pair(cpplib::geometry::GradtoRad(value_t.first), cpplib::geometry::GradtoRad(value_t.second)));

	const auto raws = raw.size();
	std::string res;
	for (size_t i = 0; i < raws; i++)
	{
		res += std::to_string(std::get<0>(raw[i])) + ':';
		res += std::to_string(std::get<1>(raw[i])) + ':';
		res += std::to_string(std::get<2>(raw[i])) + ':';
		res += std::to_string(std::get<3>(raw[i])) + ':';
		res += std::to_string(cpplib::geometry::RadtoGrad(std::get<4>(raw[i]))) + ";\n";
	}
	return res;
}

// Single thread function
static void ChildThreadFunc(const SearchGraphType::RequestGraphType& input, const SearchGraphType::AtomIndex MaxAtom, SearchDataInterfaceType& dataInterface, const bool exact) {
	SearchGraphType graph;
	while (true) {
		auto next = dataInterface.getNext();
		if (next == nullptr) {
			return;
		}
		graph.setupInput(input.makeCopy());
		SearchGraphType::DatabaseGraphType molData = SearchGraphType::DatabaseGraphType::ReadData(next, dataInterface.getMulty());
		auto id = molData.getID();
		graph.setupData(move(molData));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}