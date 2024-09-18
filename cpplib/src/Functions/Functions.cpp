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
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple&& dat, const cpplib::currents::FAMStructType& fs);


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

std::tuple<std::string, std::string, FindMoleculesType::RightType> FindMoleculesInCell(const std::array<float, 6>& unit_cell,
																		 std::vector<const char*>& symm, 
																		 cpplib::currents::FAMStructType::AtomContainerType& types,
																		 cpplib::currents::FAMStructType::PointConteinerType& points) {
	auto& distances = *(p_distances);
	if (p_distances->isReady() == false) {
		return std::make_tuple(std::string(),std::string("Error! Could not open BondLength.ini"),
							  FindMoleculesType::RightType());
	}
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, std::move(types), std::move(points));

	fc.CreateSupercell(fs.points, fs.findCutoff(distances));

	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, [fc](const PointType& p1, const PointType& p2) {return fc.distanceInCell(p1, p2); });
	
	FindMoleculesType fm(std::move(fs));

	auto ret = fm.findMolecules(distances, res.first, res.second, errorMsg);

	// For Petr's enjoyment
	auto mol_s = std::get<2>(ret).size();
	for (size_t i = 0; i < mol_s; i++)
	{
		auto j_s = std::get<2>(std::get<2>(ret)[i]).size();
		for (size_t j = 0; j < j_s; j++) {
			auto ai = std::get<2>(std::get<2>(ret)[i])[j].first;
			auto& a = std::get<0>(std::get<0>(std::get<2>(ret)[i])[ai]);
			auto bi = std::get<2>(std::get<2>(ret)[i])[j].second;
			auto& b = std::get<0>(std::get<0>(std::get<2>(ret)[i])[bi]);
			auto d = a - b;
			auto dd = d.round();
			if (d.r() > 0.5) {
				d += 0;
			}
			b += dd;
		}
	}

	for (size_t i = 0; i < mol_s; i++)
	{
		auto j_s = std::get<0>(std::get<2>(ret)[i]).size();
		for (size_t j = 0; j < j_s; j++)
		{
			auto& point = std::get<0>(std::get<0>(std::get<2>(ret)[i])[j]);
			point = fc.fracToCart() * point;
		}
	}
	return ret;
}
std::tuple<std::string, std::string, FindMoleculesType::RightType>  FindMoleculesWithoutCell(cpplib::currents::FAMStructType::AtomContainerType& types,
																			  cpplib::currents::FAMStructType::PointConteinerType& points) {
	auto& distances = *(p_distances);

	if (p_distances->isReady() == false) {
		return std::make_tuple(std::string(), std::string("Error! Could not open BondLength.ini"),
							   FindMoleculesType::RightType());
	}

	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, PointType::distance);

	FindMoleculesType fm(std::move(fs));
	
	return fm.findMolecules(distances, res.first, res.second, errorMsg);
}

std::string FindDistanceWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						   cpplib::currents::FAMStructType::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<float, float>& value)
{
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
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
						   cpplib::currents::FAMStructType::AtomContainerType& types,
						   cpplib::currents::FAMStructType::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<float, float>& value) {
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, std::move(types), std::move(points));

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

std::string FindAngleWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						cpplib::currents::FAMStructType::PointConteinerType& points,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a)
{
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
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
						cpplib::currents::FAMStructType::AtomContainerType& types,
						cpplib::currents::FAMStructType::PointConteinerType& points,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a)
{
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, std::move(types), std::move(points));

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

std::string FindTorsionWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						  cpplib::currents::FAMStructType::PointConteinerType& points,
							const std::array<int, 4>& type,
							const std::array<std::pair<float, float>, 3>& value_d,
							const std::array<std::pair<float, float>, 2>& value_a,
							const std::pair<float, float>& value_t)
{
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
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
						  cpplib::currents::FAMStructType::AtomContainerType& types,
						  cpplib::currents::FAMStructType::PointConteinerType& points,
						  const std::array<int, 4>& type,
						  const std::array<std::pair<float, float>, 3>& value_d,
						  const std::array<std::pair<float, float>, 2>& value_a,
						  const std::pair<float, float>& value_t) {
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, std::move(types), std::move(points));

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

cpplib::DATTuple FindDAT_IC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							cpplib::currents::FAMStructType::AtomContainerType& types,
							cpplib::currents::FAMStructType::PointConteinerType& points) {
	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	ParseDataType(fs, fc, symm, std::move(types), std::move(points));


	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++)
	{
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometryType fg(fs);

	return ConvertDATTuple(fg.findMolDAT_Rad(*p_distances), fs);
}

cpplib::DATTuple FindDAT_WC(cpplib::currents::FAMStructType::AtomContainerType& types,
							cpplib::currents::FAMStructType::PointConteinerType& points) {
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
	FindGeometryType fg(fs);
	return ConvertDATTuple(fg.findMolDAT_Rad(*p_distances), fs);
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
		graph.setupData(std::move(molData));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}

static void reorder(cpplib::FindGeometry::tupleDistance& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = fs.parseIndex[std::get<0>(d)];
	std::get<1>(d) = fs.parseIndex[std::get<1>(d)];

	if (std::get<0>(d) > std::get<1>(d))
		std::swap(std::get<0>(d), std::get<1>(d));
}
static void reorder(cpplib::FindGeometry::tupleAngle& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = fs.parseIndex[std::get<0>(d)];
	std::get<1>(d) = fs.parseIndex[std::get<1>(d)];
	std::get<2>(d) = fs.parseIndex[std::get<2>(d)];

	if (std::get<0>(d) > std::get<2>(d))
		std::swap(std::get<0>(d), std::get<2>(d));
}
static void reorder(cpplib::FindGeometry::tupleTorsion& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = fs.parseIndex[std::get<0>(d)];
	std::get<1>(d) = fs.parseIndex[std::get<1>(d)];
	std::get<2>(d) = fs.parseIndex[std::get<2>(d)];
	std::get<3>(d) = fs.parseIndex[std::get<3>(d)];

	if ((std::get<0>(d) > std::get<3>(d)) || (std::get<0>(d) == std::get<3>(d) && std::get<1>(d) > std::get<2>(d))) {
		std::swap(std::get<0>(d), std::get<3>(d));
		std::swap(std::get<1>(d), std::get<2>(d));
	}
}
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple&& dat, const cpplib::currents::FAMStructType& fs) {
	auto& dists = std::get<0>(dat);
	auto s_dists = dists.size();
	auto& angles = std::get<1>(dat);
	auto s_angles = angles.size();
	auto& tors = std::get<2>(dat);
	auto s_tors = tors.size();

	for (size_t i = 0; i < s_dists; i++)
	{
		reorder(dists[i], fs);
	}
	for (size_t i = 0; i < s_angles; i++)
	{
		reorder(angles[i], fs);
	}
	for (size_t i = 0; i < s_tors; i++)
	{
		reorder(tors[i], fs);
	}
	// erase dublicates
	std::sort(dists.begin(), dists.end());
	
	for (auto it2 = (++dists.begin()), it = dists.begin(); it2 != dists.end(); it++, it2++)
	{
		if ((std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (::std::abs(std::get<2>(*it) - std::get<2>(*it2)) < 0.0001)) {
			dists.erase(it2);
			it2 = it;
			it2++;
		}
	}
	std::sort(angles.begin(), angles.end());
	for (auto it2 = (++angles.begin()), it = angles.begin(); it2 != angles.end(); it++, it2++)
	{
		if ((std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (std::get<2>(*it) == std::get<2>(*it2)) && (::std::abs(std::get<3>(*it) - std::get<3>(*it2)) < 0.0001)) {
			angles.erase(it2);
			it2 = it;
			it2++;
		}
	}
	std::sort(tors.begin(), tors.end());
	for (auto it2 = (++tors.begin()), it = tors.begin(); it2 != tors.end();)
	{
		if ((std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) &&
			(std::get<2>(*it) == std::get<2>(*it2)) && (std::get<3>(*it) == std::get<3>(*it2)) && (::std::abs(std::get<4>(*it) - std::get<4>(*it2)) < 0.0001)) {
			tors.erase(it2);
			it2 = it;
		}
		else
			it++;
		it2++;
	}

	return dat;
}