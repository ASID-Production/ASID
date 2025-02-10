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
#include <map>
#include <set>
#include <functional>

using namespace cpplib::currents;

static void ChildThreadFunc(const SearchGraphType::RequestGraphType& input, const SearchGraphType::AtomIndex MaxAtom, SearchDataInterfaceType& dataInterface, const bool exact); 
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple& dat, const cpplib::currents::FAMStructType& fs);


const DistancesType* p_distances = nullptr;

bool CompareGraph(const char* search1, const char* search2, const bool exact) {
	deb_write("CompareGraph start");
	SearchGraphType graph;

	deb_write("CompareGraph CurrentSearchGraph start ReadInput");
	auto&& inputpair = SearchGraphType::RequestGraphType::ReadInput(search1);
	auto map = inputpair.first.getTypeMap();
	graph.setupInput(std::move(inputpair.first));
	deb_write("CompareGraph CurrentSearchGraph start ReadData");
	auto&& d_pair = SearchGraphType::DatabaseGraphType::ReadData(search2, inputpair.second, map);
	if (!d_pair.second) return false;
	graph.setupData(std::move(d_pair.first));
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

std::tuple<std::string, std::string, FindMoleculesType::RightType> FindMoleculesInCell(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
																		 std::vector<const char*>& symm, 
																		 cpplib::currents::FAMStructType::AtomContainerType& types,
																		 cpplib::currents::FAMStructType::PointConteinerType& points) {
	auto& distances = *p_distances;
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
	auto& distances = *p_distances;

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

std::vector<FindGeometryType::tupleDistance> FindDistanceWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						   cpplib::currents::FAMStructType::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value) {
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
	FindGeometryType fg(fs);
	const auto raw = fg.findDistance(static_cast<FindGeometryType::AtomType>(type[0]),
										static_cast<FindGeometryType::AtomType>(type[1]),
										value);
	return raw;
}

std::vector<FindGeometryType::tupleDistance> FindDistanceIC(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
						   std::vector<const char*>& symm,
						   cpplib::currents::FAMStructType::AtomContainerType& types,
						   cpplib::currents::FAMStructType::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value) {
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
	return raw;
}

std::vector<FindGeometryType::tupleAngle> FindAngleWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						                              cpplib::currents::FAMStructType::PointConteinerType& points,
						                              const std::array<int, 3>& type,
						                              const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 2>& value_d,
						                              const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value_a) {
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
	return raw;
}

std::vector<FindGeometryType::tupleAngle> FindAngleIC(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
						                              std::vector<const char*>& symm, 
						                              cpplib::currents::FAMStructType::AtomContainerType& types,
						                              cpplib::currents::FAMStructType::PointConteinerType& points,
						                              const std::array<int, 3>& type,
						                              const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 2>& value_d,
						                              const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value_a)
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
	auto raw = fg.findAngle(raw12, raw23, std::make_pair(cpplib::geometry::GradtoRad(value_a.first), cpplib::geometry::GradtoRad(value_a.second)));
	return raw;
}

std::vector<FindGeometryType::tupleTorsion> FindTorsionWC(cpplib::currents::FAMStructType::AtomContainerType& types,
						                                  cpplib::currents::FAMStructType::PointConteinerType& points,
							                              const std::array<int, 4>& type,
							                              const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 3>& value_d,
							                              const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 2>& value_a,
							                              const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value_t)
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
	auto raw = fg.findTorsion(raw123, raw234,
							  std::make_pair(cpplib::geometry::GradtoRad(value_t.first), cpplib::geometry::GradtoRad(value_t.second)));
	return raw;
}
std::vector<FindGeometryType::tupleTorsion> FindTorsionIC(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
						                                  std::vector<const char*>& symm,
						                                  cpplib::currents::FAMStructType::AtomContainerType& types,
						                                  cpplib::currents::FAMStructType::PointConteinerType& points,
						                                  const std::array<int, 4>& type,
						                                  const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 3>& value_d,
						                                  const std::array<std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>, 2>& value_a,
						                                  const std::pair<cpplib::currents::FloatingPointType, cpplib::currents::FloatingPointType>& value_t) {
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
	return raw;
}

cpplib::DATTuple FindDAT_IC(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
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
	auto Moldat = fg.findMolDAT_Rad(*p_distances);
	return ConvertDATTuple(Moldat, fs);
}

cpplib::DATTuple FindDAT_WC(cpplib::currents::FAMStructType::AtomContainerType& types,
							cpplib::currents::FAMStructType::PointConteinerType& points) {
	FAMStructType fs;
	ParseDataType(fs, std::move(types), std::move(points));
	FindGeometryType fg(fs);
	auto Moldat = fg.findMolDAT_Rad(*p_distances);
	return ConvertDATTuple(Moldat, fs);
}

// Single thread function
static void ChildThreadFunc(const SearchGraphType::RequestGraphType& input, const SearchGraphType::AtomIndex MaxAtom, SearchDataInterfaceType& dataInterface, const bool exact) {
	SearchGraphType graph;
	while (true) {
		auto next = dataInterface.getNext();
		if (next == nullptr) {
			return;
		}
		const auto& multi = dataInterface.getMulty();
		auto map = input.getTypeMap();
		graph.setupInput(input.makeCopy()); 
		auto && molData = SearchGraphType::DatabaseGraphType::ReadData(next, multi, map);
		if (!molData.second) continue;
		auto id = molData.first.getID();
		graph.setupData(std::move(molData.first));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}

static void reorder(cpplib::FindGeometry::tupleDistance& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);

	if (std::get<0>(d) > std::get<1>(d))
		std::swap(std::get<0>(d), std::get<1>(d));
}
static void reorder(cpplib::FindGeometry::tupleAngle& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);
	std::get<2>(d) = std::get<0>(fs.parseIndex[std::get<2>(d)]);

	if (std::get<0>(d) > std::get<2>(d))
		std::swap(std::get<0>(d), std::get<2>(d));
}
static void reorder(cpplib::FindGeometry::tupleTorsion& d, const cpplib::currents::FAMStructType& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);
	std::get<2>(d) = std::get<0>(fs.parseIndex[std::get<2>(d)]);
	std::get<3>(d) = std::get<0>(fs.parseIndex[std::get<3>(d)]);

	if ((std::get<0>(d) > std::get<3>(d)) || (std::get<0>(d) == std::get<3>(d) && std::get<1>(d) > std::get<2>(d))) {
		std::swap(std::get<0>(d), std::get<3>(d));
		std::swap(std::get<1>(d), std::get<2>(d));
	}
}
template<class T,class I>
static void eraseDoubles(std::vector<T>& vec, 
						 typename std::function<bool(I, I)> comp) {
	if (vec.empty()) return;
	std::sort(vec.begin(), vec.end());
	I it = vec.begin();
	I it2 = (++(vec.begin()));
	while (it2 != vec.end())
	{
		if (comp(it, it2)) {
			vec.erase(it2);
			it2 = it;
			it2++;
		}
		else {
			it++;
			it2++;
		}
	}
}
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple& dat, const cpplib::currents::FAMStructType& fs) {
	auto& dists = std::get<0>(dat);
	auto s_dists = dists.size();
	auto& angles = std::get<1>(dat);
	auto s_angles = angles.size();
	auto& tors = std::get<2>(dat);
	auto s_tors = tors.size();

	for (size_t i = 0; i < s_dists; i++)
	{
		deb_write("ConvertDATTuple: reorder dist: ", i);
		reorder(dists[i], fs);
	}
	for (size_t i = 0; i < s_angles; i++)
	{
		deb_write("ConvertDATTuple: reorder angl: ", i);
		reorder(angles[i], fs);
	}
	for (size_t i = 0; i < s_tors; i++)
	{
		deb_write("ConvertDATTuple: reorder tors: ", i);
		reorder(tors[i], fs);
	}
	// erase dublicates

	deb_write("ConvertDATTuple: erase dist");
	eraseDoubles(dists, std::function<bool(std::remove_reference<decltype(dists)>::type::iterator, std::remove_reference< decltype(dists)>::type::iterator)>(
		[](std::remove_reference< decltype(dists)>::type::iterator it, std::remove_reference< decltype(dists)>::type::iterator it2)
		{return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (::std::abs(std::get<2>(*it) - std::get<2>(*it2)) < 0.0001); }));
	deb_write("ConvertDATTuple: erase angl");
	eraseDoubles(angles, std::function<bool(std::remove_reference< decltype(angles)>::type::iterator, std::remove_reference< decltype(angles)>::type::iterator)>(
				 [](std::remove_reference< decltype(angles)>::type::iterator it, std::remove_reference< decltype(angles)>::type::iterator it2)
				 {return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (std::get<2>(*it) == std::get<2>(*it2)) && (::std::abs(std::get<3>(*it) - std::get<3>(*it2)) < 0.0001); }));
	deb_write("ConvertDATTuple: erase tors"); 
	eraseDoubles(tors, std::function<bool(std::remove_reference< decltype(tors)>::type::iterator, std::remove_reference< decltype(tors)>::type::iterator)>(
				 [](std::remove_reference< decltype(tors)>::type::iterator it, std::remove_reference< decltype(tors)>::type::iterator it2)
				 {return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) &&
				 (std::get<2>(*it) == std::get<2>(*it2)) && (std::get<3>(*it) == std::get<3>(*it2)) && (::std::abs(std::get<4>(*it) - std::get<4>(*it2)) < 0.0001);}));
	return dat;
}

std::tuple<std::vector<cpplib::currents::PointType>, std::list<std::string>> Compaq(const std::array<cpplib::currents::FloatingPointType, 6>& unit_cell,
												const std::vector<const char*>& symm,
												cpplib::currents::FAMStructType::AtomContainerType& types,
												cpplib::currents::FAMStructType::PointConteinerType& points) {
	deb_write("Compaq invoked");
	auto& distances = *p_distances;
	if (p_distances->isReady() == false) {
		return std::make_tuple(std::vector<cpplib::currents::PointType>(), std::list<std::string>(1, "Error!Could not open BondLength.ini"));
	}
	deb_write("Compaq p_distances prepared");
	FAMStructType fs;
	deb_write("Compaq FAMCellType creation start");
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	deb_write("Compaq call constructor ParseDataType");
	ParseDataType(fs, fc, symm, std::move(types), std::move(points), false);

	const auto su = fs.sizeUnique;
	std::string errorMsg;
	std::list<std::string> res_errors;
	deb_write("Compaq call fs.findBonds");
	auto res = fs.findBonds(distances, errorMsg, [fc](const PointType& p1, const PointType& p2) {return fc.distanceInCell(p1, p2); });
	if (!errorMsg.empty()) res_errors.emplace_back(std::move(errorMsg));

	deb_write("Compaq create fm");
	FindMoleculesType fm(std::move(fs));
	deb_write("Compaq call fm.compaq");
	auto & compaqed = fm.compaq(distances, res.first);
	compaqed.resize(su);
	deb_write("Compaq return");
	return std::make_tuple(std::move(compaqed), res_errors);
}

std::vector<std::tuple<cpplib::currents::PointType, cpplib::currents::AtomIndex, long, cpplib::FAM_Struct::ShiftType>> 
	ClusterCreate(std::array<cpplib::currents::FloatingPointType, 6> unit_cell,
				  const std::vector<const char*>& symm,
				  cpplib::currents::FAMStructType::AtomContainerType& types,
				  cpplib::currents::FAMStructType::PointConteinerType& points,
				  const std::vector<std::pair<cpplib::currents::PointType, cpplib::currents::FloatingPointType>>& anchors,
				  cpplib::currents::FloatingPointType over_radius,
				  bool& hasPolymer) {
	deb_write("ClusterCreate invoked");
	
	using ShiftType = FAMStructType::ShiftType;

	FAMStructType fs;
	FAMCellType fc(FAMCellType::base(unit_cell, true));
	
	ParseDataType(fs, fc, symm, std::move(types), std::move(points), false);

	deb_write("ClusterCreate: call findMoleculesForCluster");
	auto molecules = fc.findMoleculesForCluster(fs, *p_distances, hasPolymer); // [0] is empty
	deb_write("ClusterCreate: findMoleculesForCluster returns");
	for(const auto &m:molecules) {
		if(m.second) deb_write("ClusterCreate: poly true");
	}

	cpplib::Cluster cluster(fs, fc, anchors);
	cpplib::Cluster::BoxType box;
	cluster.CreateBox(box);
	deb_write("ClusterCreate: CreateBox returns");
	decltype(box) newbox;
	std::vector<bool> polyflags(molecules.size(),false);
	for (auto& pair : box) {
		for (FAMStructType::size_type i = 0; i < fs.sizePoints; i++) {
			if (pair.second[i]) continue;
			for (const auto& anch: anchors) {
				auto da = anch.first - fs.points[i];
				auto db = da - pair.first;
				if (anch.second >= (fc.fracToCart() * db).r()) {
					pair.second[i] = true;
					// find molecule
					MoleculeIndex m = 1;
					for (; m < molecules.size(); m++)
					{
						if (molecules[m].first[i].empty() == false) {
							break;
						}
					}
					_ASSERT(m != molecules.size());
					if(m == molecules.size()) {
						deb_write("ClusterCreate: m == molecules.size(): m = ", m);						
					}
					if (molecules[m].second == false) {
						cluster.Grow(box, newbox, molecules[m].first, i, pair.first);
					}
					else {
						polyflags[m] = true;
					}
				}
			}
		}
	}

	deb_write("ClusterCreate: concatination start");
	// concatinate boxes
	box.insert(newbox.begin(), newbox.end());
	
	cluster.GrowPoly(box, molecules, polyflags, over_radius);

	deb_write("ClusterCreate: concatination end");

	std::vector<std::tuple<PointType, cpplib::currents::AtomIndex, long, ShiftType>> ret; // Vector of (AtomIndex, SymmIndex, dx, dy, dz)
	//prepare ret
	for (auto& pair : box) {
		for (AtomIndex i = 0; i < fs.sizePoints; i++) {
			if (pair.second[i] == false) continue;
			auto point = pair.first + fs.points[i];
			auto shift = pair.first + std::get<2>(fs.parseIndex[i]);
			ret.emplace_back(pair.first + fs.points[i],
							 std::get<0>(fs.parseIndex[i]),
							 static_cast<long>(std::get<1>(fs.parseIndex[i])),
							 pair.first + std::get<2>(fs.parseIndex[i]));
		}
	}

	deb_write("ClusterCreate return");
	return ret;
}
