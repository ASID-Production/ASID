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
#include "Functions.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <functional>
#include <list>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../BaseHeaders/Currents.h"
#include "../BaseHeaders/DebugMes.h"

#include "../Classes/Cluster.h"
#include "../Classes/Distances.h"
#include "../Classes/FindGeometry.h"
#include "../Classes/FindMolecules.h"
#include "../Classes/Geometry.h"
#include "../Classes/Interfaces.h"
#include "../Classes/MoleculeGraph.h"
#include "../Classes/SearchGraph.h"



using namespace cpplib;
using namespace cpplib::currents;
using namespace cpplib::basic_types;
using PointType = FAM_Cell::PointType;

static void ChildThreadFunc(const SearchGraph::RequestGraphType& input, const SearchGraph::AtomIndex MaxAtom, SearchDataInterface& dataInterface, const bool exact);
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple& dat, const cpplib::FAM_Struct& fs);


const Distances* p_distances = nullptr;

/**
 * @brief Compares two molecule inputs by constructing search graphs and running a full graph search.
 *
 * Builds a query graph from `search1` and a data graph from `search2`, prepares them for searching,
 * and executes a full graph search using the specified matching strictness.
 *
 * @param search1 C-style string identifying the input used to build the query graph.
 * @param search2 C-style string identifying the input used to build the data graph.
 * @param exact When `true`, requires an exact match; when `false`, allows inexact matching.
 * @return true if the full search reports a successful match, `false` otherwise (also `false` if parsing `search2` fails).
 */
bool CompareGraph(const char* search1, const char* search2, const bool exact) {
	SearchGraph graph;

	auto&& inputpair = WITH_LOG(cpplib::MoleculeParser<AtomTypeRequest>::Read<AtomTypeRequest>, search1);
	auto map = inputpair.first.getTypeMap();
	graph.setupInput(std::move(inputpair.first));
	auto&& d_pair = WITH_LOG(cpplib::MoleculeParser<AtomTypeData>::Read<AtomTypeData>, search2, inputpair.second, map);

	if (!d_pair.second) return false;
	graph.setupData(std::move(d_pair.first));
	graph.prepareToSearch();
	return graph.startFullSearch(exact);
}
/**
 * @brief Execute a graph-based search over multiple input molecules, possibly in parallel.
 *
 * Parses the provided search specification, processes the input molecule data with a SearchDataInterface,
 * dispatches worker threads that run ChildThreadFunc to perform the search, and collects matching result ids.
 *
 * @param search Null-terminated string containing the search specification to parse.
 * @param data Vector of C-style string pointers representing input molecule entries; ownership of the container is moved.
 * @param np Requested maximum number of threads to use (clamped to hardware concurrency and input size).
 * @param exact If `true`, require exact matches during the search; otherwise allow inexact matches.
 * @return std::vector<int> Vector of result identifiers produced by the search (order as provided by SearchDataInterface).
 */
std::vector<int> SearchMain(const char* search, std::vector<const char*>&& data, const int np, const bool exact) {

	auto&& inputpair = WITH_LOG(cpplib::MoleculeParser<AtomTypeRequest>::Read<AtomTypeRequest>, search);
	SearchDataInterface databuf(std::move(data), std::move(inputpair.second));
	std::vector<std::thread> threads;
	const size_t nThreads = std::min(std::min(static_cast<unsigned int>(np), std::thread::hardware_concurrency()),
									 static_cast<unsigned int>(databuf.size())) - 1;
	threads.reserve(nThreads);
	auto ma = 1;

	for (size_t i = 0; i < nThreads; i++) {
		threads.emplace_back(ChildThreadFunc, std::cref(inputpair.first), ma, std::ref(databuf), exact);
	}
	ChildThreadFunc(inputpair.first, ma, databuf, exact);

	for (size_t i = 0; i < nThreads; i++) {
		threads[i].join();
	}

	return databuf.getAllResults();
}

std::tuple<std::string, std::string, FindMolecules::RightType> FindMoleculesInCell(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
																		 std::vector<const char*>& symm,
																		 cpplib::FAM_Struct::AtomContainerType& types,
																		 cpplib::FAM_Struct::PointConteinerType& points) {
	auto& distances = *p_distances;
	if (p_distances->isReady() == false) {
		return std::make_tuple(std::string(), std::string("Error! Could not open BondLength.ini"),
							  FindMolecules::RightType());
	}
	FAM_Struct fs;
	FAM_Cell fc(geometry::Cell(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points));

	fc.CreateSupercell(fs.points, fs.findCutoff(distances));

	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, [fc](const PointType& p1, const PointType& p2) {return fc.distanceInCell(p1, p2); });

	FindMolecules fm(std::move(fs));

	auto ret = fm.findMolecules(res.first, res.second, errorMsg);

	// For Petr's enjoyment
	auto mol_s = std::get<2>(ret).size();
	for (size_t i = 0; i < mol_s; i++) {
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

	for (size_t i = 0; i < mol_s; i++) {
		auto j_s = std::get<0>(std::get<2>(ret)[i]).size();
		for (size_t j = 0; j < j_s; j++) {
			auto& point = std::get<0>(std::get<0>(std::get<2>(ret)[i])[j]);
			point = fc.fracToCart() * point;
		}
	}
	return ret;
}
std::tuple<std::string, std::string, FindMolecules::RightType>  FindMoleculesWithoutCell(cpplib::FAM_Struct::AtomContainerType& types,
																			  cpplib::FAM_Struct::PointConteinerType& points) {
	auto& distances = *p_distances;

	if (p_distances->isReady() == false) {
		return std::make_tuple(std::string(), std::string("Error! Could not open BondLength.ini"),
							   FindMolecules::RightType());
	}

	FAM_Struct fs;
	ParseData(fs, std::move(types), std::move(points));
	std::string errorMsg;
	auto res = fs.findBonds(distances, errorMsg, PointType::distance);

	FindMolecules fm(std::move(fs));

	return fm.findMolecules(res.first, res.second, errorMsg);
}

std::vector<FindGeometry::tupleDistance> FindDistanceWC(cpplib::FAM_Struct::AtomContainerType& types,
						   cpplib::FAM_Struct::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value) {
	FAM_Struct fs;
	ParseData(fs, std::move(types), std::move(points));
	FindGeometry fg(fs);
	const auto raw = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[0]),
										static_cast<FindGeometry::AtomTypeBase>(type[1]),
										value);
	return raw;
}

std::vector<FindGeometry::tupleDistance> FindDistanceIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
						   std::vector<const char*>& symm,
						   cpplib::FAM_Struct::AtomContainerType& types,
						   cpplib::FAM_Struct::PointConteinerType& points,
						   const std::array<int, 2>& type,
						   const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value) {
	FAM_Struct fs;
	FAM_Cell fc(FAM_Cell::base(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points));

	fc.CreateSupercell(fs.points, static_cast<basic_types::FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++) {
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometry fg(fs);
	const auto raw = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(std::get<0>(type)),
										static_cast<FindGeometry::AtomTypeBase>(std::get<1>(type)),
										value);
	return raw;
}

std::vector<FindGeometry::tupleAngle> FindAngleWC(cpplib::FAM_Struct::AtomContainerType& types,
													  cpplib::FAM_Struct::PointConteinerType& points,
													  const std::array<int, 3>& type,
													  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_d,
													  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_a) {
	FAM_Struct fs;
	ParseData(fs, std::move(types), std::move(points));
	FindGeometry fg(fs);
	const auto raw12 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(std::get<0>(type)),
										static_cast<FindGeometry::AtomTypeBase>(std::get<1>(type)),
										value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(std::get<1>(type)),
										static_cast<FindGeometry::AtomTypeBase>(std::get<2>(type)),
										value_d[1]);
	const auto raw = fg.findAngle(raw12, raw23, std::make_pair(cpplib::geometry::GradtoRad(value_a.first), cpplib::geometry::GradtoRad(value_a.second)));
	return raw;
}

std::vector<FindGeometry::tupleAngle> FindAngleIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
													  std::vector<const char*>& symm,
													  cpplib::FAM_Struct::AtomContainerType& types,
													  cpplib::FAM_Struct::PointConteinerType& points,
													  const std::array<int, 3>& type,
													  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_d,
													  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_a) {
	FAM_Struct fs;
	FAM_Cell fc(FAM_Cell::base(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points));

	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++) {
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometry fg(fs);
	const auto raw12 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(std::get<0>(type)),
										static_cast<FindGeometry::AtomTypeBase>(std::get<1>(type)),
										value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(std::get<1>(type)),
										static_cast<FindGeometry::AtomTypeBase>(std::get<2>(type)),
										value_d[1]);
	auto raw = fg.findAngle(raw12, raw23, std::make_pair(cpplib::geometry::GradtoRad(value_a.first), cpplib::geometry::GradtoRad(value_a.second)));
	return raw;
}

std::vector<FindGeometry::tupleTorsion> FindTorsionWC(cpplib::FAM_Struct::AtomContainerType& types,
														  cpplib::FAM_Struct::PointConteinerType& points,
														  const std::array<int, 4>& type,
														  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 3>& value_d,
														  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_a,
														  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_t) {
	FAM_Struct fs;
	ParseData(fs, std::move(types), std::move(points));
	FindGeometry fg(fs);

	const auto raw12 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[0]), static_cast<FindGeometry::AtomTypeBase>(type[1]), value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[1]), static_cast<FindGeometry::AtomTypeBase>(type[2]), value_d[1]);
	const auto raw34 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[2]), static_cast<FindGeometry::AtomTypeBase>(type[3]), value_d[2]);
	const auto raw123 = fg.findAngle(raw12, raw23,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[0].first), cpplib::geometry::GradtoRad(value_a[0].second)));
	const auto raw234 = fg.findAngle(raw23, raw34,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[1].first), cpplib::geometry::GradtoRad(value_a[1].second)));
	auto raw = fg.findTorsion(raw123, raw234,
							  std::make_pair(cpplib::geometry::GradtoRad(value_t.first), cpplib::geometry::GradtoRad(value_t.second)));
	return raw;
}
std::vector<FindGeometry::tupleTorsion> FindTorsionIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
														  std::vector<const char*>& symm,
														  cpplib::FAM_Struct::AtomContainerType& types,
														  cpplib::FAM_Struct::PointConteinerType& points,
														  const std::array<int, 4>& type,
														  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 3>& value_d,
														  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_a,
														  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_t) {
	FAM_Struct fs;
	FAM_Cell fc(FAM_Cell::base(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points));

	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++) {
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometry fg(fs);

	const auto raw12 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[0]), static_cast<FindGeometry::AtomTypeBase>(type[1]), value_d[0]);
	const auto raw23 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[1]), static_cast<FindGeometry::AtomTypeBase>(type[2]), value_d[1]);
	const auto raw34 = fg.findDistance(static_cast<FindGeometry::AtomTypeBase>(type[2]), static_cast<FindGeometry::AtomTypeBase>(type[3]), value_d[2]);
	const auto raw123 = fg.findAngle(raw12, raw23,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[0].first), cpplib::geometry::GradtoRad(value_a[0].second)));
	const auto raw234 = fg.findAngle(raw23, raw34,
									 std::make_pair(cpplib::geometry::GradtoRad(value_a[1].first), cpplib::geometry::GradtoRad(value_a[1].second)));
	const auto raw = fg.findTorsion(raw123, raw234,
									std::make_pair(cpplib::geometry::GradtoRad(value_t.first), cpplib::geometry::GradtoRad(value_t.second)));
	return raw;
}

cpplib::DATTuple FindDAT_IC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
							std::vector<const char*>& symm,
							cpplib::FAM_Struct::AtomContainerType& types,
							cpplib::FAM_Struct::PointConteinerType& points) {
	FAM_Struct fs;
	FAM_Cell fc(FAM_Cell::base(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points));


	fc.CreateSupercell(fs.points, static_cast<FloatingPointType>(8.5), 2);

	for (size_t i = 0; i < fs.sizePoints; i++) {
		fs.points[i] = fc.fracToCart() * fs.points[i];
	}

	FindGeometry fg(fs);
	auto Moldat = fg.findMolDAT_Rad(*p_distances);
	return ConvertDATTuple(Moldat, fs);
}

cpplib::DATTuple FindDAT_WC(cpplib::FAM_Struct::AtomContainerType& types,
							cpplib::FAM_Struct::PointConteinerType& points) {
	FAM_Struct fs;
	ParseData(fs, std::move(types), std::move(points));
	FindGeometry fg(fs);
	auto Moldat = fg.findMolDAT_Rad(*p_distances);
	return ConvertDATTuple(Moldat, fs);
}

/**
 * @brief Worker loop that processes molecule requests and runs graph searches.
 *
 * Continuously obtains the next input item from dataInterface, parses it into a molecule,
 * configures a SearchGraph from the provided input template and parsed molecule data,
 * runs a full search constrained by MaxAtom and the exact flag, and pushes the molecule ID
 * to dataInterface when the search succeeds. The loop exits when dataInterface has no more items.
 *
 * @param input Request graph template used to configure each SearchGraph instance.
 * @param MaxAtom Maximum atom index/count used to constrain the search.
 * @param dataInterface Source of work items and sink for successful result IDs.
 * @param exact When true, perform an exact search match; otherwise allow inexact matching.
 */
static void ChildThreadFunc(const SearchGraph::RequestGraphType& input, const SearchGraph::AtomIndex MaxAtom, SearchDataInterface& dataInterface, const bool exact) {
	SearchGraph graph;
	while (true) {
		auto next = dataInterface.getNext();
		if (next == nullptr) {
			return;
		}
		const auto& multi = dataInterface.getMulty();
		auto map = input.getTypeMap();
		auto tempinput = input;
		graph.setupInput(std::move(tempinput));
		auto&& molData = WITH_LOG(cpplib::MoleculeParser<AtomTypeData>::Read<AtomTypeData>, next, multi, map);
		if (!molData.second) continue;
		auto id = molData.first.getID();
		graph.setupData(std::move(molData.first));
		graph.prepareToSearch();
		if (graph.startFullSearch(exact, MaxAtom)) {
			dataInterface.push_result(id);
		}
	}
}

static void reorder(cpplib::FindGeometry::tupleDistance& d, const cpplib::FAM_Struct& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);

	if (std::get<0>(d) > std::get<1>(d))
		std::swap(std::get<0>(d), std::get<1>(d));
}
static void reorder(cpplib::FindGeometry::tupleAngle& d, const cpplib::FAM_Struct& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);
	std::get<2>(d) = std::get<0>(fs.parseIndex[std::get<2>(d)]);

	if (std::get<0>(d) > std::get<2>(d))
		std::swap(std::get<0>(d), std::get<2>(d));
}
static void reorder(cpplib::FindGeometry::tupleTorsion& d, const cpplib::FAM_Struct& fs) {
	std::get<0>(d) = std::get<0>(fs.parseIndex[std::get<0>(d)]);
	std::get<1>(d) = std::get<0>(fs.parseIndex[std::get<1>(d)]);
	std::get<2>(d) = std::get<0>(fs.parseIndex[std::get<2>(d)]);
	std::get<3>(d) = std::get<0>(fs.parseIndex[std::get<3>(d)]);

	if ((std::get<0>(d) > std::get<3>(d)) || (std::get<0>(d) == std::get<3>(d) && std::get<1>(d) > std::get<2>(d))) {
		std::swap(std::get<0>(d), std::get<3>(d));
		std::swap(std::get<1>(d), std::get<2>(d));
	}
}
template<class T, class I>
static void eraseDoubles(std::vector<T>& vec,
						 typename std::function<bool(I, I)> comp) {
	if (vec.empty()) return;
	std::sort(vec.begin(), vec.end());
	I it = vec.begin();
	I it2 = (++(vec.begin()));
	while (it2 != vec.end()) {
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
/**
 * @brief Reorders index tuples in a DATTuple to parsed indices and removes duplicate entries.
 *
 * Rewrites distance, angle, and torsion tuples in-place so their atom indices match the
 * parsed indexing defined by fs, and then eliminates duplicate tuples. For numeric
 * components (distance/angle/torsion values) duplicates are determined with a small
 * tolerance (1e-4).
 *
 * @param dat The DATTuple to normalize; its distance, angle, and torsion containers are modified.
 * @param fs  FAM_Struct providing the mapping from internal to parsed indices used for reordering.
 * @return cpplib::DATTuple& Reference to the modified input `dat`.
 */
static cpplib::DATTuple& ConvertDATTuple(cpplib::DATTuple& dat, const cpplib::FAM_Struct& fs) {
	auto& dists = std::get<0>(dat);
	auto s_dists = dists.size();
	auto& angles = std::get<1>(dat);
	auto s_angles = angles.size();
	auto& tors = std::get<2>(dat);
	auto s_tors = tors.size();

	for (size_t i = 0; i < s_dists; i++) {
		reorder(dists[i], fs);
	}
	for (size_t i = 0; i < s_angles; i++) {
		reorder(angles[i], fs);
	}
	for (size_t i = 0; i < s_tors; i++) {
		reorder(tors[i], fs);
	}
	// erase dublicates

	eraseDoubles(dists, std::function<bool(std::remove_reference<decltype(dists)>::type::iterator, std::remove_reference< decltype(dists)>::type::iterator)>(
		[](std::remove_reference< decltype(dists)>::type::iterator it, std::remove_reference< decltype(dists)>::type::iterator it2)
		{return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (::std::abs(std::get<2>(*it) - std::get<2>(*it2)) < 0.0001); }));

	eraseDoubles(angles, std::function<bool(std::remove_reference< decltype(angles)>::type::iterator, std::remove_reference< decltype(angles)>::type::iterator)>(
		[](std::remove_reference< decltype(angles)>::type::iterator it, std::remove_reference< decltype(angles)>::type::iterator it2)
		{return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) && (std::get<2>(*it) == std::get<2>(*it2)) && (::std::abs(std::get<3>(*it) - std::get<3>(*it2)) < 0.0001); }));

	eraseDoubles(tors, std::function<bool(std::remove_reference< decltype(tors)>::type::iterator, std::remove_reference< decltype(tors)>::type::iterator)>(
		[](std::remove_reference< decltype(tors)>::type::iterator it, std::remove_reference< decltype(tors)>::type::iterator it2)
		{return (std::get<0>(*it) == std::get<0>(*it2)) && (std::get<1>(*it) == std::get<1>(*it2)) &&
		(std::get<2>(*it) == std::get<2>(*it2)) && (std::get<3>(*it) == std::get<3>(*it2)) && (::std::abs(std::get<4>(*it) - std::get<4>(*it2)) < 0.0001); }));
	return dat;
}

/**
 * @brief Compacts atom coordinates into molecule fragments and reports processing errors.
 *
 * Builds a FAM structure from the provided unit cell, symmetry operations, atom types, and fractional points;
 * finds bonds using global distance data; groups atoms into compacted molecular fragments; and returns their
 * Cartesian coordinates along with any error messages produced during processing.
 *
 * @param unit_cell Six-element array describing the unit cell (a, b, c, alpha, beta, gamma).
 * @param symm Vector of symmetry operator strings defining the space-group symmetry.
 * @param types Container of atom type descriptors; ownership is taken by this function (moved).
 * @param points Container of fractional atomic positions; ownership is taken by this function (moved).
 *
 * @return std::tuple<std::vector<cpplib::geometry::Point<FloatingPointType>>, std::list<std::string>>
 *   - First element: vector of compacted atom coordinates in Cartesian space (one entry per atom in the compacted output).
 *   - Second element: list of error messages produced during processing. If the global BondLength.ini is unavailable,
 *     the first element is an empty vector and the list contains a single message "Error!Could not open BondLength.ini".
 */
std::tuple<std::vector<cpplib::geometry::Point<FloatingPointType>>, std::list<std::string>> Compaq(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
												const std::vector<const char*>& symm,
												cpplib::FAM_Struct::AtomContainerType& types,
												cpplib::FAM_Struct::PointConteinerType& points) {
	auto& distances = *p_distances;
	if (p_distances->isReady() == false) {
		return std::make_tuple(std::vector<cpplib::geometry::Point<FloatingPointType>>(), std::list<std::string>(1, "Error!Could not open BondLength.ini"));
	}
	FAM_Struct fs;
	FAM_Cell fc(FAM_Cell::base(unit_cell, true));
	ParseData(fs, fc, symm, std::move(types), std::move(points), false);

	const auto su = fs.sizeUnique;
	std::string errorMsg;
	std::list<std::string> res_errors;
	auto res = fs.findBonds(distances, errorMsg, [fc](const PointType& p1, const PointType& p2) {return fc.distanceInCell(p1, p2); });
	if (!errorMsg.empty()) res_errors.emplace_back(std::move(errorMsg));

	FindMolecules fm(std::move(fs));
	auto& compaqed = fm.compaq(res.first);
	compaqed.resize(su);
	return std::make_tuple(std::move(compaqed), res_errors);
}

/**
 * @brief Build and execute a cluster search for the given cell, symmetry, atoms, and points.
 *
 * Validates that global distance data is available, constructs a Cluster from the provided
 * unit cell, symmetry definitions, anchors, points, and atom types, and returns the cluster
 * execution result while updating polymer presence information.
 *
 * @param unit_cell Six-component unit cell array (a, b, c, alpha, beta, gamma).
 * @param symm Vector of symmetry specification strings to apply to the cell.
 * @param types Atom-type container describing species for each site.
 * @param points Point container with atomic positions (fractional coordinates for the cell).
 * @param anchors Cluster anchor definitions used to seed the cluster construction.
 * @param polymer_cutoff Distance threshold used to detect polymeric connectivity.
 * @param[out] hasPolymer Set to `true` if the executed cluster contains polymeric chains, `false` otherwise.
 * @return Cluster::ClusterData Resulting cluster data produced by the cluster execution; returns an empty ClusterData if global distance data is not ready.
 */
Cluster::ClusterData ClusterCreate(std::array<cpplib::basic_types::FloatingPointType, 6> unit_cell,
			                       const std::vector<const char*>& symm,
			                       cpplib::FAM_Struct::AtomContainerType& types,
			                       cpplib::FAM_Struct::PointConteinerType& points,
			                       std::vector<cpplib::Cluster::AnchorType>& anchors,
			                       cpplib::basic_types::FloatingPointType polymer_cutoff,
			                       bool& hasPolymer) {

	using ShiftType = Cluster::ShiftType;
	auto& distances = *p_distances;
	if (p_distances->isReady() == false) {
		{
			return {};
		}
	}

	cpplib::geometry::Cell cell(unit_cell);
	std::vector<geometry::Symm<FloatingPointType>> symms;
	symms.reserve(symm.size());
	for (const auto& s : symm) {
		symms.emplace_back(s);
	}

	Cluster cluster(cell, symms, std::move(anchors), std::move(points), std::move(types), polymer_cutoff);

	auto ret = cluster.execute(distances, hasPolymer);

	return ret;
}
