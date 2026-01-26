#pragma once
#include <array>
#include <list>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"

#include "../Classes/Cluster.h"
#include "../Classes/FindGeometry.h"
#include "../Classes/FindMolecules.h"
#include "../Classes/Geometry.h"

namespace cpplib {
	using DATTuple = ::std::tuple<::std::vector<FindGeometry::tupleDistance>, ::std::vector<FindGeometry::tupleAngle>, ::std::vector<FindGeometry::tupleTorsion>>;
}

bool CompareGraph(const char* search1,
					const char* search2,
					const bool exact);

std::vector<int> SearchMain(const char* search,
							std::vector <const char*>&& data,
							const int np,
							const bool exact);

std::tuple<std::string, std::string, cpplib::FindMolecules::RightType>
FindMoleculesInCell(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
					std::vector<const char*>& symm,
					cpplib::FAM_Struct::AtomContainerType& types,
					cpplib::FAM_Struct::PointConteinerType& points);
std::tuple<std::string, std::string, cpplib::FindMolecules::RightType>
FindMoleculesWithoutCell(cpplib::FAM_Struct::AtomContainerType& types,
						 cpplib::FAM_Struct::PointConteinerType& points);


std::vector<cpplib::FindGeometry::tupleDistance>
FindDistanceWC(cpplib::FAM_Struct::AtomContainerType& types,
			   cpplib::FAM_Struct::PointConteinerType& points,
			   const std::array<int, 2>& type,
			   const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value);

std::vector<cpplib::FindGeometry::tupleDistance>
FindDistanceIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
			   std::vector<const char*>& symm,
			   cpplib::FAM_Struct::AtomContainerType& types,
			   cpplib::FAM_Struct::PointConteinerType& points,
			   const std::array<int, 2>& type,
			   const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value);

std::vector<cpplib::FindGeometry::tupleAngle>
FindAngleWC(cpplib::FAM_Struct::AtomContainerType& types,
			cpplib::FAM_Struct::PointConteinerType& points,
			const std::array<int, 3>& type,
			const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_d,
			const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_a);

std::vector<cpplib::FindGeometry::tupleAngle>
FindAngleIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
			std::vector<const char*>& symm,
			cpplib::FAM_Struct::AtomContainerType& types,
			cpplib::FAM_Struct::PointConteinerType& points,
			const std::array<int, 3>& type,
			const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_d,
			const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_a);

std::vector<cpplib::FindGeometry::tupleTorsion>
FindTorsionWC(cpplib::FAM_Struct::AtomContainerType& types,
			  cpplib::FAM_Struct::PointConteinerType& points,
			  const std::array<int, 4>& type,
			  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 3>& value_d,
			  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_a,
			  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_t);

std::vector<cpplib::FindGeometry::tupleTorsion>
FindTorsionIC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
			  std::vector<const char*>& symm,
			  cpplib::FAM_Struct::AtomContainerType& types,
			  cpplib::FAM_Struct::PointConteinerType& points,
			  const std::array<int, 4>& type,
			  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 3>& value_d,
			  const std::array<std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>, 2>& value_a,
			  const std::pair<cpplib::basic_types::FloatingPointType, cpplib::basic_types::FloatingPointType>& value_t);

cpplib::DATTuple FindDAT_IC(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
							std::vector<const char*>& symm,
							cpplib::FAM_Struct::AtomContainerType& types,
							cpplib::FAM_Struct::PointConteinerType& points);
cpplib::DATTuple FindDAT_WC(cpplib::FAM_Struct::AtomContainerType& types,
							cpplib::FAM_Struct::PointConteinerType& points);

std::tuple<std::vector<cpplib::geometry::Point<cpplib::basic_types::FloatingPointType>>, std::list<std::string>> Compaq(const std::array<cpplib::basic_types::FloatingPointType, 6>& unit_cell,
												const std::vector<const char*>& symm,
												cpplib::FAM_Struct::AtomContainerType& types,
												cpplib::FAM_Struct::PointConteinerType& points);

cpplib::cluster_detail::ClusterData ClusterCreate(std::array<cpplib::basic_types::FloatingPointType, 6> unit_cell,
														const std::vector<const char*>& symm,
														cpplib::FAM_Struct::AtomContainerType& types,
														cpplib::FAM_Struct::PointConteinerType& points,
														std::vector<cpplib::Cluster::AnchorType>& anchors,
														cpplib::basic_types::FloatingPointType polymer_cutoff,
														bool& hasPolymer);