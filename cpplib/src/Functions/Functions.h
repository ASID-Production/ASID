#pragma once
#include <vector>
#include <string>
#include "AllInOneAndCurrent.h"

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
	FindMoleculesInCell(const std::array<float, 6>& unit_cell,
						std::vector<const char*>& symm,
						cpplib::currents::FAMStructType::AtomContainerType& types,
						cpplib::currents::FAMStructType::PointConteinerType& points);
std::tuple<std::string, std::string, cpplib::FindMolecules::RightType>
	FindMoleculesWithoutCell(cpplib::currents::FAMStructType::AtomContainerType& types,
							 cpplib::currents::FAMStructType::PointConteinerType& points);


std::vector<cpplib::currents::FindGeometryType::tupleDistance> 
	FindDistanceWC(cpplib::currents::FAMStructType::AtomContainerType& types,
				   cpplib::currents::FAMStructType::PointConteinerType& points,
				   const std::array<int, 2>& type,
				   const std::pair<float, float>& value);

std::vector<cpplib::currents::FindGeometryType::tupleDistance>
	FindDistanceIC(const std::array<float, 6>& unit_cell,
				   std::vector<const char*>& symm,
				   cpplib::currents::FAMStructType::AtomContainerType& types,
				   cpplib::currents::FAMStructType::PointConteinerType& points,
				   const std::array<int, 2>& type,
				   const std::pair<float, float>& value);

std::vector<cpplib::currents::FindGeometryType::tupleAngle> 
	FindAngleWC(cpplib::currents::FAMStructType::AtomContainerType& types,
				cpplib::currents::FAMStructType::PointConteinerType& points,
				const std::array<int, 3>& type,
				const std::array<std::pair<float, float>, 2>& value_d,
				const std::pair<float, float>& value_a);

std::vector<cpplib::currents::FindGeometryType::tupleAngle> 
	FindAngleIC(const std::array<float, 6>& unit_cell,
				std::vector<const char*>& symm, 
				cpplib::currents::FAMStructType::AtomContainerType& types,
				cpplib::currents::FAMStructType::PointConteinerType& points,
				const std::array<int, 3>& type,
				const std::array<std::pair<float, float>, 2>& value_d,
				const std::pair<float, float>& value_a);

std::vector<cpplib::currents::FindGeometryType::tupleTorsion> 
	FindTorsionWC(cpplib::currents::FAMStructType::AtomContainerType& types,
				  cpplib::currents::FAMStructType::PointConteinerType& points,
	              const std::array<int, 4>& type,
				  const std::array<std::pair<float, float>, 3>& value_d,
				  const std::array<std::pair<float, float>, 2>& value_a,
				  const std::pair<float, float>& value_t);

std::vector<cpplib::currents::FindGeometryType::tupleTorsion> 
	FindTorsionIC(const std::array<float, 6>& unit_cell,
				  std::vector<const char*>& symm, 
				  cpplib::currents::FAMStructType::AtomContainerType& types,
				  cpplib::currents::FAMStructType::PointConteinerType& points,
				  const std::array<int, 4>& type,
				  const std::array<std::pair<float, float>, 3>& value_d,
				  const std::array<std::pair<float, float>, 2>& value_a,
				  const std::pair<float, float>& value_t);


cpplib::DATTuple FindDAT_IC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							cpplib::currents::FAMStructType::AtomContainerType& types,
							cpplib::currents::FAMStructType::PointConteinerType& points);
cpplib::DATTuple FindDAT_WC(cpplib::currents::FAMStructType::AtomContainerType& types,
							cpplib::currents::FAMStructType::PointConteinerType& points);

std::tuple<std::vector<cpplib::currents::PointType>, std::list<std::string>> Compaq(const std::array<float, 6>& unit_cell,
												const std::vector<const char*>& symm,
												cpplib::currents::FAMStructType::AtomContainerType& types,
												cpplib::currents::FAMStructType::PointConteinerType& points);


std::vector<std::tuple<cpplib::currents::AtomIndex,long,long,long,long>> ClusterCreate();