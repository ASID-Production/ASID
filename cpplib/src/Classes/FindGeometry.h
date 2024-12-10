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
#pragma once
#include "Geometry.h"
#include "FindMolecules.h"
#include "Distances.h"
#include <vector>
#include <tuple>

namespace cpplib {
	class FindGeometry {
	public:
		// Definitions
		using FAMStructType = currents::FAMStructType;
		
		using AtomType = FAMStructType::AtomType;
		using FloatingPointType = FAMStructType::FloatingPointType;
		using AtomIndex = FAMStructType::AtomIndex;
		using size_type = FAMStructType::size_type;
		using PointType = FAMStructType::PointType;
		using DistancesType = FAMStructType::DistancesType;

		using tupleDistance = ::std::tuple<AtomIndex, AtomIndex, FloatingPointType>;
		using tupleAngle = ::std::tuple<AtomIndex, AtomIndex, AtomIndex, FloatingPointType>;
		using tupleTorsion = ::std::tuple<AtomIndex, AtomIndex, AtomIndex, AtomIndex, FloatingPointType>;
		using MinMaxType = ::std::pair<FloatingPointType, FloatingPointType>; 

		static_assert(::std::is_same<typename DistancesType::AtomType, AtomType>::value, "AtomTypes of FAM_Struct and Distances should be the same");
	private:
		// Data
		const FAMStructType& fs;

	public:
		FindGeometry() = delete;
		constexpr explicit FindGeometry(const FAMStructType& famstr) noexcept : fs(famstr) {};
		auto findDistance(AtomType t1, AtomType t2, MinMaxType d12) const {
			std::vector<tupleDistance> res;
			const bool mirror = (t1 == t2);
			for (size_type i = 0; i < fs.sizePoints; i++)
			{
				if (t1 != fs.types[i]) continue;
				const size_type start = mirror ? i + 1 : 0;
				for (size_type j = start; j < fs.sizePoints; j++)
				{
					if (i == j) continue;
					if (t2 != fs.types[j]) continue;
					auto R = (fs.points[i] - fs.points[j]).r();
					if (inRange(d12, R))
						res.emplace_back(i, j, R);
				}
			}
			return res;
		}
		auto findAngle(const std::vector<tupleDistance>& v1, const std::vector<tupleDistance>& v2, MinMaxType A123) const {
			size_type v1s = v1.size();
			size_type v2s = v2.size();
			std::vector<tupleAngle> res;
			if (v1s == 0 || v2s == 0)
				return res;
			// flags
			const bool is_v1_mirror = (std::get<0>(v1[0]) == std::get<1>(v1[0]));
			const bool is_v2_mirror = (std::get<0>(v2[0]) == std::get<1>(v2[0]));

			// comparation loop
			for (size_type i = 0; i < v1s; i++) {
				for (size_type j = 0; j < v2s; j++) {
					if (std::get<1>(v1[i]) != std::get<0>(v2[j]))
						continue;
					if (std::get<0>(v1[i]) == std::get<1>(v2[j]))
						continue;
					auto rad = PointType::angleRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v1[i])], fs.points[std::get<1>(v2[j])]);
					if (!inRange(A123, rad))
						continue;
					res.emplace_back(std::get<0>(v1[i]), std::get<1>(v1[i]), std::get<1>(v2[j]), rad);
				}
			}
			if (is_v1_mirror) {
				for (size_type i = 0; i < v1s; i++) {
					for (size_type j = 0; j < v2s; j++) {
						if (std::get<0>(v1[i]) != std::get<0>(v2[j]))
							continue;
						if (std::get<1>(v1[i]) == std::get<1>(v2[j]))
							continue;
						auto rad = PointType::angleRad(fs.points[std::get<1>(v1[i])], fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v2[j])]);
						if (!inRange(A123, rad))
							continue;
						res.emplace_back(std::get<1>(v1[i]), std::get<0>(v1[i]), std::get<1>(v2[j]), rad);
					}
				}
			}
			if (is_v2_mirror) {
				for (size_type i = 0; i < v1s; i++) {
					for (size_type j = 0; j < v2s; j++) {
						if (std::get<1>(v1[i]) != std::get<1>(v2[j]))
							continue;
						if (std::get<0>(v1[i]) == std::get<0>(v2[j]))
							continue;
						auto rad = PointType::angleRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v2[j])], fs.points[std::get<0>(v2[j])]);
						if (!inRange(A123, rad))
							continue;
						res.emplace_back(std::get<0>(v1[i]), std::get<1>(v2[i]), std::get<0>(v2[j]), rad);
					}
				}
			}
			return res;
		}
		auto findTorsion(std::vector<tupleAngle> v1, std::vector<tupleAngle> v2, MinMaxType t1234) const {
			auto v1s = v1.size();
			auto v2s = v2.size();
			std::vector<tupleTorsion> res;
			if (v1s == 0 || v2s == 0)
				return res;
			for (size_type i = 0; i < v1s; i++)
			{
				for (size_type j = 0; j < v2s; j++)
				{
					if (std::get<1>(v1[i]) != std::get<0>(v2[j]) || std::get<2>(v1[i]) != std::get<1>(v2[j]) || std::get<0>(v1[i]) >= std::get<2>(v2[j]))
						continue;
					auto tor = PointType::torsionRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v1[i])], fs.points[std::get<1>(v2[j])], fs.points[std::get<2>(v2[j])]);
					if (inRange(t1234, tor))
						continue;
					res.emplace_back(std::get<0>(v1[i]), std::get<1>(v1[i]), std::get<1>(v2[j]), std::get<2>(v2[j]), tor);
				}
			}
			return res;
		}

		auto findMolD(const DistancesType& dist) const {
			std::vector<tupleDistance> res;
			for (size_type i = 0; i < fs.sizePoints; i++)
			{
				for (size_type j = i + 1; j < fs.sizePoints; j++)
				{
					auto R = (fs.points[i] - fs.points[j]).r();
					if (R < dist.maxDistance(fs.types[fs.parseIndex[i]], fs.types[fs.parseIndex[j]]))
						res.emplace_back(i, j, R);
				}
			}
			return res;
		}
		auto findMolDA_Rad(const DistancesType& dist) const {
			std::vector<tupleAngle> res;
			std::vector<tupleDistance> alldist = findMolD(dist);
			const size_type s = static_cast<size_type>(alldist.size());
			for (size_type i = 0; i < s; i++)
			{
				size_type j = i + 1;
				const auto indexA = ::std::get<0>(alldist[i]);
				const auto indexB = ::std::get<1>(alldist[i]);
				for (; (j < s) && indexA == ::std::get<0>(alldist[j]); j++) //a-b,a-c
				{
					const auto indexC = ::std::get<1>(alldist[j]);
					const auto rad = PointType::angleRad(fs.points[indexB], fs.points[indexA], fs.points[indexC]);
					res.emplace_back(indexB, indexA, indexC, rad);
				}
				for (; j < s; j++) //a-b,b-c and a-b,c-b
				{
					std::remove_const_t<decltype(indexA)> indexC;
					if (indexB == ::std::get<0>(alldist[j])) {
						indexC = ::std::get<1>(alldist[j]);
					}
					else if (indexB == ::std::get<1>(alldist[j])) {
						indexC = ::std::get<0>(alldist[j]);
					}
					else continue;

					const auto rad = PointType::angleRad(fs.points[indexA], fs.points[indexB], fs.points[indexC]);
					res.emplace_back(indexA, indexB, indexC, rad);
				}
			}
			return ::std::make_tuple(::std::move(alldist), ::std::move(res));
		}
		auto findMolDAT_Rad(const DistancesType& dist) const {
			std::vector<tupleTorsion> res;
			auto allAngAndDist = findMolDA_Rad(dist);
			auto& alldist = ::std::get<0>(allAngAndDist);
			auto& allang = ::std::get<1>(allAngAndDist);

			const size_type sd = static_cast<size_type>(alldist.size());
			const size_type sa = static_cast<size_type>(allang.size());
			for (size_type i = 0; i < sa; i++)
			{
				const auto indexA = ::std::get<0>(allang[i]);
				const auto indexB = ::std::get<1>(allang[i]);
				const auto indexC = ::std::get<2>(allang[i]);
				for (size_type j = 0; j < sd; j++) //a-b-c-d and d-a-b-c
				{
					const auto indexD0 = ::std::get<0>(alldist[j]);
					const auto indexD1 = ::std::get<1>(alldist[j]);
					AtomIndex indexD = 0;
					if (indexD0 == indexB || indexD1 == indexB)
						continue;
					if (indexD0 == indexA) {
						auto tor = PointType::torsionRad(fs.points[indexD1], fs.points[indexA], fs.points[indexB], fs.points[indexC]);
						res.emplace_back(indexD1, indexA, indexB, indexC, tor);
						continue;
					}
					if (indexD1 == indexA) {
						auto tor = PointType::torsionRad(fs.points[indexD0], fs.points[indexA], fs.points[indexB], fs.points[indexC]);
						res.emplace_back(indexD0, indexA, indexB, indexC, tor);
						continue;
					}
					if (indexD0 == indexC) {
						auto tor = PointType::torsionRad(fs.points[indexA], fs.points[indexB], fs.points[indexC], fs.points[indexD1]);
						res.emplace_back(indexA, indexB, indexC, indexD1, tor);
						continue;
					}
					if (indexD1 == indexC) {
						auto tor = PointType::torsionRad(fs.points[indexA], fs.points[indexB], fs.points[indexC], fs.points[indexD0]);
						res.emplace_back(indexA, indexB, indexC, indexD0, tor);
						continue;
					}
				}
			}
			return ::std::make_tuple(::std::move(alldist), ::std::move(allang), ::std::move(res));
		}
	private:
		inline bool inRange(const MinMaxType pa, const FloatingPointType f) const noexcept {
			return ((f > pa.first) && (f < pa.second));
		}
	};
}