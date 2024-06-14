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
#include <vector>
#include <tuple>

namespace cpplib {
	class FindGeometry {
	public:
		// Definitions
		using AtomType = currents::AtomTypeData;
		using FloatingPointType = currents::FloatingPointType;
		using AtomIndex = currents::AtomIndex;

		using size_type = FAM_Struct::size_type;
		using tupleDistance = ::std::tuple<AtomIndex, AtomIndex, FloatingPointType>;
		using tupleAngle = ::std::tuple<AtomIndex, AtomIndex, AtomIndex, FloatingPointType>;
		using tupleTorsion = ::std::tuple<AtomIndex, AtomIndex, AtomIndex, AtomIndex, FloatingPointType>;
		using PointType = geometry::Point<FloatingPointType>;
		using MinMaxType = ::std::pair<FloatingPointType, FloatingPointType>; 


	private:
		// Data
		const FAM_Struct& fs;

	public:
		FindGeometry() = delete;
		constexpr explicit FindGeometry(const FAM_Struct& famstr) noexcept : fs(famstr) {};
		auto findDistance(AtomType t1, AtomType t2, MinMaxType d12) const {
			std::vector<tupleDistance> res;
			const bool mirror = (t1 == t2);
			for (size_type i = 0; i < fs.sizePoints; i++)
			{
				if (fs.types[i] != t1) continue;
				const size_type start = mirror ? i + 1 : 0;
				for (size_type j = start; j < fs.sizePoints; j++)
				{
					if (i == j) continue;
					if (fs.types[j] != t2) continue;
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
	private:
		inline bool inRange(const ::std::pair<FloatingPointType, FloatingPointType> pa, const FloatingPointType f) const noexcept {
			return ((f > pa.first) && (f < pa.second));
		}
	};
}