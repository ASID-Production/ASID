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

template <class A, class AI, class T>
class FindGeometry {
public:
	// Definitions
	using PointType = geometry::Point<T>;

private:
	// Data
	const FAM_Struct<A,AI,T>& fs;
	
public:
	FindGeometry() = delete;
	explicit FindGeometry(const FAM_Struct<A, AI, T>& famstr) noexcept : fs(famstr) {};
	auto findDistance(A t1, A t2, T d12_l, T d12_h) const {
		std::vector<std::tuple<size_t, size_t, T>> res;
		const bool mirror = t1 == t2;
		for (size_t i = 0; i < fs.sizePoints; i++)
		{
			if (fs.types[i] != t1) continue;
			const size_t start = mirror ? i+1 : 0;
			for (size_t j = start; j < fs.sizePoints; j++)
			{
				if (i == j) continue;
				if (fs.types[j] != t2) continue;
				auto R = (fs.points[i] - fs.points[j]).r();
				if (R > d12_l && R < d12_h)
					res.emplace_back(i, j, R);
			}
		}
		return res;
	}
	auto findAngle(const std::vector<std::tuple<size_t, size_t, T>>& v1, const std::vector<std::tuple<size_t, size_t, T>>& v2, T A123_l, T A123_h) const {
		auto v1s = v1.size();
		auto v2s = v2.size();
		std::vector<std::tuple<size_t, size_t, size_t, T>> res;
		if (v1s == 0 || v2s == 0)
			return res;
		// flags
		const bool is_v1_mirror = (std::get<0>(v1[0]) == std::get<1>(v1[0]));
		const bool is_v2_mirror = (std::get<0>(v2[0]) == std::get<1>(v2[0]));


		// comparation loop
		for (size_t i = 0; i < v1s; i++)
		{
			for (size_t j = 0; j < v2s; j++)
			{
				if (std::get<1>(v1[i]) != std::get<0>(v2[j]))
					continue;
				auto rad = PointType::angleRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v1[i])], fs.points[std::get<1>(v2[i])]);
				if (rad < A123_l || rad > A123_h)
					continue;
				res.emplace_back(std::get<0>(v1[i]), std::get<1>(v1[i]), std::get<1>(v2[i]), rad);
			}
		}
		if (is_v1_mirror) {
			for (size_t i = 0; i < v1s; i++)
			{
				for (size_t j = 0; j < v2s; j++)
				{
					if (std::get<0>(v1[i]) != std::get<0>(v2[j]))
						continue;
					auto rad = PointType::angleRad(fs.points[std::get<1>(v1[i])], fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v2[i])]);
					if (rad < A123_l || rad > A123_h)
						continue;
					res.emplace_back(std::get<1>(v1[i]), std::get<0>(v1[i]), std::get<1>(v2[i]), rad);
				}
			}
		}
		if (is_v2_mirror) {
			for (size_t i = 0; i < v1s; i++)
			{
				for (size_t j = 0; j < v2s; j++)
				{
					if (std::get<1>(v1[i]) != std::get<1>(v2[j]))
						continue;
					auto rad = PointType::angleRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v2[i])], fs.points[std::get<0>(v2[i])]);
					if (rad < A123_l || rad > A123_h)
						continue;
					res.emplace_back(std::get<0>(v1[i]), std::get<1>(v2[i]), std::get<0>(v2[i]), rad);
				}
			}
		}
		return res;
	}
	auto findTorsion(std::vector<std::tuple<size_t, size_t, size_t, T>> v1, std::vector<std::tuple<size_t, size_t, size_t, T>> v2, T t1234_l, T t1234_h) const {
		auto v1s = v1.size();
		auto v2s = v2.size();
		std::vector<std::tuple<size_t, size_t, size_t, size_t, T>> res;
		if (v1s == 0 || v2s == 0)
			return res;
		for (size_t i = 0; i < v1s; i++)
		{
			for (size_t j = 0; j < v2s; j++)
			{
				if (std::get<1>(v1[i]) != std::get<0>(v2[j]) || std::get<2>(v1[i]) != std::get<1>(v2[j]) || std::get<0>(v1[i]) >= std::get<2>(v2[j]))
					continue;
				auto tor = PointType::torsionRad(fs.points[std::get<0>(v1[i])], fs.points[std::get<1>(v1[i])], fs.points[std::get<1>(v2[i])], fs.points[std::get<2>(v2[i])]);
				if (tor < t1234_l || tor > t1234_h) 
					continue;
				res.emplace_back(std::get<0>(v1[i]), std::get<1>(v1[i]), std::get<1>(v2[i]), std::get<2>(v2[i]), tor);	
			}
		}
		return res;
	}
};