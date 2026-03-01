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
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/Geometry.h"

namespace cpplib {
	class Distances {
		// Order of values
		// 1/1,1/2,1/3,1/4,1/5, 5/5, 2/2,2/3,2/4,2/5, 4/4,4/5, 3/3,3/4,3/5
	public:
		static constexpr size_t MAX_TYPE = cpplib::constants::mend_size - 1;

		struct DistancesException : public ::std::runtime_error {
			using ::std::runtime_error::runtime_error;
		};

		using FloatingPointType = basic_types::FloatingPointType;
		using size_type = int_fast8_t;
		using AtomTypeBase = basic_types::AtomTypeBase;
		using PointType = geometry::Point<FloatingPointType>;
		using DataArray = ::std::array<::std::array<::std::array<FloatingPointType, 2>, MAX_TYPE + 1>, MAX_TYPE + 1>;
		using AI = basic_types::AtomIndex;
	private:

		DataArray data_{}; // [i][j][0] = min, [i][j][1] = max
		bool isReady_ = false;
	public:
		Distances() = delete;
		Distances(Distances&&) = delete;
		Distances(const Distances&) = delete;

		explicit Distances(const ::std::string& filename) {
			::std::ifstream in(filename);
			int mt_temp;
			if (!(in >> mt_temp))
				throw DistancesException("Error while reading file " + filename);
			if (mt_temp > MAX_TYPE) throw DistancesException("MAX_TYPE is greater than " + ::std::to_string(MAX_TYPE));
			int i = 0;
			int j = 0;
			FloatingPointType lmin = 0.0f;
			FloatingPointType lmax = 0.0f;
			while (in >> i >> j >> lmin >> lmax) {
				if (i < 1 || j < 1 || i > MAX_TYPE || j > MAX_TYPE)
					throw DistancesException("File contains incorrect data.");
				data_[i][j] = {lmin, lmax};
				data_[j][i] = {lmin, lmax}; // symmetric storage
			}
			isReady_ = true;
		}
		inline bool isReady() const {
			return isReady_;
		}

		/// <summary>
		/// Check if the distance is in the bond range
		/// </summary>
		/// <param name="i">Atom type of first atom</param>
		/// <param name="j">Atom type of second atom</param>
		/// <param name="length">Distance between atoms</param>
		/// <returns> 0  if r() more than max, 1  if r() is betweet min and max, -1 if r() less than min </returns>
		constexpr char isBond(AtomTypeBase i, AtomTypeBase j, FloatingPointType length) const noexcept {
			assert(i <= MAX_TYPE && j <= MAX_TYPE);
			assert(i > 0 && j > 0);
			const auto& [min, max] = data_[i][j];
			if (length < max) {
				if (min < length) {
					return 1;
				} else {
					return -1;
				}
			} else {
				return 0;
			}
		}
		inline FloatingPointType minDistance(AtomTypeBase a1, AtomTypeBase a2) const noexcept {
			assert(a1 <= MAX_TYPE && a2 <= MAX_TYPE);
			assert(a1 > 0 && a2 > 0);
			return data_[a1][a2][0];
		}
		constexpr FloatingPointType maxDistance(AtomTypeBase a1, AtomTypeBase a2) const noexcept {
			assert(a1 <= MAX_TYPE && a2 <= MAX_TYPE);
			assert(a1 > 0 && a2 > 0);
			return data_[a1][a2][1];
		}

		template<typename Func>
		void filter_bond_list(std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift>& bondlist,
							  const ::std::vector<AtomTypeBase>& types,
							  const ::std::vector<PointType>& points,
							  Func dist) const noexcept {
			auto iter = ::std::begin(bondlist);

			while (iter != ::std::end(bondlist)) {
				const auto l1 = iter->first;
				const auto l2 = iter->second;
				const auto shiftcode = iter->shiftcode;

				PointType moved_point2 = points[l2];
				if (shiftcode.get_code() != 13) {
					moved_point2 += shiftcode.get_shift();
				}

				char is_real_bond = isBond(types[l1], types[l2], dist(points[l1], moved_point2));

				if (is_real_bond == 0) {
					iter->first = 0;
					iter->second = 0;
				}

				iter++;
			}
		}
	};
}