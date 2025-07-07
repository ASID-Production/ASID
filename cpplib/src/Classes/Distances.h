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
#include <fstream>
#include <string>
#include <vector>
#include "../BaseHeaders/Currents.h"

namespace cpplib {
	class Distances {
		// Order of values
		// 1/1,1/2,1/3,1/4,1/5, 5/5, 2/2,2/3,2/4,2/5, 4/4,4/5, 3/3,3/4,3/5
	public:
		struct DistancesException : public ::std::runtime_error { 
			using ::std::runtime_error::runtime_error;
		};
		using FloatingPointType = currents::FloatingPointType;
		using base = ::std::vector<FloatingPointType>;
		using size_type = int_fast8_t;
		using AtomTypeBase = currents::AtomTypeBase;
	private:
		using internal_size_type = int_fast16_t;
		static_assert (sizeof(size_type) * 2 <= sizeof(internal_size_type), "Internal_size_type should be at least 2 times bigger than size_type");

		bool isReady_ = false;
		size_type maxType_ = 0;
	public:
		Distances() = delete;
		Distances(Distances&&) = delete;
		Distances(const Distances&) = delete;

		explicit Distances(const ::std::string& filename) {
			::std::ifstream in(filename);
			int mt_temp;
			if (!(in >> mt_temp))
				throw DistancesException("Error while reading file " + filename);
			maxType_ = mt_temp;
			data_.assign(maxType_ * (maxType_ + 1), 0.0f);
			int i = 0;
			int j = 0;
			FloatingPointType lmin = 0.0f;
			FloatingPointType lmax = 0.0f;
			while (in >> i) {
				if ((((i == 0) || !(in >> j)) || !(in >> lmin)) || !(in >> lmax))
					return;
				if (i > j) ::std::swap(i, j);
				auto t = indexBond(i, j);
				data_[t] = lmin;
				data_[t + 1] = lmax;
			}
			isReady_ = true;
		}
		inline bool isReady() const {
			return isReady_;
		}

		inline char isBond(AtomTypeBase i, AtomTypeBase j, const FloatingPointType length) const noexcept {
			_ASSERT(i <= maxType_ && j <= maxType_);
			_ASSERT(i > 0 && j > 0);
			if (i > j) 
				::std::swap(i, j);
			internal_size_type n = indexBond(i, j);
			if (length < data_[n + 1]) {
				if (data_[n] < length) return 1; // Usual bond
				else return -1; // Invalid bond
			}
			else return 0; // Not a bond
		}
		inline FloatingPointType minDistance(AtomTypeBase a1, AtomTypeBase a2) const noexcept {
			_ASSERT(a1 <= maxType_ && a2 <= maxType_);
			_ASSERT(a1 > 0 && a2 > 0);
			if (a1 <= a2) return data_[indexBond(a1, a2)];
			else return data_[indexBond(a2, a1)];
		}
		inline FloatingPointType maxDistance(AtomTypeBase a1, AtomTypeBase a2) const noexcept {
			_ASSERT(a1 <= maxType_ && a2 <= maxType_);
			_ASSERT(a1 > 0 && a2 > 0);
			if (a1 <= a2) return data_[indexBond(a1, a2) + 1];
			else return data_[indexBond(a2, a1) + 1];
		}
	private:
		constexpr internal_size_type indexBond(size_type i, size_type j) const {
			auto a1 = ::std::min(i - 1, maxType_ - i);
			if (i - a1 == 1) {
				return (a1 * (maxType_ + 1) + j - i) << 1;
			}
			else {
				return ((a1 + 1) * (maxType_ + 1) - 1 - j + i) << 1;
			}
		}
	private:
		// Data
		std::vector<FloatingPointType> data_;
	};
}