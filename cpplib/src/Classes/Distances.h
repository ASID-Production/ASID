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
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/Geometry.h"
#include "../Classes/Singleton.h"

namespace cpplib {

	/// <summary>
	/// Singleton holder for bond length limits loaded from BondLength.ini.
	/// Implements the standardized cpplib singleton interface (SingletonBase):
	/// the only instance is created lazily via Distances::instance(relpath),
	/// never eagerly at program startup; the RELATIVE path is resolved against
	/// the starting module directory (get_module_directory()).
	/// Copy/move construction and assignment are deleted (inherited from
	/// SingletonBase) to guarantee uniqueness.
	/// </summary>
	class Distances : public SingletonBase<Distances> {
		// Order of values
		// 1/1,1/2,1/3,1/4,1/5, 5/5, 2/2,2/3,2/4,2/5, 4/4,4/5, 3/3,3/4,3/5
	public:
		static constexpr size_t MAX_TYPE = cpplib::constants::mend_size - 1;

		using FloatingPointType = basic_types::FloatingPointType;
		using size_type = int_fast8_t;
		using AtomTypeBase = basic_types::AtomTypeBase;
		using PointType = geometry::Point<FloatingPointType>;
		using DataArray = ::std::array<::std::array<::std::array<FloatingPointType, 2>, MAX_TYPE + 1>, MAX_TYPE + 1>;
		using AI = basic_types::AtomIndex;

	private:
		DataArray data_{}; // [i][j][0] = min, [i][j][1] = max
		bool isReady_ = false;

		// Reference to the STATE owned by SingletonBase (atomic there). Bound in the
		// ctor init list BEFORE any read of it below — safe for reload(), which
		// rebinds the reference to its own local variable.
		::std::reference_wrapper<::std::atomic<STATE>> state_;
		// Self-contained target for the fallback empty instance (never updated).
		::std::atomic<STATE> state_holder_{STATE::INVALID_ARGUMENT};

		// Only the standardized singleton interface (SingletonBase<Distances>) may
		// create the object. The constructor never throws: it reports the outcome
		// through the output STATE reference (OK / FILE_NOT_FOUND / BAD_SIZE / ...).
		Distances(const std::filesystem::path& filename, STATE& state)
			: SingletonBase(filename), data_{}, state_(::std::ref(state)) {
			state_.get() = STATE::OK;
			::std::ifstream in(filename);
			if (!in.is_open()) {
				state_.get() = STATE::FILE_NOT_FOUND;
				return;
			}
			int mt_temp;
			if (!(in >> mt_temp)) {
				state_.get() = in.eof() ? STATE::BAD_SIZE : STATE::IO_ERROR;
				return;
			}
			if (mt_temp > static_cast<int>(MAX_TYPE)) {
				state_.get() = STATE::BAD_SIZE;
				return;
			}
			int i = 0;
			int j = 0;
			FloatingPointType lmin = 0.0f;
			FloatingPointType lmax = 0.0f;
			while (in >> i >> j >> lmin >> lmax) {
				if (i < 1 || j < 1 || i > static_cast<int>(MAX_TYPE) || j > static_cast<int>(MAX_TYPE)) {
					state_.get() = STATE::BAD_CONTENT;
					data_ = {};
					return;
				}
				data_[i][j] = {lmin, lmax};
				data_[j][i] = {lmin, lmax}; // symmetric storage
			}
			if (!in.eof()) { // partial record at the end of the file
				state_.get() = STATE::BAD_SIZE;
				data_ = {};
				return;
			}
			isReady_ = true;
		}

		// Internal plumbing of the standardized singleton interface: rebinds the
		// object's STATE reference (bound to a scratch variable during
		// construction) to the long-lived shared atomic. Called by SingletonBase
		// right after the object is created / reloaded.
		friend class SingletonBase<Distances>;
		void swap_state_target(::std::atomic<STATE>& target) noexcept {
			state_ = target; // rebinds the reference_wrapper itself
		}

	public:
		// Default-constructed EMPTY instance: value-initializes all data members,
		// reads nothing. Used only as the safe fallback by the singleton base when
		// no data could ever be loaded (INVALID_ARGUMENT). Regular users go through
		// Distances::instance()/try_call() — direct creation stays impossible.
		Distances() : SingletonBase(std::filesystem::path{}), data_{}, isReady_(false), state_(state_holder_) {}

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
			assert(length >= 0);

			const auto& [min, max] = data_[i][j];

			if (length > max) {
				return 0;
			}
			if (length < min) {
				return -1;
			}
			return 1;
		}

		/// <summary>
		/// Check if the distance is in the bond range
		/// </summary>
		/// <param name="i">Atom type of first atom</param>
		/// <param name="j">Atom type of second atom</param>
		/// <param name="lengthSq">Square of distance between atoms</param>
		/// <returns> 0  if r() more than max, 1  if r() is betweet min and max, -1 if r() less than min </returns>
		constexpr char isBondSq(AtomTypeBase i, AtomTypeBase j, FloatingPointType lengthSq) const noexcept {
			assert(i <= MAX_TYPE && j <= MAX_TYPE);
			assert(i > 0 && j > 0);
			assert(lengthSq >= 0);

			const auto& [min, max] = data_[i][j];

			double max_sq = max * max;
			double min_sq = min * min;

			if (lengthSq > max_sq) {
				return 0;
			}
			if (lengthSq < min_sq) {
				return -1;
			}
			return 1;
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

		std::vector<geometry::SG<FloatingPointType>::BondWithShift> filter_bond_list(
			const std::vector<geometry::SG<FloatingPointType>::BondWithShift>& bondlist,
			const ::std::vector<AtomTypeBase>& types,
			const ::std::vector<PointType>& points) const
		{
			std::vector<geometry::SG<FloatingPointType>::BondWithShift> ret;
			ret.reserve(points.size() << 3);

			for (const auto& bond : bondlist){
				const auto l1 = bond.first;
				const auto l2 = bond.second;

				char is_real_bond = isBond(types[l1], types[l2], sqrt(bond.dist_sq));

				if (is_real_bond != 0) {
					ret.push_back(bond);
				}
			}
			return ret;
		}

		template<typename Func>
		std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift> filter_bond_list(
			const std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift>& bondlist,
			const ::std::vector<AtomTypeBase>& types,
			const ::std::vector<PointType>& points,
			Func dist) const noexcept 
		{
			auto iter = ::std::begin(bondlist);

			std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift> ret;
			ret.reserve(points.size() << 3);

			while (iter != ::std::end(bondlist)) {
				const auto l1 = iter->first;
				const auto l2 = iter->second;
				const auto shiftcode = iter->shiftcode;

				PointType moved_point2 = points[l2];
				if (shiftcode.get_code() != 13) {
					moved_point2 += shiftcode.get_shift();
				}

				char is_real_bond = isBond(types[l1], types[l2], dist(points[l1], moved_point2));

				if (is_real_bond != 0) {
					ret.push_back(*iter);
				}

				iter++;
			}
			return ret;
		}
	};

// Standardized interface parameter: default relative data path for this singleton.
template <> inline constexpr std::string_view default_relpath<Distances> = "module/BondLength.ini";
}