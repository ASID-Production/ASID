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
#include <string>

#include "../BaseHeaders/BaseTypes.h"
#include "../BaseHeaders/Concepts.h"
/**
		 * Initialize a Bond with the given atom indices.
		 * @param a1 Index of the first atom.
		 * @param a2 Index of the second atom.
		 */
		
		/**
		 * Ensure the bond's atom indices are ordered so that `first <= second`.
		 */
		
		/**
		 * Return a textual representation of the bond in the form "(first, second)".
		 * @returns A string formatted as "(<first>, <second>)".
		 */
		
		/**
		 * Default-construct a BondEx.
		 */
		
		/**
		 * Initialize a BondEx from a Bond and a length value; normalizes atom order.
		 * @param bond The base Bond providing atom indices.
		 * @param len The bond length.
		 */
		
		/**
		 * Initialize a BondEx with two atom indices and a length; atom indices are stored in ascending order.
		 * @param a1 Index of the first atom.
		 * @param a2 Index of the second atom.
		 * @param l  Bond length.
		 */
		
		/**
		 * Initialize a BondEx with two atom indices.
		 * @param a1 Index of the first atom.
		 * @param a2 Index of the second atom.
		 */
		
		/**
		 * Compare two BondEx objects by their atom indices only, ignoring `length`.
		 * @param other The BondEx to compare against.
		 * @returns `true` if the base Bond portions are equal, `false` otherwise.
		 */
		
		/**
		 * Return a textual representation of the BondEx in the form "(first, second, {"distance": length})".
		 * @returns A string formatted as "(<first>, <second>, {\"distance\": <length>})".
		 */
		namespace cpplib {
	struct Bond {
	public:
		// Declarations
		using AtomIndex = basic_types::AtomIndex;

		// Data
		AtomIndex first = 0;
		AtomIndex second = 0;

		// Constructors
		constexpr Bond() noexcept = default;
		constexpr Bond(const AtomIndex a1, const AtomIndex a2) noexcept : first(a1), second(a2) {
		};

		// Operators
		constexpr auto operator<=>(const Bond& other) const noexcept = default;

		// Functions
		constexpr void validate() noexcept {
			if (first > second) ::std::swap(first, second);
		}
		::std::string ToStr() const {
			::std::string res("(");
			res += ::std::to_string(this->first);
			res += ", ";
			res += ::std::to_string(this->second);
			res += ")";
			// (1, 2)
			return res;
		}
	};
	struct BondEx : public Bond {
	public:
		// Declarations
		using LengthType = basic_types::FloatingPointType; // Float or Double
		using base = Bond; // Bond

		// Data
		LengthType length = 0.0; // Bond Length 

		// Constructors
		constexpr BondEx() = default;

		// Always uses Bond::validate function
		constexpr BondEx(const base& bond, float len) noexcept : base(bond), length(len) {
			base::validate();
		}
		constexpr BondEx(AtomIndex a1, AtomIndex a2, LengthType l) noexcept
			: length(l) {
			if (a1 < a2) {
				base::first = a1;
				base::second = a2;
			} else {
				base::first = a2;
				base::second = a1;
			}
		}

		constexpr BondEx(AtomIndex a1, AtomIndex a2) noexcept : Bond(a1, a2) {
		}

		// Compares only "base". Ignores length.
		constexpr bool operator==(const BondEx& other) const noexcept {
			return base::operator==(other);
		}

		// Compares in next order:
		// 1. "base" 
		// 2. length
		constexpr auto operator<=>(const BondEx& other) const noexcept = default;

		// to_string constant function
		::std::string ToStr() const {
			::std::string res("(");
			res += ::std::to_string(first);
			res += ", ";
			res += ::std::to_string(second);
			res += ", {\"distance\": ";
			res += ::std::to_string(length);
			res += "})";
			return res;
			// (1, 2, {"distance": 1.0})
		}
	};


	template<class PT>
	struct BondWithPoint : public Bond {
	public:
		constexpr BondWithPoint() = default;
		using Bond::Bond;
		using PointType = PT;

		/**
 * Positional offset associated with the bond.
 *
 * Holds the point/vector value for the bond's shift; default-initialized.
 */
		PointType shift{};

		/**
		 * Determine whether two BondWithPoint objects have equal bond indices (compare only the base Bond).
		 * @param other BondWithPoint to compare against.
		 * @returns `true` if the base Bond members (`first` and `second`) are equal, `false` otherwise.
		 */
		constexpr bool operator==(const BondWithPoint& other) const noexcept {
			return Bond::operator==(other);
		}
		constexpr auto operator<=>(const BondWithPoint& other) const noexcept = default;
	};

	static_assert(BondConcept<Bond>, "Bond must satisfy BondConcept");
	static_assert(BondConcept<BondEx>, "BondEx must satisfy BondConcept");
	static_assert(BondConcept<BondWithPoint<std::array<char,3>>>, "BondWithShift must satisfy BondConcept");

} // namespace cpplib