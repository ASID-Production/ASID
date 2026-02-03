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

namespace cpplib {
	/// @brief Represents a bond between two atoms identified by their indices
	///
	/// A Bond stores two atom indices in canonical form (first <= second after validate()).
	/// Bonds are comparable and support validation to ensure consistent ordering.
	/// This is the base class for bond representations in the cpplib library.
	struct Bond {
	public:
		/// @brief Alias for atom index type
		using AtomIndex = basic_types::AtomIndex;

		/// @brief Index of the first atom in the bond
		AtomIndex first = 0;
		/// @brief Index of the second atom in the bond
		AtomIndex second = 0;

		/// @brief Default constructor, initializes both indices to 0
		constexpr Bond() noexcept = default;

		/// @brief Construct a bond from two atom indices
		/// @param a1 First atom index
		/// @param a2 Second atom index
		constexpr Bond(const AtomIndex a1, const AtomIndex a2) noexcept : first(a1), second(a2) {};

		// Operators
		/// @brief Three-way comparison operator for bond ordering
		/// @param other Bond to compare with
		/// @return Comparison result (less, equal, or greater)
		constexpr auto operator<=>(const Bond& other) const noexcept = default;

		/// @brief Validates and normalizes bond ordering (ensures first <= second)
		///
		/// This function should be called after construction or modification to ensure
		/// the bond indices are in canonical form.
		constexpr void validate() noexcept {
			if (first > second) ::std::swap(first, second);
		}

		/// @brief Convert bond to string representation
		/// @return String in format "(first, second)", e.g., "(1, 2)"
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

	/// @brief Extended bond with length information
	///
	/// BondEx extends Bond by adding a length field to store the physical distance
	/// between bonded atoms. Inherits validation and ordering from Bond.
	/// The equality operator compares only the base Bond (atom indices), while the
	/// three-way comparison includes the length for complete ordering.
	struct BondEx : public Bond {
	public:
		// Declarations
		/// @brief Alias for floating point type used for length
		using LengthType = basic_types::FloatingPointType;
		/// @brief Alias for base Bond type
		using base = Bond; // Bond

		// Data
		/// @brief Physical length/distance of the bond
		LengthType length = 0.0; // Bond Length 

		// Constructors
		/// @brief Default constructor
		constexpr BondEx() = default;

		/// @brief Construct from base Bond and length
		/// @param bond Base bond with atom indices
		/// @param len Bond length
		///
		/// Always uses Bond::validate() to ensure canonical ordering.
		constexpr BondEx(const base& bond, float len) noexcept : base(bond), length(len) {
			base::validate();
		}

		/// @brief Construct from two atoms and length
		/// @param a1 First atom index
		/// @param a2 Second atom index
		/// @param l Bond length
		///
		/// Automatically orders indices to ensure canonical form (smaller index first).
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

		/// @brief Construct from two atoms (length = 0)
		/// @param a1 First atom index
		/// @param a2 Second atom index
		constexpr BondEx(AtomIndex a1, AtomIndex a2) noexcept : Bond(a1, a2) {}

		/// @brief Equality comparison (compares only base Bond, not length)
		/// @param other BondEx to compare with
		/// @return True if base bonds (atom indices) are equal
		///
		/// Note: This comparison ignores the length field for equality testing.
		constexpr bool operator==(const BondEx& other) const noexcept {
			return base::operator==(other);
		}

		/// @brief Three-way comparison (includes length)
		/// @param other BondEx to compare with
		/// @return Comparison result
		///
		/// Compares in order:
		/// 1. Base bond (atom indices)
		/// 2. Length (if base bonds are equal)
		constexpr auto operator<=>(const BondEx& other) const noexcept = default;

		/// @brief Convert bond to string representation with distance
		/// @return String in format "(first, second, {"distance": length})"
		///
		/// Example output: "(1, 2, {"distance": 1.5})"
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


	/// @brief Bond with an associated shift/offset point
	/// @tparam PT Point type for storing the shift vector
	///
	/// BondWithPoint extends Bond by adding a shift field to track periodic boundary
	/// crossings or spatial offsets associated with the bond in periodic systems.
	/// The shift typically represents which unit cell image the second atom is in
	/// relative to the first atom.
	template<class PT>
	struct BondWithPoint : public Bond {
	public:
		/// @brief Default constructor
		constexpr BondWithPoint() = default;

		/// @brief Inherit constructors from Bond
		using Bond::Bond;

		/// @brief Alias for the point/shift type
		using PointType = PT;

		// Data
		/// @brief Shift or offset vector associated with this bond
		///
		/// In periodic systems, this represents the unit cell translation needed
		/// to correctly position the second atom relative to the first.
		PointType shift{};

		// Operators
		/// @brief Equality comparison
		/// @param other BondWithPoint to compare with
		/// @return True if both bond indices are equal (ignores shift)
		///
		/// Note: Only compares the base Bond (atom indices), not the shift.
		constexpr bool operator==(const BondWithPoint& other) const noexcept {
			return Bond::operator==(other);
		}

		/// @brief Three-way comparison operator
		/// @param other BondWithPoint to compare with
		/// @return Comparison result
		///
		/// Compares both the base bond and the shift for complete ordering.
		constexpr auto operator<=>(const BondWithPoint& other) const noexcept = default;
	};

	/// @brief Compile-time check: Bond satisfies BondConcept
	static_assert(BondConcept<Bond>, "Bond must satisfy BondConcept");
	/// @brief Compile-time check: BondEx satisfies BondConcept
	static_assert(BondConcept<BondEx>, "BondEx must satisfy BondConcept");
	/// @brief Compile-time check: BondWithPoint satisfies BondConcept
	static_assert(BondConcept<BondWithPoint<std::array<char, 3>>>, "BondWithShift must satisfy BondConcept");

} // namespace cpplib