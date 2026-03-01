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
#include <algorithm>
#include <array>
#include <string>
#include <type_traits>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../BaseHeaders/Concepts.h"
#include "../BaseHeaders/DebugMes.h"

/**
	 * Represents a simple atom type that stores a single AtomTypeBase identifier.
	 *
	 * Stores one AtomTypeBase value and provides queries and conversions for that
	 * identifier: `get_bitset()` yields a TypeBitset with the atom's bit set when the
	 * identifier is greater than zero, `contains(t)` reports whether the stored
	 * identifier equals `t`, `intersect`/`operator==` test equality with another
	 * SimpleAtom, and an explicit conversion returns the underlying AtomTypeBase.
	 */
	namespace cpplib {

	class SimpleAtom {
	private:
		using ConstRef = SimpleAtom;
	public:
		using AtomTypeBase = basic_types::AtomTypeBase;
		using TypeBitset = basic_types::TypeBitset;

		// Constructors
		SimpleAtom() = default;
		constexpr explicit SimpleAtom(AtomTypeBase input) : type_(input) {
			assert(input < TypeBitset().size());
		}

		constexpr TypeBitset get_bitset() const noexcept {
			TypeBitset bits;
			if (type_ > 0) bits.set(type_);
			return bits;
		}
		constexpr bool contains(const AtomTypeBase t) const noexcept {
			assert(t > 0 && TypeBitset().size());
			return type_ == t;
		}
		// Same as operator==
		constexpr bool intersect(ConstRef other) const {
			return type_ == other.type_;
		}

		// operators
		constexpr bool operator==(ConstRef other) const noexcept {
			return type_ == other.type_;
		}
		constexpr std::strong_ordering operator<=>(ConstRef other) const noexcept {
			return type_ <=> other.type_;
		}

		// Converts to AtomTypeBase
		constexpr explicit operator AtomTypeBase() const {
			return type_;
		}
	private:
		AtomTypeBase type_ = 0;
	};

	static_assert(AtomTypeConcept<SimpleAtom>, "SimpleAtom must satisfy AtomTypeConcept");

	/**
		 * Add an atom type to this composite atom's type set.
		 * @param t Atom type index to add; must be greater than 0 and less than mend_size.
		 */
		
		/**
		 * Checks whether this composite atom includes the specified atom type.
		 * @param t Atom type index to check; must be greater than 0 and less than the bitset size.
		 * @returns `true` if the bit for `t` is set in the internal type bitset, `false` otherwise.
		 */
		
		/**
		 * Obtain the bitset representing all atom types contained by this composite atom.
		 * @returns The `TypeBitset` with bits set for each included atom type.
		 */
		
		/**
		 * Determine whether this composite atom and another share at least one atom type.
		 * @param other CompositeAtom to test against.
		 * @returns `true` if the internal bitsets have any common set bit, `false` otherwise.
		 */
		
		/**
		 * Compare two CompositeAtom instances for equality based on their basetype value.
		 * @param other CompositeAtom to compare with.
		 * @returns `true` if both instances have the same `basetype`, `false` otherwise.
		 */
		
		/**
		 * Perform a three-way comparison of `basetype` values between two CompositeAtom instances.
		 * @param other CompositeAtom to compare with.
		 * @returns A `std::strong_ordering` value reflecting the ordering of the two `basetype` values.
		 */
		class CompositeAtom {
	private:
		using ConstRef = const CompositeAtom&;
	public:
		using AtomTypeBase = basic_types::AtomTypeBase;
		using TypeBitset = basic_types::TypeBitset;

		// Constructors
		CompositeAtom() = default;
		constexpr explicit CompositeAtom(AtomTypeBase input) : basetype(input) {
			if (input > 0) {
				types.set(input);
			}
		}

		void AddType(const AtomTypeBase t) {
			assert(t > 0);
			assert(t < constants::mend_size);
			types.set(t);
		}
		constexpr bool contains(const AtomTypeBase t) const {
			assert(t > 0 && t < types.size());
			return types[t];
		}
		constexpr TypeBitset get_bitset() const noexcept {
			return types;
		}
		inline bool intersect(ConstRef other) const noexcept {
			return (types & other.types).any();
		}
		// operators
		constexpr bool operator==(ConstRef other) const noexcept {
			return basetype == other.basetype;
		}
		constexpr std::strong_ordering operator<=>(ConstRef other) const noexcept {
			return basetype <=> other.basetype;
		}

		/**
		 * Yield the stored base atom type identifier for this CompositeAtom.
		 *
		 * @returns `AtomTypeBase` value representing the CompositeAtom's base type.
		 */
		constexpr explicit operator AtomTypeBase() const noexcept {
			return basetype;
		}

	private:
		AtomTypeBase basetype = 0;
		TypeBitset types = {0};
	};

	static_assert(AtomTypeConcept<CompositeAtom>, "SimpleAtom must satisfy AtomTypeConcept");

	class Coord {
	public:
		using innerType = int_fast8_t;
		using argumentType = innerType;
	private:
		innerType low = 0;
		innerType high = 0;
	public:
		constexpr Coord() noexcept = default;
		/**
		 * Construct a Coord representing a degenerate interval whose low and high bounds are the same.
		 * @param mono Value used to initialize both `low` and `high`.
		 */
		constexpr explicit Coord(argumentType mono) noexcept : low(mono), high(mono) {
		}
		constexpr Coord(argumentType first, argumentType second) noexcept : low(first), high(second) {
		}
		inline bool intersect(const Coord other) const noexcept {
			return first() <= other.second() && other.first() <= second();
		}
		inline argumentType getLow() const {
			return first();
		}
		inline argumentType getHigh() const {
			return second();
		}
	private:
		inline innerType first() const {
			return low;
		}
		inline innerType second() const {
			return high;
		}
	};

	/**
	 * Container of neighbour index shifts with a fixed maximum capacity.
	 *
	 * Stores a compact, index-shift based list of neighbour references and provides
	 * insertion, removal, replacement, simple sort, and iterator access for the
	 * active range.
	 *
	 * @note Indices and sizes are bounds-checked with assertions in debug builds.
	 */
	
	/**
	 * Represents an atom/node in the graph parameterized by an atom type.
	 *
	 * Provides node identity, atom-type, hydrogen-like count, coordinate interval,
	 * and a neighbour list (as index shifts). Supports neighbour queries, adding
	 * and removing mutual bonds, sorting neighbours, coordinate calculation, and a
	 * swap that updates neighbour references to keep internal shift values correct.
	 *
	 * @tparam A Atom type satisfying the AtomTypeConcept (e.g., SimpleAtom or CompositeAtom).
	 */
	class NeighboursType {
	public:
		static constexpr size_t maxNeighbours = constants::maxNeighbours;
		using ShiftType = basic_types::AtomIndex;
	private:
		::std::array<ShiftType, maxNeighbours> data_{0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		size_t size_ = 0;
	public:
		NeighboursType() noexcept : data_{0} {
			data_.fill(0);
		}
		inline void push_back(const ShiftType obj) {
			assert(size_ < maxNeighbours);
			data_[size_] = obj;
			size_++;
		}
		inline auto size() const {
			return size_;
		}
		inline ShiftType operator[](int8_t i) const {
			assert(i < size_);
			return data_[i];
		}
		inline ShiftType& operator[](int8_t i) {
			assert(i < size_);
			return data_[i];
		}
		void erase(int8_t i) {
			assert(i < size_);
			size_--;
			for (; i < size_; i++) {
				data_[i] = data_[i + 1];
			}
		}
		bool exchange(const ShiftType cur, const ShiftType next) {

			for (int8_t i = 0; i < size_; i++) {
				if (data_[i] == cur) {
					data_[i] = next;
					return true;
				}
			}
			return false;
		}
		constexpr void simpleSort() {
			std::sort(data_.begin(), data_.begin() + size_);
		}
		constexpr void addShift(ShiftType add) {
			for (char i = 0; i < size_; i++)
			{
				data_[i] += add;
			}
		}

		auto begin() const {
			return data_.begin();
		}
		auto end() const {
			return data_.begin() + size_;
		}
	};

	template<AtomTypeConcept A>
	class Node {
	public:
		// Declarations
		using NeighbourValueType = Node*;
		using HType = basic_types::HType;
		using AtomIndex = basic_types::AtomIndex;
		using ShiftType = NeighboursType::ShiftType;
		using AtomType = A;

		template <AtomTypeConcept X>
		friend class Node;

		friend constexpr bool ExactCompare(const Node<CompositeAtom>& a, const Node< SimpleAtom>& b);
		friend constexpr bool NotExactCompare(const Node<CompositeAtom>& a, const Node< SimpleAtom>& b);

	private:
		// Data
		NeighboursType neighbours_{};
		AtomType type_{0};
		HType hAtoms_ = 0;
		AtomIndex id_ = 0;
		Coord coord_{0};
	public:
		// Constructors
		constexpr Node() = default;
		constexpr Node(const AtomType& t1, const HType& h1, const AtomIndex& id)
			: type_(t1), hAtoms_(h1), id_(id), coord_(h1) {
		}

		// Operators
		inline bool operator==(const Node& other) const noexcept {
			return (type_ == other.type_) &&
				(hAtoms_ == other.hAtoms_) &&
				(neighbours_.size() == other.neighbours_.size()) &&
				(coord_.intersect(other.coord_));
		}
		// Raw comparision
		inline bool RawLess(const Node& other) const noexcept {
			if (type_ != other.type_)
				return type_ < other.type_;
			if (hAtoms_ + neighbours_.size() != other.hAtoms_ + other.neighbours_.size())
				return hAtoms_ + neighbours_.size() < other.hAtoms_ + other.neighbours_.size();
			if (coord_.getLow() != other.coord_.getLow())
				return coord_.getLow() < other.coord_.getLow();
			if (coord_.getHigh() != other.coord_.getHigh())
				return coord_.getHigh() < other.coord_.getHigh();
			return id_ > other.id_;
		}
		inline bool RawMore(const Node& other) const noexcept {
			if (type_ != other.type_)
				return type_ > other.type_;
			if (hAtoms_ + neighbours_.size() != other.hAtoms_ + other.neighbours_.size())
				return hAtoms_ + neighbours_.size() > other.hAtoms_ + other.neighbours_.size();
			if (coord_.getLow() != other.coord_.getLow())
				return coord_.getLow() > other.coord_.getLow();
			if (coord_.getHigh() != other.coord_.getHigh())
				return coord_.getHigh() > other.coord_.getHigh();
			return id_ < other.id_;
		}
		// simple sorting
		inline bool operator<(const Node& other) const noexcept {
			return RawLess(other);
		}
		inline bool operator>(const Node& other) const noexcept {
			return RawMore(other);
		}

		// Methods Neighbours
		constexpr bool isNeighbour(const Node& node) const noexcept(noexcept(neighbours_.operator[](0)) && noexcept(neighbours_.size())) {
			const auto s = neighbours_.size();
			const auto node_id = node.getID();
			ShiftType shift = (&node) - this;
			for (size_t i = 0; i < s; i++)
				if (neighbours_[i] == shift) return true;
			return false;
		}
		constexpr AtomIndex neighboursSize() const noexcept {
			return static_cast<AtomIndex>(neighbours_.size());
		}
		constexpr bool hasNeighbours() const noexcept {
			return neighbours_.size() != 0;
		}
		constexpr NeighbourValueType getNeighbour(AtomIndex neighbour_iterator) const {
			return const_cast<NeighbourValueType>(this + neighbours_[neighbour_iterator]);
		}

		inline AtomIndex getID() const noexcept {
			return id_;
		}
		inline void setID(const AtomIndex& id) noexcept {
			id_ = id;
		}
		inline A getType() const noexcept {
			return this->type_;
		}
		inline void setType(const A& type) noexcept {
			this->type_ = type;
		}
		inline HType getHAtoms() const noexcept {
			return this->hAtoms_;
		}
		inline void setHAtoms(const HType& hAtoms) noexcept {
			this->hAtoms_ = hAtoms;
		}
		inline Coord getCoord() const noexcept {
			return this->coord_;
		}
		inline void setCoord(const Coord& c) {
			coord_ = c;
		}

		// Algorithms
		inline void calculateCoord() {
			auto n = neighboursSize() + hAtoms_;
			coord_ = Coord(n, n);
		}
		inline void sortNeighbours() {
			neighbours_.simpleSort();
		}
		constexpr void addBondWithSort(Node& other) {
			neighbours_.push_back(&other - this);
			other.neighbours_.push_back(this - &other);

			const auto ns = neighboursSize();
			for (AtomIndex i = 0; i < ns; i++) {
				(this + neighbours_[i])->sortNeighbours();
			}

			const auto nso = other.neighboursSize();
			for (AtomIndex i = 0; i < nso; i++) {
				((&other) + other.neighbours_[i])->sortNeighbours();
			}
		}
		constexpr void addBondSimple(Node& other) {
			neighbours_.push_back(&other - this);
			other.neighbours_.push_back(this - &other);
		}
		constexpr void deleteBond(Node& other) {
			deleteNeighbour(other);
			other.deleteNeighbour(this);

			const auto ns = neighboursSize();
			for (AtomIndex i = 0; i < ns; i++) {
				(this + neighbours_[i])->sortNeighbours();
			}

			const auto nso = other.neighboursSize();
			for (AtomIndex i = 0; i < nso; i++) {
				((&other) + other.neighbours_[i])->sortNeighbours();
			}
		}
		AtomIndex findNeighbour(NeighbourValueType other) const noexcept {
			const auto s = neighboursSize();
			const ShiftType shift = other - this;
			for (AtomIndex i = 0; i < s; i++) {
				if (neighbours_[i] == shift)
					return i;
			}
			return AtomIndex(-1);
		}
		void addNeighboursVector(NeighboursType&& other) {
			neighbours_ = std::move(other);
		}
		NeighboursType getNeighboursVector() const {
			return neighbours_;

		}
		void exchangeNeighbour(const Node* cur, const Node* next) noexcept {
			auto curshift = cur - this;
			auto nextshift = next - this;
			bool ret = neighbours_.exchange(curshift, nextshift);
			assert(ret);
		}
		constexpr void swap(Node& other) noexcept {
			::std::swap(type_, other.type_);
			::std::swap(hAtoms_, other.hAtoms_);
			::std::swap(coord_, other.coord_);

			const auto pother = &other;
			this->changeNeigboursOfNeighbours(pother);
			pother->changeNeigboursOfNeighbours(this);

			::std::swap(neighbours_, other.neighbours_);
			NeighboursType::ShiftType shift = &other - this;
			neighbours_.addShift(shift);
			other.neighbours_.addShift(-shift);

			// Aditional shift if bond exists
			bool bondExists = neighbours_.exchange(0, shift);
			if (bondExists == true) other.neighbours_.exchange(0, -shift);
		}
	private:
		void changeNeigboursOfNeighbours(Node* other) noexcept {
			const auto n1size = this->neighboursSize();
			for (AtomIndex i = 0; i < n1size; i++) {
				this->getNeighbour(i)->exchangeNeighbour(this, other);
			}
		}
		void deleteNeighbour(const Node& node) noexcept {
			return deleteNeighbour(&node);
		}
		constexpr void deleteNeighbour(const Node* pnode) noexcept {
			assert(pnode != this);
			auto s = neighbours_.size();
			ShiftType shift = pnode - this;
			decltype(s) i = 0;
			for (; i < s; i++) {
				if (neighbours_[i] == shift) {
					neighbours_.erase(i);
					return;
				}
			}
		}
	};

	constexpr bool ExactCompare(const Node<CompositeAtom>& a, const Node<SimpleAtom>& b) {
		static_assert(std::is_same_v<SimpleAtom::AtomTypeBase, CompositeAtom::AtomTypeBase>);
		return a.type_.contains(static_cast<typename CompositeAtom::AtomTypeBase>(b.type_)) &&
			a.hAtoms_ == b.hAtoms_ &&
			a.neighbours_.size() == b.neighbours_.size() &&
			a.coord_.intersect(b.coord_);
	}
	/**
	 * Determines whether a composite-node pattern `a` can match a simple-node candidate `b`
	 * under relaxed (non-exact) constraints.
	 * @param a Composite-node pattern to test.
	 * @param b Simple-node candidate to test against the pattern.
	 * @returns `true` if `a`'s type contains `b`'s type, `a`'s hydrogen count is less than or equal to `b`'s,
	 * `a`'s neighbour count is less than or equal to `b`'s, and their coordinate intervals intersect; `false` otherwise.
	 */
	constexpr bool NotExactCompare(const Node<CompositeAtom>& a, const Node<SimpleAtom>& b) {
		static_assert(std::is_same_v<SimpleAtom::AtomTypeBase, CompositeAtom::AtomTypeBase>);
		return a.type_.contains(static_cast<typename CompositeAtom::AtomTypeBase>(b.type_)) &&
			a.hAtoms_ <= b.hAtoms_ &&
			a.neighbours_.size() <= b.neighbours_.size() &&
			a.coord_.intersect(b.coord_);
	}
}
/**
 * Swap the complete state of two cpplib::Node instances.
 *
 * @param n1 First node whose contents will be exchanged.
 * @param n2 Second node whose contents will be exchanged.
 */
namespace std {
	// std::swap extention for Node class
	template<cpplib::AtomTypeConcept A>
	constexpr void swap(cpplib::Node<A>& n1, cpplib::Node<A>& n2) noexcept {
		n1.swap(n2);
	}
}