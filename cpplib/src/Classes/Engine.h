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
	 * Lightweight representation of a single atom type.
	 *
	 * Stores a single atom type identifier and exposes simple queries and comparisons
	 * suitable for use where an atom type is represented by a single base value.
	 */
	
	/**
	 * Construct a SimpleAtom with the given atom type.
	 *
	 * @param input Atom type identifier to store.
	 */
	
	/**
	 * Return a bitset with the bit for this atom's type set.
	 *
	 * @returns A TypeBitset with the bit corresponding to this atom's type set if the type is greater than zero; otherwise an empty bitset.
	 */
	
	/**
	 * Check whether this atom's type equals the given type.
	 *
	 * @param t Atom type identifier to compare.
	 * @returns `true` if this atom's type equals `t`, `false` otherwise.
	 */
	
	/**
	 * Determine whether this atom and another represent the same type.
	 *
	 * @param other Other SimpleAtom to compare.
	 * @returns `true` if both atoms have the same type, `false` otherwise.
	 */
	
	/**
	 * Compare equality with another SimpleAtom.
	 *
	 * @param other Other SimpleAtom to compare.
	 * @returns `true` if both atoms have the same type, `false` otherwise.
	 */
	
	/**
	 * Three-way ordering comparison by atom type.
	 *
	 * @param other Other SimpleAtom to compare.
	 * @returns The strong ordering result of comparing atom types.
	 */
	
	/**
	 * Convert this SimpleAtom to its underlying AtomTypeBase value.
	 *
	 * @returns The stored AtomTypeBase value.
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
		 * Representation of a composite atom type that can contain multiple constituent types.
		 *
		 * Stores a primary basetype and a bitset of constituent types. Provides operations to
		 * add and query constituent types, obtain the underlying bitset, test intersection
		 * with another CompositeAtom, and compare by basetype.
		 *
		 * Methods:
		 *  - CompositeAtom(AtomTypeBase input): sets basetype and marks `input` in the bitset when > 0.
		 *  - void AddType(AtomTypeBase t): marks type `t` in the bitset; requires 0 < t < constants::mend_size.
		 *  - bool contains(AtomTypeBase t) const: returns whether type `t` is present; requires 0 < t < bitset size.
		 *  - TypeBitset get_bitset() const noexcept: returns the internal type bitset.
		 *  - bool intersect(const CompositeAtom& other) const noexcept: returns whether any constituent type is shared.
		 *  - operator== / operator<=>: compare or order by the stored basetype.
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
		 * Determine whether this coordinate interval overlaps with another interval.
		 *
		 * Intervals are inclusive: [low, high].
		 * @param other Other Coord to test for intersection.
		 * @returns `true` if the intervals overlap (have any point in common), `false` otherwise.
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
		 * Initialize the coordinate interval to a single value.
		 * @param mono Value to set for both `low` and `high`.
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
	 * Fixed-capacity container of neighbour index shifts for a node.
	 *
	 * Stores up to `maxNeighbours` neighbour shifts in a compact array and tracks the current active size.
	 *
	 * @note All index/shift arguments are interpreted as relative shifts used by Node pointers; callers must ensure values are valid and within `maxNeighbours`.
	 *
	 * @param push_back.obj Shift value to append to the active neighbour list; appends at the current end.
	 * @returns size The current number of active neighbours.
	 *
	 * Member methods:
	 * - push_back(const ShiftType obj): append a neighbour shift (asserts capacity).
	 * - size() const: return the number of active neighbours.
	 * - operator[](int8_t i) / operator[](int8_t i) const: access the i-th active shift (asserts i < size()).
	 * - erase(int8_t i): remove the element at index i and shift subsequent elements left.
	 *   @param i Index of the element to erase (must be < size()).
	 * - exchange(const ShiftType cur, const ShiftType next): replace the first occurrence of `cur` with `next`.
	 *   @param cur Shift value to find and replace.
	 *   @param next Shift value to write in place of `cur`.
	 *   @returns `true` if a replacement occurred, `false` otherwise.
	 * - simpleSort(): sort the active range of shifts in ascending order.
	 * - addShift(ShiftType add): add `add` to each active shift (used to adjust shifts after node relocation).
	 *   @param add Value added to each active shift (may be negative).
	 * - begin() / end(): iterators over the underlying array; `end()` points to `begin() + size()`.
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
	 * Determines whether composite node `a` non-strictly matches simple node `b`.
	 *
	 * @param a Composite-node candidate that may contain multiple atom types.
	 * @param b Simple-node candidate with a single atom type.
	 * @returns `true` if `a`'s type contains `b`'s type, `a.hAtoms_` <= `b.hAtoms_`, `a.neighbours_.size()` <= `b.neighbours_.size()`, and `a.coord_` intersects `b.coord_`; `false` otherwise.
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
 * Exchange all core members and neighbour references between two Node objects.
 *
 * @param n1 First node to swap.
 * @param n2 Second node to swap.
 */
namespace std {
	// std::swap extention for Node class
	template<cpplib::AtomTypeConcept A>
	constexpr void swap(cpplib::Node<A>& n1, cpplib::Node<A>& n2) noexcept {
		n1.swap(n2);
	}
}