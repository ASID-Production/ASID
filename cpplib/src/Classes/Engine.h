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
#include <type_traits> // for std::conditional
#include <vector>
#include <string>

#include <bitset> // for std::array in XAtom

#include "../BaseHeaders/Support.h"
#include "../BaseHeaders/DebugMes.h"
#include "../BaseHeaders/Currents.h"
namespace cpplib {
	class XAtom {
	public:
		using SimpleAtomType = currents::AtomTypeData;
	private:
		SimpleAtomType simple_representation = 0;
		currents::TypeBitset types = {0};
	public:
		// Constructors
		XAtom() = default;
		constexpr XAtom(SimpleAtomType input) : simple_representation(input) {
			if (input > 0) {
				types.set(input);
			}
		}

		inline void AddType(const SimpleAtomType t) {
			_ASSERT(t > 0);
			_ASSERT(t < mend_size);
			types.set(t);
		}
		inline bool include(const SimpleAtomType t) const {
			_ASSERT(t > 0 && t < types.size());
			return types.test(t);
		}
		//inline bool include(const XAtom& t) const {
		//	return !((types ^ t.types) & t.types).any();
		//}
		inline bool simple_eq(const SimpleAtomType other) const {
			return simple_representation == other;
		}
		inline bool simple_eq(const XAtom& other) const {
			return simple_representation == other.simple_representation;
		}
		inline const currents::TypeBitset& get_bitset() const {
			return types;
		}
		inline SimpleAtomType get_simple() const {
			return simple_representation;
		}
		inline bool intersect(const XAtom& other) const {
			return (types & other.types).any();
		}

		// operators
		inline bool operator==(const currents::AtomTypeData other) const noexcept {
			return include(other);
		}
		// operator for sorting
		inline bool operator<(const XAtom& other) const noexcept {
			return simple_representation < other.simple_representation;
		}
		// operator for sorting
		inline bool operator>(const XAtom& other) const noexcept {
			return simple_representation > other.simple_representation;
		}
		// operator for sorting
		inline bool operator==(const XAtom& other) const noexcept {
			return simple_eq(other);
		}
		// operator for sorting
		inline bool operator!=(const XAtom& other) const noexcept {
			return !simple_eq(other);
		}
		inline explicit operator SimpleAtomType() const { return simple_representation; }
	};

	class Coord {
	public:
		using innerType = int_fast8_t;
		using argumentType = innerType;
		static constexpr innerType max = 100;
	private:
		innerType low = 0;
		innerType high = 0;
	public:
		constexpr inline Coord() noexcept {};
		constexpr inline Coord(argumentType mono) noexcept { low = mono; high = mono; };
		constexpr inline Coord(argumentType first, argumentType second) noexcept { low = first; high = second; };
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

	class NeighboursType {
	public:
		static constexpr int8_t maxNeighbours = 100;
		using ShiftType = currents::AtomIndex;
	private:
		::std::array<ShiftType, maxNeighbours> data_{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
		int8_t size_ = 0;
	public:
		NeighboursType() noexcept : data_{0} { data_.fill(0); }
		inline void push_back(const ShiftType obj) {
			_ASSERT(size_ < maxNeighbours);
			data_[size_] = obj;
			size_++;
		}
		inline auto size() const {
			return size_;
		}
		inline ShiftType operator[](int8_t i) const {
			_ASSERT(i < size_);
			return data_[i];
		}
		inline ShiftType& operator[](int8_t i) {
			_ASSERT(i < size_);
			return data_[i];
		}
		void erase(int8_t i) {
			_ASSERT(i < size_);
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
		inline void simpleSort() {
			std::sort(&(data_[0]), &(data_[size_]));
		}
		void addShift(ShiftType add) {
			for (char i = 0; i < size_; i++)
			{
				data_[i] += add;
			}
		}
	};

	template<class A>
	class Node {
	public:
		// Declarations
		using NeighbourValueType = Node*;
		using HType = currents::HType;
		using AtomIndex = currents::AtomIndex;

		template <class X>
		friend class Node;
	private:
		// Data
		NeighboursType neighbours_; // should be first
		A type_ = 0;
		HType hAtoms_ = 0;
		AtomIndex id_ = 0;
		Coord coord_ = 0;
	public:
		// Constructors
		constexpr Node(): neighbours_() {}
		constexpr Node(const A& t1, const HType& h1, const AtomIndex& id)
			: type_(t1), hAtoms_(h1), id_(id), coord_(h1), neighbours_() {}

		// Operators
		template <class X> 
		inline bool operator==(const Node<X>& other) const noexcept {
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
			//return id_ > other.id_;
			return RawLess(other);
		}
		inline bool operator>(const Node& other) const noexcept {
			//return id_ < other.id_;
			return RawMore(other);
		}
		inline bool notExactCompare(const Node<currents::AtomTypeData>& other) const noexcept {
			return type_ == other.type_ &&
				hAtoms_ <= other.hAtoms_ &&
				neighbours_.size() <= other.neighbours_.size() &&
				coord_.intersect(other.coord_);
		}

		// Methods Neighbours
		constexpr bool isNeighbour(const Node& node) const noexcept(noexcept(neighbours_.operator[](0)) && noexcept(neighbours_.size())) {
			const auto s = neighbours_.size();
			const auto node_id = node.getID();
			using ShiftType = NeighboursType::ShiftType;
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
			for (AtomIndex i = 0; i < s; i++) {
				if (neighbours_[i] == other)
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
			_ASSERT(ret);
		}
		constexpr void swap(Node& other) noexcept {
			::std::swap(type_, other.type_);
			::std::swap(hAtoms_, other.hAtoms_);
			::std::swap(coord_, other.coord_);
			const auto n1size = this->neighboursSize();
			const auto n2size = other.neighboursSize();
			for (AtomIndex i = 0; i < n1size; i++) {
				this->getNeighbour(i)->exchangeNeighbour(this, &other);
			}
			for (AtomIndex i = 0; i < n2size; i++) {
				other.getNeighbour(i)->exchangeNeighbour(&other, this);
			}
			::std::swap(neighbours_, other.neighbours_);
			NeighboursType::ShiftType shift = &other - this;
			neighbours_.addShift(shift);
			other.neighbours_.addShift(-shift);

			// Aditional shift if bond exists
			bool bondExists = neighbours_.exchange(0, shift);
			if(bondExists == true) other.neighbours_.exchange(0, -shift);
		}
	private:
		inline void deleteNeighbour(const Node& node) noexcept {
			return deleteNeighbour(&node);
		}
		constexpr void deleteNeighbour(const Node* pnode) noexcept {
			_ASSERT(pnode != this);
			using ShiftType = NeighboursType::ShiftType;
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

	struct Bond {
	public:
		// Declarations
		using AtomIndex = currents::AtomIndex;

		// Data
		AtomIndex first = 0;
		AtomIndex second = 0;

		// Constructors
		constexpr Bond() noexcept = default;
		constexpr Bond(const AtomIndex a1, const AtomIndex a2) noexcept : first(a1), second(a2) {};

		// Operators
		constexpr auto operator==(const Bond& other) const noexcept {
			return first == other.first &&
				second == other.second;
		}
		constexpr auto operator!=(const Bond& other) const noexcept {
			return first != other.first ||
				second != other.second;
		}
		constexpr auto operator<(const Bond& other) const noexcept {
			if (first != other.first)
				return first < other.first;
			return second < other.second;
		}

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
		using LengthType = currents::FloatingPointType; // Float or Double
		using base = Bond; // Bond

		// Data
		LengthType length = 0.0; // Bond Length 

		// Constructors

		constexpr BondEx() = default;
		constexpr BondEx(const base& bond, float len) noexcept : base(bond), length(len) {
			base::validate();
		}
		constexpr BondEx(AtomIndex a1, AtomIndex a2, LengthType l) noexcept
			: length(l) {
			if (a1 < a2) {
				base::first = a1;
				base::second = a2;
			}
			else {
				base::first = a2;
				base::second = a1;
			}
		}

		// Compares only "base". Ignores length.
		constexpr bool operator==(const BondEx& other) const noexcept {
			return base::operator==(other);
		}

		// Compares in next order:
		// 1. "base" 
		// 2. length
		constexpr bool operator<(const BondEx& other) const noexcept {
			if (base::operator!=(other))
				return base::operator<(other);
			return length < other.length;
		}

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
}

namespace std {
	// std::swap extention for Node class
	template<class A>
	constexpr void swap(cpplib::Node<A>& n1, cpplib::Node<A>& n2) noexcept {
		n1.swap(n2);
	}
}