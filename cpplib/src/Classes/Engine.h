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

class XAtom {
private:
	char simple_representation = 0;
	std::bitset<mend_size> types = {0};
public:
	// Constructors
	constexpr XAtom() = default;
	explicit constexpr XAtom(char input) : simple_representation(input) {
		if (input > 0) {
			types.set(input);
		}
	}

	inline void AddType(const char t) {
		types.set(t);
	}
	inline bool include(const char t) const {
		return types.test(t);
	}
	inline bool include(const XAtom& t) const {
		return !((types ^ t.types) & t.types).any();
	}
	inline bool simple_eq(const char other) const {
		return simple_representation == other;
	}
	inline const std::bitset<mend_size>& get_bitset() const {
		return types;
	}

	// operators
	inline bool operator==(const char other) const noexcept {
		return include(other);
	}
	inline bool operator!=(const char other) const noexcept {
		return !include(other);
	}

	// operator for sorting
	inline bool operator<(const XAtom& other) const noexcept {
		return simple_representation < other.simple_representation;
	}
	// operator for sorting
	inline bool operator<=(const XAtom& other) const noexcept {
		return simple_representation <= other.simple_representation;
	}
	// operator for sorting
	inline bool operator>(const XAtom& other) const noexcept {
		return simple_representation > other.simple_representation;
	}
	// operator for sorting
	inline bool operator>=(const XAtom& other) const noexcept {
		return simple_representation >= other.simple_representation;
	}
	// operator for sorting
	inline bool operator==(const XAtom& other) const noexcept {
		return simple_representation == other.simple_representation;
	}
	// operator for sorting
	inline bool operator!=(const XAtom& other) const noexcept {
		return simple_representation != other.simple_representation;
	}
	inline explicit operator char() const { return simple_representation; }
};

template<class A, class H>
// H - integral type
class BaseNode {
public:
	// Data
	A type;
	H hAtoms;

	template <class A1,class H1>
	friend class BaseNode;

	//Constructors
	constexpr BaseNode() = default;
	constexpr BaseNode(const A type1, const H hAtoms1) noexcept : type(type1), hAtoms(hAtoms1) {}
	constexpr BaseNode(const BaseNode<A, H>&) noexcept = default;
	constexpr BaseNode(BaseNode<A, H>&&) noexcept = default;

	// Operators
	template <class A1, class H1>
	constexpr bool operator==(const BaseNode<A1, H1>& other) const noexcept {
		return type == other.type &&
			hAtoms == other.hAtoms;
	}
	template <class A1, class H1>
	constexpr bool operator<(const BaseNode<A1, H1>& other) const noexcept {
		if (type != other.type)
			return type < other.type;
		return hAtoms < other.hAtoms;
	}
	template <class A1, class H1>
	constexpr bool operator>(const BaseNode<A1, H1>& other) const noexcept {
		if (type != other.type)
			return type > other.type;
		return hAtoms > other.hAtoms;
	}
};

template<class A, class H, class AI>
class Node : private BaseNode<A, H> {
public:
	// Declarations
	using base = BaseNode<A, H>;
	using NeighbourValueType = Node*;
	using NeighboursType = std::vector<NeighbourValueType>;

	template <class X, class H1, class AI1>
	friend class Node;
private:
	// Data
	AI id_ = 0;
	NeighboursType neighbours_;
public:
	// Constructors
	constexpr Node() : base() {
		initializeNeighbours();
	}
	constexpr Node(const base& bn, const AI& id)
		: base(bn), id_(id) {
		initializeNeighbours();
	}
	constexpr Node(const A& t1, const H& h1, const AI& id)
		: base(t1, h1), id_(id) {
		initializeNeighbours();
	}
	constexpr Node(const A& t1, const H& h1, const AI& id, NeighboursType&& neighbours)
		: base(t1, h1), neighbours_(std::move(neighbours)), id_(id) {}

	// Operators
	template <class X>
	bool operator==(const Node<X,H,AI>& other) const noexcept {
		return (((base&)(*this)) == ((const struct Node<X, H, AI>::base&)other)) &&
			(neighbours_.size() == other.neighbours_.size());
	}

	bool operator<(const Node& other) const noexcept {
		if (!((base&)(*this) == (const base&)other))
			return (base&)(*this) < (const base&)other;
		if (neighbours_.size() != other.neighbours_.size())
			return neighbours_.size() < other.neighbours_.size();
		return id_ < other.id_;
	}
	bool operator>(const Node& other) const noexcept {
		if (!((base&)(*this) == (const base&)other))
			return (base&)(*this) > (const base&)other;
		if (neighbours_.size() != other.neighbours_.size())
			return neighbours_.size() > other.neighbours_.size();
		return id_ > other.id_;
	}
	bool operator<=(const Node& other) const noexcept {
		return this->operator==(other) || this->operator<(other);
	}
	template<class X> bool notExactCompare(const Node<X,H,AI>& other) const noexcept {
		return ((const base&)(*this)).type == ((const struct Node<X, H, AI>::base&)other).type &&
			((const base&)(*this)).hAtoms <= ((const struct Node<X, H, AI>::base&)other).hAtoms &&
			(neighbours_.size() <= other.neighbours_.size());
	}

	// Methods Neighbours
	constexpr bool isNeighbour(const Node& node) const noexcept(noexcept(neighbours_.operator[](0)) && noexcept(neighbours_.size())) {
		const auto s = neighbours_.size();
		const auto node_id = node.getID();
		for (size_t i = 0; i < s; i++)
			if (neighbours_[i] == &node) return true;
		return false;
	}
	constexpr AI neighboursSize() const noexcept {
		return static_cast<AI>(neighbours_.size());
	}
	constexpr bool hasNeighbours() const noexcept {
		return !neighbours_.empty();
	}
	constexpr NeighbourValueType getNeighbour(AI neighbour_iterator) const {
		return neighbours_[neighbour_iterator];
	}

	// Interface Const ID(ref)
	constexpr const AI& getID() const noexcept {
		return id_;
	}
	// Interface Non-const ID(ref)
	constexpr void setID(const AI& id) noexcept {
		id_ = id;
	}
	// Interface Const type(ref)
	constexpr const A& getType() const noexcept {
		return this->base::type;
	}
	// Interface Non-const type(ref)
	constexpr void setType(const A& type) noexcept {
		this->base::type = type;
	}	
	// Interface Const hAtoms(ref)
	constexpr const H& getHAtoms() const noexcept {
		return this->base::hAtoms;
	}
	constexpr void setHAtoms(const H& hAtoms) noexcept {
		this->base::hAtoms = hAtoms;
	}

	// Algorithms
	constexpr void sortNeighbours() {
		std::sort(neighbours_.begin(), neighbours_.end(),
						  [](const Node* a1, const Node* a2) {
							  return (*a1) > (*a2);
						  });
	}
	constexpr void addBondWithSort(Node& other) {
		neighbours_.emplace_back(&other);
		other.neighbours_.emplace_back(this);

		const auto ns = neighboursSize();
		for (AI i = 0; i < ns; i++) {
			neighbours_[i]->sortNeighbours();
		}

		const auto nso = other.neighboursSize();
		for (AI i = 0; i < nso; i++) {
			other.neighbours_[i]->sortNeighbours();
		}
	}
	constexpr void deleteBond(Node& other) {
		deleteNeighbour(other);
		other.deleteNeighbour(*this);

		const auto ns = neighboursSize();
		for (AI i = 0; i < ns; i++) {
			neighbours_[i]->sortNeighbours();
		}

		const auto nso = other.neighboursSize();
		for (AI i = 0; i < nso; i++) {
			other.neighbours_[i]->sortNeighbours();
		}
	}
	AI findNeighbour(NeighbourValueType other) const noexcept {
		const auto s = neighboursSize();
		for (AI i = 0; i < s; i++) {
			if (neighbours_[i] == other)
				return i;
		}
		return AI(-1);
	}
	void addNeighboursVector(NeighboursType&& other) {
		neighbours_ = std::move(other);
	}

	// Does not change id or neighbours of neighbours
	constexpr void swap(Node& other) noexcept {
		std::swap(base::type, other.type);
		std::swap(base::hAtoms, other.hAtoms);
		neighbours_.swap(other.neighbours_);

		// Does not change id or neighbours of neighbours
	}
	constexpr void exchangeNeighbour(const NeighbourValueType oldNei, const NeighbourValueType newNei) {
		*(std::find(neighbours_.begin(), neighbours_.end(), oldNei)) = newNei;
	}
	constexpr NeighbourValueType core_neighbours_pop_back() noexcept {
		NeighbourValueType res = std::move(neighbours_.back());
		neighbours_.pop_back();
		return res;
	}
private:
	constexpr void initializeNeighbours() {
		neighbours_.reserve(4);
	}
	constexpr void deleteNeighbour(const Node& node) noexcept {
		size_t s = neighbours_.size();
		size_t i = 0;
		for (; i < s; i++) {
			if (neighbours_[i] == (&node)) break;
		}
		s--;
		for (; i < s; i++) {
			neighbours_[i] = neighbours_[i + 1];
		}
		neighbours_.pop_back();
	}
};

namespace std {
	template<class A, class H, class AI>
	constexpr void swap(Node<A, H, AI>& n1, Node<A, H, AI>& n2) noexcept {
		n1.swap(n2);
	}
}

template<class AI>
// AI is integral type
struct Bond {
public:
	// Data
	AI first = 0;
	AI second = 0;

	// Constructors
	constexpr Bond() noexcept = default;
	constexpr Bond(const AI a1, const AI a2) noexcept : first(a1), second(a2) {};

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
		if (first > second) std::swap(first, second);
	}
	constexpr std::string ToStr() const {
		std::string res("(");
		res += std::to_string(this->first);
		res += ", ";
		res += std::to_string(this->second);
		res += ")";
		// (1, 2)
		return res;
	}
};

template<class AI, class LengthType>
// AI - integral type,
// LengthType - floating point type
struct BondEx : public Bond<AI> {
	LengthType length = 0.0;
	using base = Bond<AI>;

	constexpr BondEx() = default;

	constexpr BondEx(const Bond<AI>& bond, float len) noexcept : Bond<AI>(bond), length(len) {
		Bond<AI>::validate();
	}

	constexpr BondEx(AI a1, AI a2, LengthType l) noexcept
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

	constexpr bool operator==(const BondEx& other) const noexcept {
		return Bond<AI>::operator==(other);
	}
	constexpr auto operator<(const BondEx& other) const noexcept {
		if (!(Bond<AI>::operator==(other)))
			return Bond<AI>::operator<(other);
		return this->length < other.length;
	}

	constexpr std::string ToStr() const {
		std::string res("(");
		res += std::to_string(this->first);
		res += ", ";
		res += std::to_string(this->second);
		res += ", {\"distance\": ";
		res += std::to_string(this->length);
		res += "})";
		return res;
		// (1, 2, {"distance": 1.0})
	}
};
