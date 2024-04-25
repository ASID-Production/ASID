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
#include "Engine.h" // for Node and Bond
#include <sstream> // for explicit string constructor
#include <vector> // for using std::vector
#include <type_traits> // for std::fundamental

template<class A, class H, class AI, class MI>
class MoleculeGraph : private std::vector<Node<A, H, AI> > {
public:
	// Declarations
	using base = std::vector<Node<A, H, AI> >;
	using NodeType = Node<A, H, AI>;
	using BondType = Bond<A>;

private:
	// Data
	MI id_ = 0;

public:
	constexpr MoleculeGraph() noexcept = default;
	constexpr MoleculeGraph(const MoleculeGraph&) = delete;
	constexpr MoleculeGraph(MoleculeGraph&&) noexcept = default;
	explicit constexpr MoleculeGraph(base&& other, MI) noexcept(std::is_nothrow_move_constructible<base>::value)
		: base(std::move(other)) {}
	static constexpr MoleculeGraph<A, H, AI, MI> ReadData(const char* str) {
		MoleculeGraph mg;
		std::stringstream ss(str);
		mg.parseMainstring(ss);
		return mg;
	}
	static constexpr MoleculeGraph<A, H, AI, MI> ReadInput(const char* str) {
		MoleculeGraph mg;
		std::stringstream ss(str);
		const auto sn = mg.parseMainstring(ss);
		mg.parseMultiatom(ss, sn);
		return mg;
	}

	constexpr AI size() const noexcept {
		return static_cast<AI>(base::size());
	}

	// Bond functions
	constexpr std::vector<Bond<AI>> getBonds() const {
		std::vector<Bond<AI>> ret;

		AI s = base::size();
		for (AI i = 1; i < s; i++) {
			const auto neigh_s = this->operator[](i).neighboursSize();
			for (AI j = 0; j < neigh_s; j++) {
				if (j > i) {
					ret.emplace_back(i, j);
				}
			}
		}
		ret.shrink_to_fit();
		return ret;
	}
	constexpr AI countBonds() const {
		AI ret(0);
		for (const auto& node : *this) {
			ret += node.neighboursSize();
		}
		return ret >> 1;
	}

	// Add and delete bonds
	constexpr void addBond(const AI a, const AI b) {
		base::operator[](a).addBondWithSort(base::operator[](b));
	}
	constexpr void deleteBond(const AI a, const AI b) {
		base::operator[](a).deleteBond(base::operator[](b));
	}
	constexpr void addBondsFromVector(const std::vector<BondType>& bond) {
		auto bs = bond.size();
		for (decltype(bs) i = 0; i < bs; i++) {
			addBond(bond[i].first, bond[i].second);
		}
	}

	// For Search
	constexpr AI findStart() const {
		const AI s = size();
		AI m = 1;
		for (AI i = 2; i < s; i++) {
			if (base::operator[](i) > base::operator[](m))
				m = i;
		}
		return m;
	}
	constexpr AI getNeighbourId(AI cur, AI neighbourIt) const noexcept {
		return base::operator[](cur).getNeighbour(neighbourIt)->getID();
	}
	constexpr const NodeType& getNeighbourReference(AI cur, AI neighbourIt) const noexcept {
		return *(base::operator[](cur).getNeighbour(neighbourIt));
	}
	constexpr const NodeType* getNeighbourPointer(AI cur, AI neighbourIt) const noexcept {
		return base::operator[](cur).getNeighbour(neighbourIt);
	}
	constexpr void exchange(AI a1, AI a2) {
		const auto n1size = base::operator[](a1).neighboursSize();
		const auto n2size = base::operator[](a2).neighboursSize();
		auto a1p = &(base::operator[](a1));
		auto a2p = &(base::operator[](a2));
		for (AI i = 0; i < n1size; i++) {
			base::operator[](getNeighbourId(a1, i)).exchangeNeighbour(a1p, a2p);
		}
		for (AI i = 0; i < n2size; i++) {
			base::operator[](getNeighbourId(a2, i)).exchangeNeighbour(a2p, a1p);
		}
		base::operator[](a1).swap(base::operator[](a2));
	}

	// Copy function
	[[nodiscard]] constexpr MoleculeGraph makeCopy() const noexcept(false) {
		// Copy
		MoleculeGraph ret;
		ret.id_ = id_;
		AI s = base::size();
		ret.reserve(s);
		// Convertion to Correct Neighbours
		for (AI i = 0; i < s; i++) {
			auto& node = base::operator[](i);
			ret.emplace_back(node.getType(), node.getHAtoms(), node.getID());
		}
		for (AI i = 0; i < s; i++) {
			auto& node = base::operator[](i);

			using NeighboursType = typename NodeType::NeighboursType;
			NeighboursType ret_nei(0);
			const auto nodeSize = node.neighboursSize();
			for (AI j = 0; j < nodeSize; j++) {
				ret_nei.push_back(&(ret[node.getNeighbour(j)->getID()]));
			}
			ret[i].addNeighboursVector(std::move(ret_nei));
		}
		return ret;
	}

	// Operators
	constexpr const NodeType& operator[](const AI s) const noexcept {
		return base::operator[](s);
	}
	constexpr NodeType& operator[](const AI s) noexcept {
		return base::operator[](s);
	}
	constexpr MoleculeGraph& operator=(const MoleculeGraph& other) noexcept = delete;
	constexpr MoleculeGraph& operator=(MoleculeGraph&& other) noexcept = default;

	// Iterators
	constexpr auto begin() noexcept {
		return base::begin();
	}
	constexpr auto end() noexcept {
		return base::end();
	}
	constexpr auto begin() const noexcept {
		return base::begin();
	}
	constexpr auto end() const noexcept {
		return base::end();
	}
	// Interface Const ID(ref)
	constexpr const auto& getID() const noexcept {
		return id_;
	}
private: 
	inline AI parseMainstring(std::stringstream& ss) {
		MI id {};
		ss >> id_;
		AI sn {};
		ss >> sn;
		AI sb {};
		ss >> sb;
		sn++;
		base::reserve(sn);
		base::emplace_back(A(0), H(0), AI(0));


		// Atomic loop
		for (AI i = 1; i < sn; i++) {
			int a;
			ss >> a;
			int b;
			ss >> b;
			base::emplace_back(A(a), H(b), AI(i));
		}
		// Bond loop
		for (AI i = 0; i < sb; i++) {
			int a;
			ss >> a;
			int b;
			ss >> b;
			base::operator[](a).addBondWithSort(base::operator[](b));
		}
		return sn;
	}
	inline void parseMultiatom(std::stringstream & ss, const AI sn) {
		int xty;

		if (!(ss >> xty))
			return;
		//-1 6 7 8 -2 8 16 -3 9 17 0
		while (xty != 0) {
			A real(static_cast<char>(xty));
			int next_xty = 0;
			auto& xtsize = A::size;
			for (char i = 0; i < xtsize; i++)
			{
				ss >> next_xty;
				if (next_xty <= 0) break;
				real.AddType(static_cast<char>(next_xty));
			}
			for (AI i = 1; i < sn; i++)
			{
				if (this->operator[](i).getType().simple_eq((static_cast<char>(xty)))) {
					this->operator[](i).setType(real);
				}
			}
			// Foolproof
			while (next_xty > 0) {
				ss >> next_xty;
			}
			xty = next_xty;
		}
	}
};
