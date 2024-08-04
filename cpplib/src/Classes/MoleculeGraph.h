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
#include "../BaseHeaders/Support.h" // for mend_size
#include <sstream> // for explicit string constructor
#include <vector> // for using std::vector
#include <list> // for using std::list
#include <type_traits> // for std::fundamental
#include <numeric> // for std::iota
#include <algorithm> // for std::stable_sort
namespace cpplib {
	template<class A> class MoleculeGraph  {
	public:
		// Declarations
		using NodeType = Node<A>;
		using NodeContainer = ::std::vector<NodeType>;
		using BondType = currents::BondType;
		using AtomIndex = currents::AtomIndex;
		using MoleculeIndex = currents::MoleculeIndex;
		using HType = typename NodeType::HType;

	private:
		// Data
		NodeContainer data_;
		MoleculeIndex id_ = 0;

	public:
		constexpr MoleculeGraph() noexcept = default;
		constexpr MoleculeGraph(const MoleculeGraph&) = delete;
		constexpr MoleculeGraph(MoleculeGraph&&) noexcept = default;
		explicit constexpr MoleculeGraph(NodeContainer&& other) noexcept(::std::is_nothrow_move_constructible<NodeContainer>::value)
			: data_(::std::move(other)) {}
		
		static constexpr MoleculeGraph ReadData(const char* str, const ::std::bitset<mend_size>& multiAtomBits) {
			MoleculeGraph mg;
			::std::stringstream ss(str);
			const auto sn = mg.parseMainstring<false>(ss);

			mg.release_HAtoms_CharAtom(multiAtomBits, sn);
			//mg.sortGraph();
			return mg;
		}

		static constexpr ::std::pair<MoleculeGraph, ::std::bitset<mend_size>> ReadInput(const char* str) {
			MoleculeGraph mg;
			::std::stringstream ss(str);
			const auto sn = mg.parseMainstring<true>(ss);
			auto multiAtomBits = mg.parseMultiatom(ss, sn);
			mg.release_HAtoms_XAtom(multiAtomBits, sn);
			//mg.sortGraph();
			return ::std::make_pair(std::move(mg), ::std::move(multiAtomBits));
		}


		constexpr AtomIndex size() const noexcept {
			return static_cast<AtomIndex>(data_.size());
		}

		// Bond functions
		constexpr ::std::vector<BondType> getBonds() const {
			std::vector<BondType> ret;

			AtomIndex s = size();
			for (AtomIndex i = 1; i < s; i++) {
				const auto neigh_s = this->operator[](i).neighboursSize();
				for (AtomIndex j = 0; j < neigh_s; j++) {
					if (j > i) {
						ret.emplace_back(i, j);
					}
				}
			}
			ret.shrink_to_fit();
			return ret;
		}
		constexpr AtomIndex countBonds() const {
			AtomIndex ret(0);
			for (const auto& node : *this) {
				ret += node.neighboursSize();
			}
			return ret >> 1;
		}

		// Add and delete bonds
		constexpr void addBond(const AtomIndex a, const AtomIndex b) {
			data_[a].addBondWithSort(data_[b]);
		}
		constexpr void deleteBond(const AtomIndex a, const AtomIndex b) {
			data_[a].deleteBond(data_[b]);
		}
		constexpr void addBondsFromVector(const ::std::vector<BondType>& bond) {
			auto bs = bond.size();
			for (decltype(bs) i = 0; i < bs; i++) {
				addBond(bond[i].first, bond[i].second);
			}
		}

		// For Search
		constexpr AtomIndex findStart() const {
			const AtomIndex s = size();
			AtomIndex m = 1;
			for (AtomIndex i = 2; i < s; i++) {
				if ((data_[i]).RawMore(data_[m]))
					m = i;
			}
			return m;
		}
		constexpr AtomIndex getNeighbourId(AtomIndex cur, AtomIndex neighbourIt) const noexcept {
			return data_[cur].getNeighbour(neighbourIt)->getID();
		}
		constexpr const NodeType& getNeighbourReference(AtomIndex cur, AtomIndex neighbourIt) const noexcept {
			return *(data_[cur].getNeighbour(neighbourIt));
		}
		constexpr const NodeType* getNeighbourPointer(AtomIndex cur, AtomIndex neighbourIt) const noexcept {
			return data_[cur].getNeighbour(neighbourIt);
		}
		constexpr void exchange(AtomIndex a1, AtomIndex a2) {
			const auto n1size = data_[a1].neighboursSize();
			const auto n2size = data_[a2].neighboursSize();
			auto a1p = &(data_[a1]);
			auto a2p = &(data_[a2]);
			for (AtomIndex i = 0; i < n1size; i++) {
				data_[getNeighbourId(a1, i)].exchangeNeighbour(a1p, a2p);
			}
			for (AtomIndex i = 0; i < n2size; i++) {
				data_[getNeighbourId(a2, i)].exchangeNeighbour(a2p, a1p);
			}
			data_[a1].swap(data_[a2]);
		}

		// Copy function
		[[nodiscard]] constexpr MoleculeGraph makeCopy() const noexcept(false) {
			// Copy
			MoleculeGraph ret;
			ret.id_ = id_;
			AtomIndex s = size();
			ret.data_.reserve(s);
			// Convertion to Correct Neighbours
			for (AtomIndex i = 0; i < s; i++) {
				const auto& node = data_[i];
				ret.data_.emplace_back(node);
			}
			for (AtomIndex i = 0; i < s; i++) {
				auto& node = data_[i];

				using NeighboursType = typename NodeType::NeighboursType;
				NeighboursType ret_nei(0);
				const auto nodeSize = node.neighboursSize();
				for (AtomIndex j = 0; j < nodeSize; j++) {
					ret_nei.push_back(&(ret[node.getNeighbour(j)->getID()]));
				}
				ret[i].addNeighboursVector(::std::move(ret_nei));
			}
			return ret;
		}

		// Operators
		constexpr const NodeType& operator[](const AtomIndex s) const noexcept {
			return data_[s];
		}
		constexpr NodeType& operator[](const AtomIndex s) noexcept {
			return data_[s];
		}
		inline MoleculeGraph& operator=(const MoleculeGraph& other) noexcept = delete;
		inline MoleculeGraph& operator=(MoleculeGraph&& other) noexcept = default;

		// Iterators
		constexpr auto begin() noexcept {
			return data_.begin();
		}
		constexpr auto end() noexcept {
			return data_.end();
		}
		constexpr auto begin() const noexcept {
			return data_.begin();
		}
		constexpr auto end() const noexcept {
			return data_.end();
		}
		// Interface Const ID(ref)
		constexpr const auto& getID() const noexcept {
			return id_;
		}

		// static section for old tests

		static ::std::string _ParseOldInputString(const char* str)  {
			::std::string ret(str);
			int na;
			int no;
			{
				::std::stringstream ss(ret);
				ss >> no;
				ss >> na;
			}
			auto pos = ret.find_first_of(' ', 0);
			pos = ret.find_first_of(' ', pos + 1);
			pos = ret.find_first_of(' ', pos + 1);
			for (int i = 0; i < na; i++)
			{
				pos = ret.find_first_of(' ', pos + 1);
				pos = ret.find_first_of(' ', pos + 1);
				if (ret.length() > pos)
					ret.insert(pos, " 0 14");
				else
					ret.append(" 0 14");
				pos += 5;
			}
			return ret;
		}

	private:
		template<bool is_request> inline void parseAtomsBlock(::std::stringstream& ss, const AtomIndex sn) {
			data_.reserve(sn);
			data_.emplace_back(A(0), HType(0), AtomIndex(0));
			for (AtomIndex i = 1; i < sn; i++) {
				int a;
				ss >> a;
				int b;
				ss >> b;
				data_.emplace_back(A(a), HType(b), AtomIndex(i));
				if (is_request) {
					int first;
					int second;
					ss >> first;
					ss >> second;
					data_.back().setCoord(Coord(static_cast<Coord::argumentType>(first), static_cast<Coord::argumentType>(second)));
				}
			}
		}
		::std::pair<AtomIndex, AtomIndex> parseInit(std::stringstream& ss) {
			ss >> id_;
			::std::pair<AtomIndex, AtomIndex> r;
			ss >> r.first;
			ss >> r.second;
			r.first++;
			return r;
		}
		template <bool is_request> inline AtomIndex parseMainstring(::std::stringstream& ss) {
			::std::pair<AtomIndex, AtomIndex>&& sn_sb = parseInit(ss);
			AtomIndex& sn = sn_sb.first;
			AtomIndex& sb = sn_sb.second;

			// Atomic loop
			parseAtomsBlock<is_request>(ss, sn);

			// Bond loop
			for (AtomIndex i = 0; i < sb; i++) {
				int a;
				ss >> a;
				int b;
				ss >> b;
				data_[a].addBondWithSort(data_[b]);
			}
			if (is_request == false) {
				for (AtomIndex i = 1; i < sn; i++) {
					Coord c = Coord(static_cast<Coord::argumentType>(data_[i].getHAtoms() + data_[i].neighboursSize()));
					data_[i].setCoord(::std::move(c));
				}
			}
			return sn;
		}
		::std::bitset<mend_size> parseMultiatom(::std::stringstream& ss, const AtomIndex sn) {
			int xty;
			::std::bitset<mend_size> multiAtomBits;

			if (!(ss >> xty))
				return multiAtomBits;
			//-1 6 7 8 -2 8 16 -3 9 17 0
			while (xty != 0) {
				A real(static_cast<char>(xty));
				int next_xty = 0;
				for (char i = 0; i < mend_size; i++)
				{
					ss >> next_xty;
					if (next_xty <= 0) break;
					real.AddType(static_cast<char>(next_xty));
				}
				for (AtomIndex i = 1; i < sn; i++)
				{
					if (this->operator[](i).getType().simple_eq((static_cast<char>(xty)))) {
						this->operator[](i).setType(real);
						if (!real.include(1)) continue;

						for (AtomIndex j = 0; j < this->operator[](i).neighboursSize(); j++)
						{
							multiAtomBits |= this->operator[](i).getNeighbour(j)->getType().get_bitset();
						}
					}
				}
				// Foolproof
				while (next_xty > 0) {
					ss >> next_xty;
				}
				xty = next_xty;
			}

			return multiAtomBits;
		}
		void release_HAtoms_XAtom(const ::std::bitset<mend_size>& bits, const AtomIndex sn) {
			if (bits.none()) return;
			size_t hs = this->size();
			::std::list< NodeType> hydrogenAtoms;

			for (AtomIndex i = 1; i < sn; i++)
			{
				if ((((XAtom)this->operator[](i).getType()).get_bitset() & bits).any()) {

					auto hAtoms = this->operator[](i).getHAtoms();
					for (HType j = 0; j < hAtoms; j++)
					{
						hydrogenAtoms.emplace_back(XAtom(1), 0, hs);
						this->operator[](i).addBondWithSort(hydrogenAtoms.back());
						hs++;
					}
					this->operator[](i).setHAtoms(0);
				}
			}
			AddNewAtoms(std::move(hydrogenAtoms), hs);

		}
		void release_HAtoms_CharAtom(const std::bitset<mend_size>& bits, const AtomIndex sn) {
			if (bits.none()) return;
			size_t hs = this->size();
			std::list<NodeType> hydrogenAtoms;

			for (AtomIndex i = 1; i < sn; i++)
			{
				if (bits[(char)(this->operator[](i).getType())]) {
					auto hAtoms = this->operator[](i).getHAtoms();
					for (HType j = 0; j < hAtoms; j++)
					{
						hydrogenAtoms.emplace_back(1, 0, static_cast<AtomIndex>(hs));
						this->operator[](i).addBondWithSort(hydrogenAtoms.back());
						hs++;
					}
					this->operator[](i).setHAtoms(0);
				}
			}
			AddNewAtoms(std::move(hydrogenAtoms), hs);
		}


		void AddNewAtoms(std::list<NodeType>&& l, AtomIndex size_plus_l) {
			std::vector<NodeType> ex(size_plus_l);

			const AtomIndex s = size();

			for (AtomIndex i = 1; i < s; i++)
			{
				ex[i].swap(data_[i]);
				ex[i].setID(i);
			}

			auto iter = l.begin();
			for (AtomIndex i = s; i < size_plus_l; i++) {

				ex[i].swap(*iter);
				ex[i].setID(i);
				iter++;
			}
			std::swap(ex, data_);
		}

		void sortGraph() {
			const auto s = this->size();
			std::vector<AtomIndex> order(s);
			std::iota(order.begin(), order.end(), 0);
			auto MoreLambda = [this](AtomIndex a, AtomIndex b) {return this->operator[](a).RawMore(this->operator[](b)); };
			std::stable_sort((++order.begin()), order.end(), MoreLambda);

			std::vector<AtomIndex> reorder(s,0);
			for (typename ::std::remove_const<decltype(s)>::type i = 0; i < s; i++)
			{
				reorder[order[i]] = i;
			}

			for (typename ::std::remove_const<decltype(s)>::type i = 1; i < s; i++)
			{
				if (reorder[i] == i)
					continue;
				decltype(i) other = order[i];
				order[reorder[i]] = other;
				
				exchange(i, other);
			}
			return;
		}
	};
}
