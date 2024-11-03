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
#include <vector> // for using std::vector
#include <list> // for using std::list
#include <type_traits> // for std::fundamental
#include <numeric> // for std::iota
#include <algorithm> // for std::stable_sort
#include <sstream>
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

		template <class OT>
		friend class MoleculeGraph;

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
		
		static MoleculeGraph ReadData(const char* str, const ::std::bitset<mend_size>& multiAtomBits) {
			MoleculeGraph mg;
			//::std::stringstream ss(str);
			const auto sn = mg.parseMainstring<false>(str);
			
			mg.release_HAtoms(multiAtomBits, sn);
			//mg.sortGraph();
			return mg;
		}

		static ::std::pair<MoleculeGraph, ::std::bitset<mend_size>> ReadInput(const char* str) {
			MoleculeGraph mg;
			const auto sn = mg.parseMainstring<true>(str);
			auto multiAtomBits = mg.parseMultiatom(str, sn);
			mg.release_HAtoms(multiAtomBits, sn);
			mg.sortGraph();
			return ::std::make_pair(std::move(mg), ::std::move(multiAtomBits));
		}


		constexpr AtomIndex size() const noexcept {
			return static_cast<AtomIndex>(data_.size());
		}

		// Bond functions
		::std::vector<BondType> getBonds() const {
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
			return 1;
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

		// Copy functions
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
			return ret;
		}
		template <class OT>
		[[nodiscard]] constexpr MoleculeGraph<OT> makeCopyEx() const noexcept(false) {
			// Copy
			MoleculeGraph<OT> ret;
			ret.id_ = id_;
			AtomIndex s = size();
			ret.data_.reserve(s);
			// Convertion to Correct Neighbours
			for (AtomIndex i = 0; i < s; i++) {
				const auto& node = data_[i];
				ret.data_.emplace_back(XAtom(node.getType()).get_simple(), node.getHAtoms(), node.getID());
				ret.data_.back().addNeighboursVector(node.getNeighboursVector());
				ret.data_.back().setCoord(std::move(node.getCoord()));
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

		void sortGraph() {
			// Seclection sort

			AtomIndex s = data_.size();
			AtomIndex best = 1;
			for (AtomIndex i = 1; i < s; i++)
			{
				best = i;
				for (AtomIndex j = i + 1; j < s; j++)
				{
					if (data_[best] < data_[j]) {
						best = j;
					}
				}
				if (best != i) {
					exchange(i, best);
				}
			}
			for (AtomIndex i = 1; i < s; i++)
			{
				data_[i].sortNeighboursSimple();
			}
		}

		static std::string ResortString(const char* str) {
			MoleculeGraph mg;
			mg.parseMainstring<false>(str);
			mg.sortGraph();
			
			return mg.writeDataString();
		}

	private:
		bool readToNext(const char*& str) {

			while (*str != '-' && *str != '\0' && (*str > '9' || *str < '0')) {
				str++;
			}
			return *str != '\0';
		}
		AtomIndex readSingleInt(const char*& str) {
			readToNext(str);
			bool minus = false;
			if (*str == '-') {
				minus = true;
				str++;
			}
			AtomIndex res = (*str)-'0';
			str++;
			while ((*str)!=' ' && (*str)!='\0')
			{
				res *= 10;
				res += (*str) - '0';
				str++;
			}
			if (minus)
				return -res;
			else
				return res;
		}
		template<bool is_request> inline void parseAtomsBlock(const char *& str, const AtomIndex sn) {
			data_.reserve(sn);
			data_.emplace_back(A(0), HType(0), AtomIndex(0));
			for (AtomIndex i = 1; i < sn; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				data_.emplace_back(A(a), HType(b), AtomIndex(i));
				if (is_request) {
					int first = readSingleInt(str);
					int second = readSingleInt(str);
					data_.back().setCoord(Coord(static_cast<Coord::argumentType>(first), static_cast<Coord::argumentType>(second)));
				}
			}
		}
		::std::pair<AtomIndex, AtomIndex> parseInit(const char*& str) {
			id_ = readSingleInt(str);
			::std::pair<AtomIndex, AtomIndex> r;
			r.first = readSingleInt(str);
			r.second = readSingleInt(str);
			r.first++;
			return r;
		}
		template <bool is_request> inline AtomIndex parseMainstring(const char*& str) {
			::std::pair<AtomIndex, AtomIndex>&& sn_sb = parseInit(str);
			AtomIndex& sn = sn_sb.first;
			AtomIndex& sb = sn_sb.second;

			// Atomic loop
			parseAtomsBlock<is_request>(str, sn);

			// Bond loop
			for (AtomIndex i = 0; i < sb; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				data_[a].addBondSimple(data_[b]);

				//data_[a].addBondWithSort(data_[b]);
			}
			if (is_request == false) {
				for (AtomIndex i = 1; i < sn; i++) {
					Coord c = Coord(static_cast<Coord::argumentType>(data_[i].getHAtoms() + data_[i].neighboursSize()));
					data_[i].setCoord(::std::move(c));
				}
			}
			for (AtomIndex i = 0; i < sn; i++) {
				data_[i].sortNeighbours();
			}

			return sn;
		}
		::std::bitset<mend_size> parseMultiatom(const char* str, const AtomIndex sn) {
			AtomIndex xty;
			::std::bitset<mend_size> multiAtomBits;
			if(readToNext(str))
				xty = readSingleInt(str);
			else {
				return multiAtomBits;
			}
			//-1 6 7 8 -2 8 16 -3 9 17 0
			while (xty != 0) {
				A real(static_cast<char>(xty));
				int next_xty = 0;
				for (char i = 0; i < mend_size; i++)
				{
					next_xty = readSingleInt(str);
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
					next_xty = readSingleInt(str);
				}
				xty = next_xty;
			}

			return multiAtomBits;
		}
		void release_HAtoms(const ::std::bitset<mend_size>& bits, const AtomIndex sn) {
			if (bits.none()) return;
			size_t hs = this->size();
			//::std::list< NodeType> hydrogenAtoms;

			for (AtomIndex i = 1; i < sn; i++)
			{
				bool check;
				if (::std::is_same<A, XAtom>::value == true) {
					check = ((static_cast<XAtom>(data_[i].getType())).get_bitset() & bits).any();
				}
				else {
					check = bits[(currents::AtomTypeData)(this->operator[](i).getType())];
				}
				if(check) {
					auto hAtoms = this->operator[](i).getHAtoms();
					for (HType j = 0; j < hAtoms; j++)
					{
						data_.emplace_back(A(1), 0, hs);
						data_[i].addBondWithSort(data_.back());
						hs++;
					}
					data_[i].setHAtoms(0);
				}
			}
		}
		constexpr void exchange(AtomIndex a1, AtomIndex a2) {
			data_[a1].swap(data_[a2]);
		}
		std::string writeDataString() const {
			AtomIndex ns = data_.size();
			AtomIndex bs = 0;
			std::vector<BondType> bonds;
			std::string bond_str; // starts with ' '
			std::string node_str;
			for (AtomIndex i = 1; i < ns; i++)
			{
				node_str += ' ';
				node_str += std::to_string(static_cast<int>(static_cast<char>(data_[i].getType())));
				node_str += ' ';
				node_str += std::to_string(static_cast<int>(data_[i].getHAtoms()));

				int8_t neis = data_[i].neighboursSize();
				for (int8_t j = 0; j < neis; j++)
				{
					AtomIndex neindex = getNeighbourId(i, j);
					if(neindex > i) {
						bond_str += ' ';
						bond_str += std::to_string(i);
						bond_str += ' ';
						bond_str += std::to_string(neindex);
					}					
				}
				bs += neis;
			}
			bs >>= 1;

			std::string res = (std::to_string(id_) + ' ')
				+ (std::to_string(ns-1) + ' ') + std::to_string(bs);
			return res + node_str + bond_str;
		}
	};
}
