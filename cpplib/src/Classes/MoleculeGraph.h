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
#include <algorithm> // for std::stable_sort
#include <cctype> // for isdigit
#include <charconv>
#include <list> 
#include <numeric> // for std::iota
#include <optional>
#include <ranges>
#include <sstream>
#include <type_traits> // for std::fundamental
#include <vector> // for using std::vector

#include "../BaseHeaders/BaseTypes.h" // for constants::mend_size
#include "../BaseHeaders/Concepts.h"
#include "../BaseHeaders/Currents.h"
#include "../Classes/Bond.h" // for Bond
#include "../Classes/Engine.h" // for Node

/**
		 * Collects all unique bonds present in the molecule and returns them as a vector of bond pairs.
		 *
		 * @param strategy Allocation hint for the returned vector; when `Exact` reserves the exact number of bonds,
		 *                 when `Estimated` reserves a heuristic amount, and when `None` performs no pre-reservation.
		 * @returns Vector of `BondType` where each element represents a bond (a, b) with `a < b` and both indices are atom IDs.
		 */
		namespace cpplib {		
	enum class ReserveStrategy :char {
		None,
		Exact,
		Estimated
	};

	class TypeMap {
	public:
		using AtomIndex = basic_types::AtomIndex;
		using indexType = int8_t;
	private:
		std::array<AtomIndex, constants::mend_size> data_{};
	public:
		static_assert (INT8_MAX >= constants::mend_size, "mend_size chould be less than INT8_MAX");
		constexpr TypeMap() {
			std::ranges::fill(data_, AtomIndex(-1));
		}
		constexpr explicit TypeMap(const AtomIndex value) {
			for (indexType i = 0; i < constants::mend_size; i++)
			{
				data_[i] = value;
			}
		}
		inline AtomIndex& operator[](const indexType i) {
			assert(i < data_.size());
			return data_[i];
		}
		inline const AtomIndex& operator[](const indexType i) const {
			assert(i < data_.size());
			return data_[i];
		}
		inline void initialize(const basic_types::TypeBitset& bits) {
			for (indexType i = 1; i < constants::mend_size; i++)
			{
				if (bits[i]) data_[i] = AtomIndex(0);
			}
		}
		constexpr indexType size() const {
			return data_.size();
		}
		inline bool isFinished() const {
			for (indexType i = 1; i < constants::mend_size; i++)
			{
				if (data_[i] > 0) return false;
			}
			return true;
		}
	};

	template<AtomTypeConcept A>
	class MoleculeParser;

	template<AtomTypeConcept A>
	class MoleculeCore {
	public:
		// Definitions
		using NodeType = Node<A>;
		using NodeContainer = ::std::vector<NodeType>;
		using BondType = Bond;
		using AtomIndex = basic_types::AtomIndex;
		using MoleculeIndex = basic_types::MoleculeIndex;
		using HType = typename NodeType::HType;
		using AtomTypeBase = typename A::AtomTypeBase;

		template<AtomTypeConcept T>
		friend class MoleculeCore;

		friend class MoleculeParser<A>;

	private:
		// Data
		NodeContainer data_;
		MoleculeIndex id_ = 0;

	public:
		constexpr MoleculeCore() noexcept = default;
		explicit constexpr MoleculeCore(NodeContainer&& other) noexcept(::std::is_nothrow_move_constructible_v<NodeContainer>)
			: data_(::std::move(other)) {
		}

		constexpr AtomIndex size() const noexcept {
			return static_cast<AtomIndex>(data_.size());
		}

		// Bond functions
		::std::vector<BondType> getBonds(ReserveStrategy strategy = ReserveStrategy::Estimated) const {
			::std::vector<BondType> ret;
			switch (strategy) {
			case ReserveStrategy::Exact:
				ret.reserve(countBonds());
				break;
			case ReserveStrategy::Estimated:
				ret.reserve(size() * 2);
				break;
			default: break;
			}


			AtomIndex s = size();
			for (AtomIndex i = 1; i < s; i++) {
				const auto& node = data_[i];
				const auto neigh_s = node.neighboursSize();

				for (AtomIndex j = 0; j < neigh_s; j++) {
					const AtomIndex neighbor_id = node.getNeighbour(j)->getID();
					if (neighbor_id > i) {
						ret.emplace_back(i, neighbor_id);
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

		// Copy functions
		template <class OT>
		[[nodiscard]] constexpr MoleculeCore<OT> makeCopyEx() const noexcept(false) {
			// Copy
			MoleculeCore<OT> ret;
			ret.id_ = id_;
			AtomIndex s = size();
			ret.data_.reserve(s);
			// Convertion to Correct Neighbours
			for (AtomIndex i = 0; i < s; i++) {
				const auto& node = data_[i];
				ret.data_.emplace_back(static_cast<OT>(static_cast<AtomTypeBase>(node.getType())), node.getHAtoms(), node.getID());
				ret.data_.back().addNeighboursVector(node.getNeighboursVector());
				ret.data_.back().setCoord(std::move(node.getCoord()));
			}
			return ret;
		}

		// Operators
		constexpr const NodeType& operator[](const AtomIndex s) const noexcept {
			assert(s < data_.size());
			return data_[s];
		}
		constexpr NodeType& operator[](const AtomIndex s) noexcept {
			assert(s < data_.size());
			return data_[s];
		}

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

		constexpr TypeMap getTypeMap() const {
			TypeMap map;
			static_assert(map.size() == constants::mend_size);
			AtomIndex s = size();
			for (AtomIndex i = 1; i < s; i++)
			{
				const auto t = static_cast<basic_types::AtomTypeBase>(data_[i].getType());
				const auto h = data_[i].getHAtoms();

				if (map[1] == AtomIndex(-1))
					map[1] = AtomIndex(h);
				else
					map[1] += h;
				if (t > 0) {
					if (map[t] == AtomIndex(-1))
						map[t] = AtomIndex(1);
					else
						map[t]++;
				}
				else {
					for (TypeMap::indexType j = 1; j < map.size(); j++)
					{
						if (data_[i].getType().contains(j) && map[j] == AtomIndex(-1))
							map[j] = AtomIndex(0);
					}
				}
			}
			return map;
		}
		void unpackHydrogens(AtomIndex index) {
			auto s = data_[index].getHAtoms();
			for (AtomIndex i = 0; i < s; i++)
			{
				AtomIndex last = data_.size();
				data_.emplace_back(A(1), 0, last);
				addBond(index, last);
				data_[last].setCoord(Coord(1, constants::maxNeighbours));
			}
			data_[index].setHAtoms(0);
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
					data_[i].swap(data_[best]);
					// IDs don't swap and stay actual
				}
			}
			for (AtomIndex i = 1; i < s; i++)
			{
				data_[i].sortNeighbours();
			}
		}

	};

	template<AtomTypeConcept A>
	class MoleculeParser {
	public:
		using GraphType = MoleculeCore<A>;
		using AtomIndex = typename GraphType::AtomIndex;
		using HType = typename GraphType::HType;


		explicit MoleculeParser(GraphType& graph) noexcept : graph_(graph) {}

		static ::std::string ResortString(const char* str) {
			GraphType gr;
			MoleculeParser p(gr);
			p.parseMainstringData(str, TypeMap(0));
			gr.sortGraph();

			return p.writeDataString();
		}

		static ::std::string _ParseOldInputString(const char* str) {
			::std::string ret(str);
			long na;
			long no;
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

		template<typename T = A> requires std::is_same_v<T, currents::AtomTypeRequest>
		static ::std::pair<GraphType, basic_types::TypeBitset> Read(const char* str) {
			GraphType mg;
			MoleculeParser p(mg);

			const auto sn = p.parseMainstringRequest(str);
			auto multiAtomBits = p.parseMultiatom(str, sn);
			p.release_HAtoms(multiAtomBits);
			mg.sortGraph();
			return { mg, multiAtomBits };
		}
		template<typename T = A> requires std::is_same_v<T, currents::AtomTypeData>
		static ::std::pair<GraphType, bool> Read(const char* str, const basic_types::TypeBitset& multiAtomBits, const TypeMap& map) {
			GraphType mg;
			MoleculeParser p(mg);

			const auto is_correct = p.parseMainstringData(str, map);
			if (is_correct == false) return { GraphType(), false };
			p.release_HAtoms(multiAtomBits);
			return { mg, true };
		}

	private:
		constexpr bool readToNext(const char*& str) const noexcept {
			while (*str != '-' && *str != '\0' && (*str > '9' || *str < '0')) {
				str++;
			}
			return *str != '\0';
		}
		AtomIndex readSingleInt(const char*& str) const {
			readToNext(str);
			AtomIndex value;
			auto [ptr, ec] = std::from_chars(str, str + 15, value);
			if (ec == std::errc()) str = ptr;
			else throw std::invalid_argument("Invalid number");
			return value;
		}
		::std::vector<bool> parseAtomsBlockData(const char*& str, const AtomIndex sn, const TypeMap& argMap) {
			graph_.data_.reserve(sn);
			::std::vector<bool> is_used(sn, false);
			TypeMap map(argMap);
			graph_.data_.emplace_back(A(0), HType(0), AtomIndex(0));

			for (AtomIndex i = 1; i < sn; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				if (a > 0) {
					if (map[a] == -1)
						continue;
					if (map[a] > 0)
						map[a]--;
				}
				if (b > 0) {
					if (map[1] > 0) {
						map[1] -= b;
						if (map[1] < 0)
							map[1] = 0;
					}
				}
				is_used[i] = true;
				graph_.data_.emplace_back(A(a), HType(b), AtomIndex(graph_.data_.size()));
			}
			if (!map.isFinished())
				is_used[0] = false;
			else
				is_used[0] = true;
			return is_used;
		}
		void parseAtomsBlockRequest(const char*& str, const AtomIndex sn) {
			graph_.data_.reserve(sn);
			graph_.data_.emplace_back(A(0), HType(0), AtomIndex(0));
			for (AtomIndex i = 1; i < sn; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				graph_.data_.emplace_back(A(a), HType(b), AtomIndex(i));
				int first = readSingleInt(str);
				int second = readSingleInt(str);
				graph_.data_.back().setCoord(Coord(static_cast<Coord::argumentType>(first), static_cast<Coord::argumentType>(second)));
			}
		}
		::std::pair<AtomIndex, AtomIndex> parseInit(const char*& str) {
			graph_.id_ = readSingleInt(str);
			::std::pair<AtomIndex, AtomIndex> r;
			r.first = readSingleInt(str);
			r.second = readSingleInt(str);
			r.first++;
			return r;
		}

		template<typename T = A> requires std::is_same_v<T, currents::AtomTypeData>
		bool parseMainstringData(const char*& str, const TypeMap& map) {

			::std::pair<AtomIndex, AtomIndex>&& sn_sb = parseInit(str);
			AtomIndex& sn = sn_sb.first;
			AtomIndex& sb = sn_sb.second;

			// Atomic loop
			auto used = parseAtomsBlockData(str, sn, map);
			if (used[0] == false) return false;
			std::vector<AtomIndex> reI(sn, 0);
			AtomIndex reI_last = 1;
			for (AtomIndex i = 1; i < sn; i++)
			{
				if (used[i]) {
					reI[i] = reI_last;
					reI_last++;
				}
			}

			std::vector<Coord::innerType> coord_counters(sn, 0);

			// Bond loop
			for (AtomIndex i = 0; i < sb; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				bool b_reIa = reI[a] != 0;
				bool b_reIb = reI[b] != 0;
				if (b_reIa && b_reIb)
					graph_[reI[a]].addBondSimple(graph_[reI[b]]);
				else {
					if (b_reIa)
						coord_counters[reI[a]]++;
					if (b_reIb)
						coord_counters[reI[b]]++;
				}
			}
			for (AtomIndex i = 1; i < reI_last; i++) {
				auto c = Coord(static_cast<Coord::argumentType>(graph_[i].getHAtoms() + graph_[i].neighboursSize() + coord_counters[i]));
				graph_[i].setCoord(::std::move(c));
			}
			for (AtomIndex i = 0; i < reI_last; i++) {
				graph_[i].sortNeighbours();
			}

			return true;
		}

		template<typename T = A> requires std::is_same_v<T, currents::AtomTypeRequest>
		AtomIndex parseMainstringRequest(const char*& str) {
			::std::pair<AtomIndex, AtomIndex>&& sn_sb = parseInit(str);
			AtomIndex& sn = sn_sb.first;
			AtomIndex& sb = sn_sb.second;

			// Atomic loop
			parseAtomsBlockRequest(str, sn);

			// Bond loop
			for (AtomIndex i = 0; i < sb; i++) {
				int a = readSingleInt(str);
				int b = readSingleInt(str);
				graph_.data_[a].addBondSimple(graph_.data_[b]);
			}
			for (AtomIndex i = 0; i < sn; i++) {
				graph_.data_[i].sortNeighbours();
			}

			return sn;
		}

		template<typename T = A> requires std::is_same_v<T, currents::AtomTypeRequest>
		basic_types::TypeBitset parseMultiatom(const char* str, const AtomIndex sn) {
			AtomIndex xty;
			basic_types::TypeBitset multiAtomBits;
			if (readToNext(str))
				xty = readSingleInt(str);
			else {
				return multiAtomBits;
			}
			//-1 6 7 8 -2 8 16 -3 9 17 0
			while (xty != 0) {
				A real(static_cast<char>(xty));
				int next_xty = 0;
				for (char i = 0; i < constants::mend_size; i++)
				{
					next_xty = readSingleInt(str);
					if (next_xty <= 0) break;
					real.AddType(static_cast<char>(next_xty));
				}
				for (AtomIndex i = 1; i < sn; i++)
				{
					if (static_cast<basic_types::AtomTypeBase>(graph_.data_[i].getType()) == static_cast<basic_types::AtomTypeBase>(xty)) {
						graph_.data_[i].setType(real);
						if (!real.contains(1)) continue;

						for (AtomIndex j = 0; j < graph_.operator[](i).neighboursSize(); j++)
						{
							multiAtomBits |= graph_.operator[](i).getNeighbour(j)->getType().get_bitset();
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

		void release_HAtoms(const basic_types::TypeBitset& bits) {
			auto& data_ = graph_.data_;
			if (bits.none()) return;
			const AtomIndex original_size = graph_.size();
			AtomIndex new_index = original_size;

			size_t totalH = 0;
			for (AtomIndex i = 1; i < original_size; i++) {
				totalH += data_[i].getHAtoms();
			}
			data_.reserve(original_size + totalH);

			for (AtomIndex i = 1; i < original_size; i++) {
				bool check = false;

				if constexpr (::std::is_same_v<A, CompositeAtom>) {
					check = (data_[i].getType().get_bitset() & bits).any();
				}
				else {
					const auto type = static_cast<basic_types::AtomTypeBase>(data_[i].getType());
					check = (type > 0) && bits[type];
				}

				if (check) {
					const auto hAtoms = data_[i].getHAtoms();
					for (HType j = 0; j < hAtoms; j++) {
						data_.emplace_back(A(1), 0, new_index);
						data_[i].addBondWithSort(data_.back());
						data_.back().setCoord(Coord(1, constants::maxNeighbours));
						new_index++;
					}
					data_[i].setHAtoms(0);
				}
			}
		}
		::std::string writeDataString() const {
			using BondType = typename MoleculeCore<A>::BondType;
			AtomIndex ns = graph_.size();
			AtomIndex bs = 0;
			::std::vector<BondType> bonds;
			::std::string bond_str; // starts with ' '
			::std::string node_str;

			bond_str.reserve(2048);
			node_str.reserve(2048);

			for (AtomIndex i = 1; i < ns; i++)
			{
				node_str += ' ';
				node_str += std::to_string(static_cast<int>(static_cast<A::AtomTypeBase>(graph_[i].getType())));
				node_str += ' ';
				node_str += std::to_string(static_cast<int>(graph_[i].getHAtoms()));

				int8_t neis = graph_[i].neighboursSize();
				for (int8_t j = 0; j < neis; j++)
				{
					AtomIndex neindex = graph_[i].getNeighbour(j)->getID();
					if (neindex > i) {
						bond_str += ' ';
						bond_str += std::to_string(i);
						bond_str += ' ';
						bond_str += std::to_string(neindex);
					}
				}
				bs += neis;
			}
			bs >>= 1;

			std::string res = (std::to_string(graph_.id_) + ' ')
				+ (std::to_string(ns - 1) + ' ') + std::to_string(bs);
			return res + node_str + bond_str;
		}

	private:
		// Data ref
		GraphType& graph_;
	};
}