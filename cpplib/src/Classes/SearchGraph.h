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
#include "MoleculeGraph.h"
#include <list>
#include <vector>
namespace cpplib {
	class SearchGraph {
	public:
		// Declarations
		using AtomIndex = currents::AtomIndex;
		using MoleculeIndex = currents::MoleculeIndex;
		using RequestGraphType = MoleculeGraph<currents::AtomTypeRequest>;
		using DatabaseGraphType = MoleculeGraph<currents::AtomTypeData>;

		using BondType = DatabaseGraphType::BondType;
		using RequestNodeType = RequestGraphType::NodeType;
		using DatabaseNodeType = DatabaseGraphType::NodeType;
		using CompareVectorType = ::std::vector<AtomIndex>;
		using Log = ::std::list<::std::pair<BondType, BondType>>;

	private:
		// Sizes
		AtomIndex inputSize_ = 0;
		AtomIndex dataSize_ = 0;
		// Data
		RequestGraphType input_;
		DatabaseGraphType data_;
		CompareVectorType comp_;
		Log log_;
		::std::vector<bool> usedInComp_;

	public:
		SearchGraph() {}
		inline void setupInput(RequestGraphType&& molGraph) noexcept {
			input_ = ::std::move(molGraph);
			inputSize_ = input_.size();
		}
		inline void setupData(DatabaseGraphType&& molGraph) noexcept {
			data_ = ::std::move(molGraph);
			//data_.sortGraph();
			dataSize_ = data_.size();
		}

		// Input and Data should be ready
		inline void prepareToSearch() {
			this->prepareHAtoms();
			comp_.assign(inputSize_, 0);
			log_.clear();
			usedInComp_.assign(dataSize_, false);
		}
		bool searchTry(AtomIndex startI, AtomIndex startD, bool exact) {
			addComp(startI, startD);
			bool result;
				if (input_[startI].hasNeighbours()) {
					result = recursiveSearchHasNeighbours(startI, exact);
				}
				else {
					result = recursiveSearchNoNeighbours(exact);
				}
			
			if(result)
				return true;
			else
				deleteComp(startI);
			return false;
		}
		// destroys all data, need reinitialization!
		bool startFullSearch(const bool exact, AtomIndex startAtom = 0) {
			if (startAtom == 0) startAtom = input_.findStart();
			for (AtomIndex i = 1; i < dataSize_; i++) {
				if (compare(input_[startAtom], data_[i], exact) == false) {
					continue;
				}
				if (searchTry(startAtom, i, exact) == true)
					return true;
				//prepareToSearch();
			}
			return false;
		}
	private:
		// Node comparision
		bool compare(const RequestNodeType& inputNode, const DatabaseNodeType& dataNode, const bool exact) const noexcept {
			if (compareLow(inputNode, dataNode, exact) == false)
				return false;

			auto si = inputNode.neighboursSize();
			auto sn = dataNode.neighboursSize();

			AtomIndex j = 0;
			for (AtomIndex i = 0; i < si; ++i) {
				if (inputNode.getNeighbour(i)->getType().get_simple() < 0)
					continue;
				bool condition = false;
				for (; j < sn; ++j) {
					condition = compareLow(*(inputNode.getNeighbour(i)), *(dataNode.getNeighbour(j)), exact);
					if (condition) {
						break;
					}
				}
				if (condition == false) return false;
			}
			return true;
		}
		inline bool compareLow(const RequestNodeType& inputNode, const DatabaseNodeType& dataNode, const bool exact) const noexcept {
			if (exact) {
				return inputNode == dataNode;
			}
			else {
				return inputNode.notExactCompare(dataNode);
			}
		}

		void prepareHAtoms() {

			currents::TypeBitset bits;

			for (AtomIndex i = 1; i < inputSize_; i++)
			{
				if (input_[i].getType().include(currents::AtomTypeData(1))) {
					for (AtomIndex j = 0; j < input_[i].neighboursSize(); j++)
					{
						auto nei = input_[i].getNeighbour(j);
						bits |= nei->getType().get_bitset();
					}
				}
			}
			for (AtomIndex i = 1; i < dataSize_; i++)
			{
				if (data_[i].getType() == currents::AtomTypeData(1)) {
					for (AtomIndex j = 0; j < data_[i].neighboursSize(); j++)
					{
						auto nei = data_[i].getNeighbour(j);
						bits[nei->getType()]= true;
					}
				}
			}

			for (AtomIndex i = 1; i < inputSize_; i++)
			{
				if ((input_[i].getType().get_bitset() & bits).any()) {
					input_.unpackHydrogens(i);
				}
			}
			input_.sortGraph();
			for (AtomIndex i = 1; i < dataSize_; i++)
			{
				if (bits[data_[i].getType()]) {
					data_.unpackHydrogens(i);

				}
			}
			//data_.sortGraph();
			inputSize_ = input_.size();
			dataSize_ = data_.size();
		}


		// Search functions
		inline void addLog(const AtomIndex i1, const AtomIndex i2, const AtomIndex d1, const AtomIndex d2) {
			log_.emplace_front(BondType(i1, i2), BondType(d1, d2));
		}
		inline void deleteLog() noexcept {
			log_.pop_front();
		}
		inline void addComp(const AtomIndex i1, const AtomIndex d1) noexcept {
			comp_[i1] = d1;
			usedInComp_[d1] = true;
		}
		inline void deleteComp(const AtomIndex i1) {
			auto d1 = comp_[i1];
			comp_[i1] = 0;
			usedInComp_[d1] = false;
		}
		inline void prepareLogAndNodes(const AtomIndex cur, const AtomIndex next) {
			input_.deleteBond(cur, next);
			data_.deleteBond(comp_[cur], comp_[next]);
			addLog(cur, next, comp_[cur], comp_[next]);
		}
		inline void reverseLogAndNodes(const AtomIndex cur, const AtomIndex next) {
			input_.addBond(cur, next);
			data_.addBond(comp_[cur], comp_[next]);
			deleteLog();
		}

		inline AtomIndex findUnvisitedAtom() const noexcept {
			using stype = typename CompareVectorType::size_type;
			const stype s = comp_.size();
			for (stype i = 1; i < s; i++) {
				if (comp_[i] == 0)
					return i;
			}
			return 0;
		}
		inline AtomIndex findAtomWithNeighbours() const noexcept {
			for (auto& l : log_) {
				if (input_[l.first.first].hasNeighbours())
					return l.first.first;
			}
			for (AtomIndex i = 1; i < inputSize_; i++) {
				if (input_[i].hasNeighbours())
					return i;
			}
			return 0;
		}
		constexpr bool isDataAtomCompared(const AtomIndex dataAtomId) const noexcept {
			for (AtomIndex i = 1; i < inputSize_; i++) {
				if (comp_[i] == dataAtomId)
					return true;
			}
			return false;
		}
		bool FinalComparision(const bool exact) {
			// Preparation reversed comp_
			std::vector<bool> isDataAtomWasCompared(dataSize_, false);
			for (AtomIndex i = 1; i < inputSize_; ++i) {
				if (comp_[i] != 0)
					isDataAtomWasCompared[comp_[i]] = true;
			}
			AtomIndex j = 1;
			for (AtomIndex i = 1; i < inputSize_; ++i) {
				if (comp_[i] != 0) {
					continue;
				}

				bool condition = false;
				for (; j < dataSize_; ++j) {
					if (isDataAtomWasCompared[j])
						continue;
					if (compareLow(input_[i], data_[j], exact)) {
						condition = true;
						break;
					}
				}
				if (condition == false) return false;
			}
			return true;
		}
		bool recursiveSearchNoNeighbours(const bool exact) {
			// Find visited atom with neighbours;
			AtomIndex nextI = findAtomWithNeighbours();
			// If failed - search completed or there is atom without neighbours
			if (nextI == 0) {
				bool completed = FinalComparision(exact);
				if (completed)
					return true;
				else {
					for (size_t i = 1; i < inputSize_; i++)
					{
						if (comp_[i] == 0) {
							nextI = i;
							break;
						}
					}
				}
			}

			// nextI is an atom with neighbours
			if (comp_[nextI] != 0) {
				return recursiveSearchHasNeighbours(nextI, exact);
			}
			// need to find Twin for nextI
			for (AtomIndex i = 1; i < dataSize_; i++) {
				if ((usedInComp_[i] == true) || compare(input_[nextI], data_[i], exact) == false)
					continue;

				addComp(nextI, i);
				bool result;
				if (input_[nextI].hasNeighbours()) {
					result = recursiveSearchHasNeighbours(nextI, exact);
				}
				else {
					result = recursiveSearchNoNeighbours(exact);
				}
				if (result) return true;
				deleteComp(nextI);
			}
			return false;
		}
		bool recursiveSearchHasNeighbours(const AtomIndex curI, const bool exact) {
			_ASSERT(data_[1].getID() == 1);
			// Ring check
			const AtomIndex neiSize = input_[curI].neighboursSize();
			const AtomIndex curD = comp_[curI];
			for (AtomIndex i = 0; i < neiSize; i++) {
				const AtomIndex neiID = input_.getNeighbourId(curI, i);

				if (comp_[neiID] == 0)
					continue;
				if (data_[curD].isNeighbour(data_[comp_[neiID]])) {
					prepareLogAndNodes(curI, neiID);
					bool result;
					if (input_[curI].hasNeighbours()) {
						result = recursiveSearchHasNeighbours(curI, exact);
					}
					else {
						result = recursiveSearchNoNeighbours(exact);
					}
					if (result) return true;
					reverseLogAndNodes(curI, neiID);
				}
				return false;
			}

			// Check neighbours
			const AtomIndex nextI = input_.getNeighbourId(curI, 0);
			const AtomIndex neiDataSize = data_[curD].neighboursSize();
			for (AtomIndex i = 0; i < neiDataSize; i++) {
				const auto& neiData = data_.getNeighbourReference(curD, i);
				const AtomIndex neiID = neiData.getID();
				if (usedInComp_[neiID] == true)
					continue;
				if (compare(input_[nextI], neiData, exact) == false)
					continue;
				addComp(nextI, neiID);
				prepareLogAndNodes(curI, nextI);
				bool result;
				if (input_[nextI].hasNeighbours()) {
					result = recursiveSearchHasNeighbours(nextI, exact);
				}
				else {
					result = recursiveSearchNoNeighbours(exact);
				}
				if (result) return true;
				reverseLogAndNodes(curI, nextI);
				deleteComp(nextI);
			}
			return false;
		}
	};
}