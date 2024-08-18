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
			dataSize_ = data_.size();
		}

		// Input and Data should be ready
		inline void prepareToSearch() {
			comp_.assign(inputSize_, 0);
			log_.clear();
			usedInComp_.assign(dataSize_, false);
		}
		bool searchTry(AtomIndex startI, AtomIndex startD, bool exact) {
			addComp(startI, startD);
			try {
				if (input_[startI].hasNeighbours()) {
					recursiveSearchHasNeighbours(startI, exact);
				}
				else {
					recursiveSearchNoNeighbours(exact);
				}
			}
			catch (bool) {
				return true;
			}
			return false;
		}
		// destroys all data, need reinitialization!
		bool startFullSearch(const bool exact, AtomIndex startAtom = 0) {
			auto ullmann = createUllmannBitMartix(exact);
			if (ullmann.empty()) 
				return false;
			const auto inputBackup = input_.makeCopy();
			const auto dataBackup = data_.makeCopy();
			if (startAtom == 0) startAtom = input_.findStart();
			for (AtomIndex i = 1; i < dataSize_; i++) {
				if (compare(input_[startAtom], data_[i], exact) == false) {
					continue;
				}
				if (searchTry(startAtom, i, exact) == true)
					return true;
				input_ = inputBackup.makeCopy();
				data_ = dataBackup.makeCopy();
				prepareToSearch();
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
				if (static_cast<char>(inputNode.getNeighbour(i)->getType()) < 0)
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
		void recursiveSearchNoNeighbours(const bool exact) {
			// Find visited atom with neighbours;
			AtomIndex nextI = findAtomWithNeighbours();
			// If failed - search completed
			if (nextI == 0) {
				bool completed = FinalComparision(exact);
				if (completed)
					throw true;
			}

			// nextI is an atom with neighbours
			if (comp_[nextI] != 0) {
				recursiveSearchHasNeighbours(nextI, exact);
				return;
			}
			// need to find Twin for nextI
			for (AtomIndex i = 1; i < dataSize_; i++) {
				if ((usedInComp_[i] == true) || compare(input_[nextI], data_[i], exact) == false)
					continue;

				addComp(nextI, i);
				if (input_[nextI].hasNeighbours()) {
					recursiveSearchHasNeighbours(nextI, exact);
				}
				else {
					recursiveSearchNoNeighbours(exact);
				}
				deleteComp(nextI);
			}
		}
		void recursiveSearchHasNeighbours(const AtomIndex curI, const bool exact) {
			// Ring check
			const AtomIndex neiSize = input_[curI].neighboursSize();
			const AtomIndex curD = comp_[curI];
			for (AtomIndex i = 0; i < neiSize; i++) {
				const AtomIndex neiID = input_.getNeighbourId(curI, i);

				if (comp_[neiID] == 0)
					continue;
				if (data_[curD].isNeighbour(data_[comp_[neiID]])) {
					prepareLogAndNodes(curI, neiID);
					if (input_[curI].hasNeighbours()) {
						recursiveSearchHasNeighbours(curI, exact);
					}
					else {
						recursiveSearchNoNeighbours(exact);
					}
					reverseLogAndNodes(curI, neiID);
				}
				return;
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
				if (input_[nextI].hasNeighbours()) {
					recursiveSearchHasNeighbours(nextI, exact);
				}
				else {
					recursiveSearchNoNeighbours(exact);
				}
				reverseLogAndNodes(curI, nextI);
				deleteComp(nextI);
			}
		}

	private:
		std::vector<std::vector<bool>> createUllmannBitMartix(const bool exact) const {
			
			// 1. Create Trees
			std::vector<std::vector<std::vector<AtomIndex>>> ullmannTreeInput(inputSize_);
			std::vector<std::vector<std::vector<AtomIndex>>> ullmannTreeData(dataSize_);
			size_t diagonal = 1;

			// 1.1 Input tree
			for (AtomIndex i = 1; i < inputSize_; i++)
			{
				std::vector<bool> ullmannUsed(inputSize_, false);

				ullmannUsed[i] = true;
				ullmannTreeInput[i].push_back(std::vector<AtomIndex>(1, i));

				size_t j = 1;
				for (;; j++) // Level of Tree
				{
					std::vector<AtomIndex> lastRow;
					for (auto nodeID : ullmannTreeInput[i].back())
					{
						AtomIndex ns = input_[nodeID].neighboursSize();
						for (AtomIndex k = 0; k < ns; k++) //check neighbours
						{
							const auto curNei = input_[nodeID].getNeighbour(k);
							const auto curNeiID = curNei->getID();
							if (ullmannUsed[curNeiID] == false) {
								ullmannUsed[curNeiID] = true;
								lastRow.push_back(curNeiID);
							}
						}
					}
					// Check last row is empty
					if (lastRow.empty())
						break;
					// if not empty - sort and add it to tree
					std::sort(lastRow.begin(), lastRow.end());
					ullmannTreeInput[i].push_back(std::move(lastRow));
				}
				if (diagonal < j) diagonal = j;
			}

			// 1.2 Data tree
			for (AtomIndex i = 1; i < dataSize_; i++)
			{
				std::vector<bool> ullmannUsed(dataSize_, false);

				ullmannUsed[i] = true;
				ullmannTreeData[i].push_back(std::vector<AtomIndex>(1, i));


				for (size_t j = 1; j < diagonal; j++) // Level of Tree
				{
					std::vector<AtomIndex> lastRow;
					for (auto& nodeID : ullmannTreeData[i].back())
					{
						AtomIndex ns = data_[nodeID].neighboursSize();
						for (AtomIndex k = 0; k < ns; k++) //check neighbours
						{
							const auto curNei = data_[nodeID].getNeighbour(k);
							const auto curNeiID = curNei->getID();
							if (ullmannUsed[curNeiID] == false) {
								ullmannUsed[curNeiID] = true;
								lastRow.push_back(curNeiID);
							}
						}
					}
					// Check last row is empty
					if (lastRow.empty())
						break;
					// if not empty - sort and add it to tree
					std::sort(lastRow.begin(), lastRow.end());
					ullmannTreeData[i].push_back(std::move(lastRow));
				}
			}


			std::vector<std::vector<bool>> ullmannBits(inputSize_, std::vector<bool>(dataSize_, false));
			// 2. First decision, second order
			for (AtomIndex i = 1; i < inputSize_; i++)
			{
				for (AtomIndex j = 1; j < dataSize_; j++)
				{
					auto is = ullmannTreeInput[i].size();
					auto js = ullmannTreeData[j].size();
					if (((!exact)&& (is < js)) || (is ==js))
						ullmannBits[i][j] = compare(input_[i], data_[j], exact);
				}
			}

			static std::vector<int> rd(30, 0);

			// 3 Refinement
			const int maxRefine = 10;
			for (int i = 0; i < maxRefine; i++)
			{
				int changes = 0;
				for (AtomIndex x = 1; x < inputSize_; x++)
				{
					AtomIndex is_active = 0;
					for (AtomIndex y = 1; y < dataSize_; y++)
					{
						if (ullmannBits[x][y] == false)
							continue;

						const auto maxDepth = ullmannTreeInput[x].size()*0.6;
						for (int depth = 0; depth < maxDepth; depth++) {
							for (auto inputID : ullmannTreeInput[x][depth]) {
								bool good = false;
								for (auto dataID : ullmannTreeData[y][depth]) {
									if (ullmannPartialComparision(ullmannBits, inputID, dataID)) {
										good = true;
										break;
									}
								}
								if (good == false) {
									rd[depth]++;
									goto wrong;
								}
							}
						}
						is_active++;
						continue;
					wrong:
						changes++;
						ullmannBits[x][y] = false;
					}
					if (is_active == 0)
						return std::vector<std::vector<bool>>();
				}
				if (changes == 0)
					break;
			}
			return ullmannBits;
		}
	private:
		bool ullmannPartialComparision(const std::vector<std::vector<bool>>& ullmannBits, AtomIndex inputID, AtomIndex dataID) const {
			if ((ullmannBits[inputID][dataID] == false))
				return false;

			//const auto nsi = input_[inputID].neighboursSize();
			//if (nsi == 0) return true;
			//const auto nsj = data_[dataID].neighboursSize();
			//if (nsj == 0) return false;
			////AtomIndex i = 0;
			//for (AtomIndex i = 0; i < nsi; i++) {
			//	auto iNeiID = input_[inputID].getNeighbour(i)->getID();
			//	bool good = false;
			//	for (AtomIndex j = 0; j < nsj; j++)
			//	{
			//		auto neiID = data_[dataID].getNeighbour(j)->getID();
			//		
			//		if(ullmannBits[iNeiID][neiID]) {
			//			good = true;
			//			break;
			//		}
			//	}
			//	if (good == false) return false;
			//}
			return true;
		}
	};
}