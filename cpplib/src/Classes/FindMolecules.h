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
#include "../BaseHeaders/Support.h"
#include "Distances.h"
#include "Hash.h"
#include "Geometry.h"
#include "SearchGraph.h"
#include "Engine.h"
#include <array>
#include <tuple>
#include <cassert>
#include <string>
#include <vector>
#include <functional>
#include <numeric>
#include <map>

namespace cpplib {
	struct FAM_Struct {
		// Definitions
		using AtomType = currents::AtomTypeData;
		using AtomTypeBase = currents::AtomTypeBase;
		using size_type = currents::size_type;
		using PointType = currents::PointType;
		using FloatingPointType = PointType::value_type;
		using NodeType = Node<AtomType>;
		using DistancesType = Distances;
		using BondType = Bond;
		using BondExType = BondEx;
		using AtomIndex = NodeType::AtomIndex;
		using AtomContainerType = ::std::vector<AtomTypeBase>;
		using PointConteinerType = ::std::vector<PointType>;
		using DistanceFunction = ::std::function<currents::FloatingPointType(const PointType& p1, const PointType& p2)>;
		using ShiftType = geometry::Point<int>;
		using SymmRef = unsigned int;
		using ParseIndexType = ::std::vector<std::tuple<AtomIndex, SymmRef, ShiftType>>;

		static_assert(::std::is_same_v<typename Bond::AtomIndex, typename NodeType::AtomIndex>, "Bond::AtomIndex and NodeType::AtomIndex should be the same");

		// Data
		AtomContainerType types;
		PointConteinerType points;
		ParseIndexType parseIndex;

		size_type sizeUnique = 0;
		size_type sizePoints = 0;

		// Constructors
		constexpr FAM_Struct() noexcept = default;
		constexpr FAM_Struct(AtomContainerType&& t, PointConteinerType&& p) noexcept : types(std::move(t)), points(std::move(p)) {
			sizeUnique = types.size();
			sizePoints = points.size();
			parseIndex.resize(sizeUnique);
			for (decltype(sizeUnique) i = 0; i < sizeUnique; i++)
			{
				std::get<0>(parseIndex[i]) = i;
				std::get<1>(parseIndex[i]) = static_cast<SymmRef>(0);
				std::get<2>(parseIndex[i]) = ShiftType(0,0,0);
			}
		};

		// Methods
		auto findBonds(const DistancesType& distances, std::string& errorMSG, const DistanceFunction& distance_f) const {
			::std::vector<BondType> res;
			::std::vector<AtomIndex> invalidAtoms;
			for (size_type i = 0; i < sizePoints; i++) {
				const auto indexI = std::get<0>(parseIndex[i]);
				const auto type_i = types[indexI];
				for (size_type j = i + 1; j < sizePoints; j++) {
					const auto indexJ = std::get<0>(parseIndex[j]);
					const auto type_j = types[indexJ];
					auto dist = distance_f(points[i], points[j]);
					char isbond = distances.isBond(type_i, type_j, dist);
					switch (isbond) {
					case  0:
						break;
					case -1:
						// Incorrect Bond marked, but also added to results
						if (errorMSG.empty()) {
							errorMSG = "Too short bond between ";
							errorMSG += mend[type_i];
							errorMSG += std::to_string(indexI);
							errorMSG += " and ";
							errorMSG += mend[type_j];
							errorMSG += std::to_string(indexJ);
							errorMSG += ", which is ";
							errorMSG += std::to_string(dist);
						}
						invalidAtoms.push_back(i);
						invalidAtoms.push_back(j);

						// [[fallthrowgh]]
					case  1:
						res.emplace_back(i, j);
						break;
					default:
						assert(false);
						// should not triger default case;
						break;
					}
				}
			}
			return std::make_pair(res, invalidAtoms);
		}
		auto findBondsEx(const DistancesType& distances, std::string& errorMSG, const DistanceFunction& distance_f) const {
			std::vector<BondExType> res;
			std::vector<AtomIndex> invalidAtoms;
			for (size_type i = 0; i < sizePoints; i++) {
				const auto indexI = std::get<0>(parseIndex[i]);
				const auto type_i = types[indexI];
				for (size_type j = i + 1; j < sizePoints; j++) {
					const auto indexJ = std::get<0>(parseIndex[j]);
					const auto type_j = types[indexJ];
					FloatingPointType dist = distance_f(points[i], points[j]);
					char isbond = distances.isBond(type_i, type_j, dist);
					switch (isbond) {
					case  0:
						break;
					case -1:
						// Incorrect Bond marked, but also added to results
						if (errorMSG.empty()) {
							errorMSG = "Too short bond between ";
							errorMSG += mend[type_i];
							errorMSG += std::to_string(indexI);
							errorMSG += " and ";
							errorMSG += mend[type_j];
							errorMSG += std::to_string(indexJ);
							errorMSG += ", which is ";
							errorMSG += std::to_string(dist);
						}
						invalidAtoms.push_back(i);
						invalidAtoms.push_back(j);

						// [[fallthrowgh]]
					case  1:
						res.emplace_back(i, j, dist);
						break;
					default:
						assert(false);
						// should not triger default case;
						break;
					}
				}
			}
			return std::make_pair(res, invalidAtoms);
		}
		inline FloatingPointType findCutoff(const DistancesType& distances) const noexcept {
			FloatingPointType ret = 0;
			for (size_type i = 0; i < sizeUnique; i++) {
				auto var = distances.maxDistance(types[i], types[i]);
				if (var > ret) ret = var;
			}
			return ::std::fma(ret, static_cast<FloatingPointType>(2), static_cast<FloatingPointType>(0.0001));
		}
	};

	struct FAM_Cell : public geometry::Cell<currents::FloatingPointType> {
		using FloatingPointType = currents::FloatingPointType;
		using base = geometry::Cell<FloatingPointType>;
		using PointType = geometry::Point<FloatingPointType>;
		using PointConteinerType = ::std::vector<PointType>;
		using SuperCellCounter = uint_fast8_t;
		using DimmentionType = uint_fast8_t;
		using SymmType = geometry::Symm<FloatingPointType>;
		using DistancesType = Distances;
		using AtomIndex = currents::AtomIndex;

		explicit FAM_Cell(base&& cell) : base(::std::move(cell)) {}
		void GenerateSymm(FAM_Struct& fs, const std::vector<SymmType>& symm, const bool intoCell, const bool force_unique) const {
			const size_t p_s = fs.points.size();
			const FAM_Struct::SymmRef s_s = symm.size();

			// Find all translated atoms in "unique" atoms
			std::vector<bool> unique(p_s, true);
			if (force_unique) unique = this->FindUnique(fs, p_s, intoCell);

			for (FAM_Struct::SymmRef s = 0; s < s_s; s++) {
				if (symm[s].is_Eq()) continue;
				for (FAM_Struct::AtomIndex p = 0; p < p_s; p++) {
					if (unique[p] == false) // skip non unique atoms
						continue;

					PointType newpoint = symm[s].GenSymm(fs.points[p]);
					PointType shift = (PointType(0.5, 0.5, 0.5) - newpoint).round();
					if (intoCell) newpoint += shift;

					if (isTheSame_(fs.points, newpoint, intoCell)) {
						continue;
					}
					if (force_unique) {
						const size_t au = isAnotherUnique(fs.points, newpoint, p, fs.sizeUnique, intoCell);
						if (au != static_cast<size_t>(-1)) {
							unique[au] = false;
						}
					}
					fs.points.push_back(newpoint);
					fs.parseIndex.emplace_back(p, s, FAM_Struct::ShiftType(static_cast<int>(shift.get(0)), 
																		   static_cast<int>(shift.get(1)),
																		   static_cast<int>(shift.get(2))));
				}
			}
			fs.sizePoints = fs.points.size();
			if (force_unique == false) return;
			// sort and recount unique atoms
			// 1. Recount unique atoms and create shifting vector
			size_t new_count = 0;
			std::vector<size_t> shift(p_s, static_cast<size_t>(-1));
			for (size_t i = 0; i < p_s; i++)
			{
				if (unique[i]) {
					shift[i] = new_count;
					if (i != new_count) {
						fs.points[new_count] = fs.points[i];
						fs.types[new_count] = fs.types[i];
					}
					new_count++;
				}
			}
			fs.sizeUnique = new_count;
			if (new_count == p_s) {
				fs.sizePoints = fs.points.size();
				return;
			}

			// Creation and unique reorder completed
			// 2. Reorder else and change parseIndex
			const auto ds = p_s - new_count;
			const auto fps = fs.points.size();
			for (size_t i = p_s; i < fps; i++)
			{
				fs.points[i - ds] = fs.points[i];
				fs.parseIndex[i - ds] = fs.parseIndex[shift[std::get<0>(fs.parseIndex[i])]];
			}
			fs.points.resize(fps - ds);
			fs.types.resize(new_count);
			fs.parseIndex.resize(fps - ds);
			fs.sizePoints = fps - ds;
		}
		void CreateSupercell(PointConteinerType& points, FloatingPointType cutoff, SuperCellCounter minimum = 1) {

			auto super = base::template findOptimalSupercell<SuperCellCounter>(cutoff, minimum);
			// super_ is ready. Next step is resizing of actual points
			SuperCellCounter mult_super_cell = super.get(2) * super.get(1) * super.get(0);
			size_t sizePoints = points.size();
			points.reserve(static_cast<typename PointConteinerType::size_type>(sizePoints) * mult_super_cell);
			for (DimmentionType i = 0; i < static_cast<DimmentionType>(3); i++) {
				SuperCellCounter cur_mult_minus1 = (super.get(i) - SuperCellCounter(1));
				if (cur_mult_minus1 == 0)
				{
					continue;
				}
				SuperCellCounter div2 = cur_mult_minus1 >> SuperCellCounter(1);
				sizePoints = points.size();
				for (SuperCellCounter j = 1; j <= div2; j++) {
					for (size_t k = 0; k < sizePoints; k++)
					{
						points.emplace_back(points[k]);
						points.back().set(i, (points.back().get(i) + static_cast<FloatingPointType>(j)));
						points.emplace_back(points[k]);
						points.back().set(i, (points.back().get(i) - static_cast<FloatingPointType>(j)));
					}
				}

				if ((cur_mult_minus1 & 1) == 1) {
					for (size_t k = 0; k < sizePoints; k++)
					{
						points.emplace_back(points[k]);
						if (points[k].get(i) < 0.5)
						{
							points.back().set(i, (points.back().get(i) + static_cast<FloatingPointType>(div2 + 1)));
						}
						else {
							points.back().set(i, (points.back().get(i) - static_cast<FloatingPointType>(div2 + 1)));
						}
					}
				}
				sizePoints = points.size();
				// Shrink cell
				for (size_t j = 0; j < sizePoints; j++) {
					points[j].set(i, (points[j].get(i) - (static_cast<FloatingPointType>(0.5))) / static_cast<FloatingPointType>(super.get(i)) + static_cast<FloatingPointType>(0.5));
				}
				base::lat_dir(i) *= super.get(i);
			}
			// Update base cell
			base::create(base::lat_dir(0), base::lat_dir(1), base::lat_dir(2), base::getAngleGrad(0), base::getAngleGrad(1), base::getAngleGrad(2), true);
		}
		FloatingPointType distanceInCell(const PointType& p1, const PointType& p2) const noexcept {
			PointType dp = (p1 - p2).MoveToCell();
			FloatingPointType ret = 0;
			for (DimmentionType i = 0; i < static_cast<DimmentionType>(3); i++) {
				FloatingPointType val = dp.get(i);
				if (val > 0.5)
					dp.set(i, val - 1);
			}
			return (base::fracToCart() * dp).r();
		}
		auto findMoleculesForCluster(const FAM_Struct& fs, const DistancesType& distances, bool& hasPoymer) const {
			using ShiftType = FAM_Struct::ShiftType;
			using MolType = std::vector<std::vector<ShiftType>>; // mol[atomIndex][0-...?]
			hasPoymer = false;
			std::vector<std::pair<Bond,ShiftType>> bonds; // Shift-type bond container
			std::list<AtomIndex> polis;

			std::vector<AtomIndex> ref(fs.sizePoints, 0);
			std::vector<std::pair<MolType,bool>> allMolecules(1, std::make_pair(std::vector<std::vector<ShiftType>>(fs.sizePoints),false)); // Molecules start from [1]. Value [0] is always empty

			for (AtomIndex i = 0; i < fs.sizePoints; i++) {
				const auto type_i = fs.types[std::get<0>(fs.parseIndex[i])];
				for (AtomIndex j = i + 1; j < fs.sizePoints; j++) {
					const auto type_j = fs.types[std::get<0>(fs.parseIndex[j])];
					FloatingPointType dist = distanceInCell(fs.points[i], fs.points[j]);
					char isbond = distances.isBond(type_i, type_j, dist);
					switch (isbond) {
					case -1:
						//[[fallthrough]]
					case  1:
					{
						bonds.emplace_back(Bond(i, j), toShift((fs.points[i] - fs.points[j]).round()));
						const auto & curbond = bonds.back();
						AtomIndex k1 = 0;
						for (; k1 < allMolecules.size(); k1++)
						{
							if (!allMolecules[k1].first[i].empty()) break;
						}
						bool k1f = k1 != allMolecules.size(); // k1 found
						AtomIndex k2 = 0;
						for (; k2 < allMolecules.size(); k2++)
						{
							if (!allMolecules[k2].first[j].empty()) break;
						}
						bool k2f = k2 != allMolecules.size(); // k2 found

						switch ((k1f ? 1 : 0) + (k2f ? 2 : 0))
						{
						case 0: // None
							allMolecules.emplace_back(std::vector<std::vector<ShiftType>>(fs.sizePoints), false);
							allMolecules.back().first[i].push_back(std::get<2>(fs.parseIndex[i]));
							allMolecules.back().first[j].push_back(std::get<2>(fs.parseIndex[i]) + curbond.second);
							break;
						case 1: // Only k1 found
							_ASSERT(!allMolecules[k1].first[i].empty());
							allMolecules[k1].first[j].push_back(allMolecules[k1].first[i][0] + curbond.second);
							_ASSERT(allMolecules[k1].first[j].size() == 1);
							break;
						case 2: // Only k2 found
							_ASSERT(!allMolecules[k2].first[j].empty());
							allMolecules[k2].first[i].push_back(allMolecules[k2].first[j][0] - curbond.second);
							_ASSERT(allMolecules[k2].first[i].size() == 1);
							break;
						case 3: // Both found
							if (k1 != k2) { // different molecules
								Merge(fs, allMolecules[k1].first, allMolecules[k2].first, curbond);
								allMolecules.erase(allMolecules.begin() + k1);
							}
							else {
								if (!(allMolecules[k1].first[i][0] + curbond.second == allMolecules[k2].first[j][0])) { // It's a polymer!
									allMolecules[k1].first[j].push_back(allMolecules[k1].first[i][0] + curbond.second);
									allMolecules[k1].second = true;
									hasPoymer = true;
								}
							}
							break;
						}
					}
						break;
					default:
						break;
					}
				}
			}
			return allMolecules;
		}
		static inline currents::FAMStructType::ShiftType toShift(const PointType& a1) {
			using namespace currents;
			return FAMStructType::ShiftType(static_cast<FAMStructType::ShiftType::value_type>((a1.get(0))),
											static_cast<FAMStructType::ShiftType::value_type>((a1.get(1))),
											static_cast<FAMStructType::ShiftType::value_type>((a1.get(2))));
		}

	private:

		inline void Merge(const FAM_Struct& fs, 
						  const std::vector<std::vector<FAM_Struct::ShiftType>>& molO /*i*/,
						  std::vector<std::vector<FAM_Struct::ShiftType>>& molN /*j*/,
						  const std::pair<Bond, FAM_Struct::ShiftType>& bond) const {
			// calculate shift
			FAM_Struct::ShiftType s_o, s_n, s;
			if (molO[bond.first.first].empty()) {
				s_n = molN[bond.first.first][0];
				s_o = molO[bond.first.second][0];
				s = s_n - s_o + bond.second;
			}
			else {
				s_o = molO[bond.first.first][0];
				s_n = molN[bond.first.second][0];
				s = s_n - s_o - bond.second;
			}

			for (AtomIndex i = 0; i < fs.sizePoints; i++)
			{
				for (AtomIndex j = 0; j < molO[i].size(); j++)
				{
					molN[i].emplace_back(molO[i][j] + s);
				}
			}
		}


		inline std::vector<bool> FindUnique(const FAM_Struct& fs, const size_t p_s, const bool intoCell) const {
			std::vector<bool> unique(p_s, true);

			// Find all translated atoms in "unique" atoms
			for (size_t i = 0; i < p_s; i++)
			{
				if (unique[i] == false) continue;
				for (size_t j = i + 1; j < p_s; j++)
				{
					if (isTheSame(fs.points[i], fs.points[j], intoCell)) {
						unique[j] = false;
					}
				}
			}
			return unique;
		}
		inline bool isTheSame(const PointType& p1, const PointType& p2, const bool intoCell) const {
			constexpr FloatingPointType distanceCutoff = 0.1;
			if (intoCell) return distanceInCell(p1, p2) < distanceCutoff;
			else return (base::fracToCart() * (p2 - p1)).r() < distanceCutoff;
		}
		inline bool isTheSame_(PointConteinerType& points, const PointType& newpoint, const bool intoCell) const {
			const size_t sp = points.size();
			for (size_t i = 0; i < sp; i++)
			{
				if (isTheSame(newpoint, points[i], intoCell)) {
					return true;
				}
			}
			return false;
		}
		inline size_t isAnotherUnique(const PointConteinerType& points, const PointType& newpoint, const size_t counter, const size_t sizeUnique, const bool intoCell) const {
			for (size_t i = counter + 1; i < sizeUnique; i++)
			{
				if (isTheSame(newpoint, points[i], intoCell)) {
					return i;
				}
			}
			return static_cast<size_t>(-1);
		}
	};
	

	class FindMolecules {
	public:
		using FAMSType = FAM_Struct;

		using AtomTypeBase = FAMSType::AtomTypeBase;
		using AtomIndex = FAMSType::AtomIndex;
		using FloatingPointType = FAMSType::FloatingPointType;
		using AtomType = FAMSType::AtomType;
		using AtomTypeConteinerType = std::vector<AtomType>;
		using PointType = geometry::Point<FloatingPointType>;
		using PointConteinerType = std::vector<PointType>;
		using size_type = AtomIndex;
		using BondType = FAM_Struct::BondType;
		using DistancesType = Distances;
		using NodeType = Node<AtomType>;
		using HType = NodeType::HType;
		using HashType = Hash<AtomType>;
		using RightType = std::vector<std::tuple<std::vector<std::tuple<PointType, AtomIndex>>, int, std::vector<BondType>>>;
		using MoleculeType = std::pair<std::vector<AtomIndex>, std::vector<BondType>>;
	private:
		FAMSType fs_;

	public:
		FindMolecules() noexcept = delete;
		explicit FindMolecules(FAMSType&& fs) : fs_(std::move(fs)) {}
		std::tuple<std::string, std::string, RightType> findMolecules(const DistancesType& distances, std::vector<BondType>& bonds, const std::vector<AtomIndex>& invalids, std::string& errorMsg) {
			std::vector<NodeType> net;
			net.reserve(fs_.sizePoints);

			// 1. Create Nodes in net
			for (size_type i = 0; i < fs_.sizePoints; i++) {
				net.emplace_back(static_cast<typename decltype(net)::value_type::AtomType>(fs_.types[i]), 0, i);
			}

			// 2. Add bonds to net
			auto bond_size = bonds.size();
			for (decltype(bond_size) i = 0; i < bond_size; i++) {
				net[bonds[i].first].addBondWithSort(net[bonds[i].second]);
			}

			// 3. Add negative atoms
			std::vector<bool> negative_atoms(fs_.sizePoints, false);
			for (size_t i = 0; i < invalids.size(); i++)
			{
				negative_atoms[invalids[i]] = true;
			}
			for (size_type i = 0; i < fs_.sizePoints; i++) {
				auto contacts = net[i].neighboursSize();
				const AtomType& type = net[i].getType();
				if ((type <= AtomType(10) && contacts > 8) ||
					(type > AtomType(10) && contacts > 14)) {
					negative_atoms[i] = true;
					if (errorMsg.empty()) {
						errorMsg = "Atom ";
						errorMsg += mend[static_cast<AtomTypeBase>(type)];
						errorMsg += std::to_string(i);
						errorMsg += " has too many bonds (";
						errorMsg += std::to_string(contacts);
						errorMsg += ')';
					}
				}
			}

			// 4. Find molecules

			std::vector<bool> seen(fs_.sizeUnique, false); // seen unique atoms
			std::vector<std::tuple<std::vector<AtomIndex>, int, std::vector<BondType>>> molecules;
			std::vector<HashType> hashs;
			for (AtomIndex i = 0; i < fs_.sizeUnique; i++) {
				if (net[i].getType() == AtomType(1) || seen[i] == true)
					continue;
				seen[i] = true;
				auto temp_pair = findNextMolecule(i, net, seen);
				auto singleTable = std::move(temp_pair.first);

				// 4.1. 
				auto singleTableSize = singleTable.size();
				bool skip = false;
				for (size_type j = 0; j < singleTableSize; j++) {
					if (negative_atoms[singleTable[j]] == true)
						skip = true;
				}
				if (skip == true)
					continue;

				// singleTable contains AIs of single molecule
				HashType hash(singleTable, net);

				bool exist = false;
				auto molecules_size = molecules.size();
				for (size_type j = 0; j < molecules_size; j++) {
					if (hash == hashs[j]) {
						exist = true;
						(std::get<1>(molecules[j]))++;
						break;
					}
				}
				if (!exist) {
					molecules.emplace_back(std::move(singleTable), 1, std::move(temp_pair.second));
					hashs.emplace_back(hash);
				}
			}
			auto molecules_size = molecules.size();
			RightType right(0);
			for (size_t i = 0; i < molecules_size; i++)
			{
				const auto mis = std::get<0>(molecules[i]).size();
				std::vector<std::tuple<PointType, AtomIndex>> oneMol;
				for (size_t j = 0; j < mis; j++)
				{
					auto realID = std::get<0>(molecules[i])[j];
					oneMol.emplace_back(fs_.points[std::get<0>(molecules[i]).operator[](j)], std::get<0>(fs_.parseIndex[realID]));
				}
				right.emplace_back(std::move(oneMol), std::get<1>(molecules[i]), std::move(std::get<2>(molecules[i])));
			}
			auto outputStr = output(molecules, net);
			auto res = MoleculeParser<AtomType>::ResortString(outputStr.c_str()).substr(2);
			return std::make_tuple(res, errorMsg, std::move(right));
		}
		std::tuple<std::vector<AtomIndex>, std::vector<MoleculeType>> separateGraphs(std::vector<BondType>& bonds) {
			std::vector<AtomIndex> refs(fs_.sizeUnique, 0);
			std::vector<MoleculeType> molecules(1);
			std::vector<Node<currents::AtomTypeData>> nodes;
			nodes.reserve(fs_.sizePoints);

			for (AtomIndex i = 0; i < fs_.sizePoints; i++)
			{
				nodes.emplace_back(currents::AtomTypeData(fs_.types[std::get<0>(fs_.parseIndex[i])]), 0, i);
			}
			for (AtomIndex i = 0; i < bonds.size(); i++)
			{
				nodes[bonds[i].first].addBondSimple(nodes[bonds[i].second]);
			}
			std::vector<bool> seen(fs_.sizeUnique, false);
			
			for (AtomIndex i = 0; i < fs_.sizeUnique; i++)
			{
				if (seen[i]) continue;
				molecules.emplace_back(findNextMolecule(i, nodes, seen));
				for (AtomIndex j = 0; j < molecules.back().first.size(); j++)
				{
					refs[molecules.back().first[j]] = molecules.size() - 1;
				}
			}
			return std::make_tuple(refs, molecules);
		}
		PointConteinerType& compaq(std::vector<BondType>& bonds) {

			// 1. Find closest atoms
			deb_write("FM::compaq Phase 1. Find closest atoms");
			std::vector<AtomIndex> closest(fs_.sizeUnique, 0);
			std::iota(closest.begin(), closest.end(), 0);
			for (AtomIndex i = 0; i < fs_.sizePoints; i++)
			{
				constexpr PointType m(0.5, 0.5, 0.5);
				auto& cur = closest[std::get<0>(fs_.parseIndex[i])];
				if ((fs_.points[i] - m).r() < (fs_.points[cur] - m).r()) {
					cur = i;
				}
			}

			// 2. Create Nodes in net
			deb_write("FM::compaq Phase 2. Create Nodes in net");
			std::vector<NodeType> net;
			net.reserve(fs_.sizePoints);
			for (size_type i = 0; i < fs_.sizePoints; i++) {
				net.emplace_back(NodeType::AtomType(fs_.types[i]), 0, i);
			}

			// 3. Add bonds to net
			deb_write("FM::compaq Phase 3. Add bonds to net");
			auto bond_size = bonds.size();
			for (decltype(bond_size) i = 0; i < bond_size; i++) {
				net[bonds[i].first].addBondSimple(net[bonds[i].second]);
			}

			// 4. Find molecules
			deb_write("FM::compaq Phase 4. Find molecules");
			std::vector<bool> seen(fs_.sizeUnique, false); // seen unique atoms

			for (AtomIndex i = 0; i < fs_.sizeUnique; i++) {
				if (seen[i] == true)
					continue;
				seen[i] = true;

				deb_write("FM::compaq Phase 4.0. invoke FM::findNextUniquePart");
				std::vector<AtomIndex> singleTable = findNextUniquePart(closest[i], net, seen);
				deb_write("FM::compaq Phase 4.0. FM::findNextUniquePart successful");

				auto singleTableSize = static_cast<size_type>(singleTable.size());

				// 4.1. Shift Center of Mass
				deb_write("FM::compaq Phase 4.1. Shift Center of Mass");
				deb_write("Center of Mass");
				PointType center(0., 0., 0.);
				for (size_type j = 0; j < singleTableSize; j++)
				{
					center += fs_.points[singleTable[j]];
				}
				center /= fs_.sizeUnique;

				PointType shift = (PointType(0.5, 0.5, 0.5) - center).round();
				if (shift.r() > 0.5) {
					for (size_type j = 0; j < singleTableSize; j++)
					{
						fs_.points[singleTable[j]] += shift;
					}
				}

				// 4.2. Swap coordinates

				deb_write("FM::compaq Phase 4.2. Swap coordinates");
				for (size_type j = 0; j < singleTableSize; j++) {
					if(singleTable[j] != std::get<0>(fs_.parseIndex[singleTable[j]]))
						std::swap(fs_.points[singleTable[j]], fs_.points[std::get<0>(fs_.parseIndex[singleTable[j]])]);
				}
			}

			return fs_.points;
		}

	private:
		std::string output(const std::vector<std::tuple<std::vector<AtomIndex>, int, std::vector<BondType>>>& molecules, const std::vector<NodeType>& net) const {
			size_type n = 1;
			size_type b = 0;
			std::string res;
			std::string atoms;
			std::string bonds;
			size_type s_net = net.size();
			std::vector<AtomIndex> revers(s_net, AtomIndex(-1));
			size_type m1s = molecules.size();

			for (size_type i = 0; i < m1s; i++) {
				size_type m2s = std::get<0>(molecules[i]).size();
				for (size_type j = 0; j < m2s; j++) {
					if (net[std::get<0>(molecules[i])[j]].getType() == AtomType(1) && net[std::get<0>(molecules[i])[j]].neighboursSize() == 1) {
						continue;
					}

					revers[std::get<0>(molecules[i])[j]] = n;
					n++;
					res += std::to_string(static_cast<int>(static_cast<AtomTypeBase>(net[std::get<0>(molecules[i])[j]].getType())));
					res += " ";
					res += std::to_string(calculateHAtoms(net[std::get<0>(molecules[i])[j]]));
					res += " ";
				}
			}
			n--;
			for (size_type i = 0; i < s_net; i++) {
				if (revers[i] == AtomIndex(-1))
					continue;
				auto nei_size = net[i].neighboursSize();
				for (size_type j = 0; j < nei_size; j++) {
					auto neiID = net[i].getNeighbour(j)->getID();
					if (neiID < i || revers[neiID] == AtomIndex(-1))
						continue;
					res += std::to_string(revers[i]);
					res += " ";
					res += std::to_string(revers[neiID]);
					res += " ";
					b++;
				}
			}

			return std::string("0 ") + std::to_string(n) + " " + std::to_string(b) + " " + res;
		}
		HType calculateHAtoms(const NodeType& node) const {
			HType res = 0;
			auto ns = node.neighboursSize();
			for (AtomIndex i = 0; i < ns; i++)
			{
				if (node.getNeighbour(i)->getType() == AtomType(1) && node.getNeighbour(i)->neighboursSize() == 1) {
					res++;
				}
			}
			return res;
		}	
		std::pair<std::vector<AtomIndex>, std::vector<BondType>> findNextMolecule(const AtomIndex start, const std::vector<NodeType>& net, std::vector<bool>& seen) const {
			// 1. Create vector with single element and define limitators
			std::vector<AtomIndex> res(1, start);
			std::vector<BondType> bonds(0);
			AtomIndex low_pos = 0;

			// 2. Breadth-first
			while (low_pos < res.size()) {
				AtomIndex nei_size = net[res[low_pos]].neighboursSize();
				for (AtomIndex j = 0; j < nei_size; j++) {
					const auto cur_id = net[res[low_pos]].getNeighbour(j)->getID();
					auto is_m = static_cast<AtomIndex>(is_member(cur_id, res));
					if (is_m != (static_cast<AtomIndex>((size_t)(-1)))) {
						if (low_pos < is_m)
							bonds.emplace_back(low_pos, is_m);
						else {
							continue;
						}
						continue;
					}
					else {
						bonds.emplace_back(low_pos, static_cast<AtomIndex>(res.size()));
					}
					res.emplace_back(cur_id);
				}
				seen[std::get<0>(fs_.parseIndex[res[low_pos]])] = true;
				low_pos++;
			}
			return std::make_pair(res, bonds);
		}
		std::vector<AtomIndex> findNextUniquePart(const AtomIndex start, const std::vector<NodeType>& net, std::vector<bool>& seen) {
			// 1. Create vector with single element and define limitators
			std::vector<AtomIndex> res(1, start);
			AtomIndex low_pos = 0;

			// 2. Breadth-first
			while (low_pos < res.size()) {
				AtomIndex nei_size = net[res[low_pos]].neighboursSize();
				for (AtomIndex j = 0; j < nei_size; j++) {
					const auto cur_id = net[res[low_pos]].getNeighbour(j)->getID();
					auto is_m = static_cast<AtomIndex>(is_member(cur_id, res));
					if (is_m == (static_cast<AtomIndex>((size_t)(-1))) && seen[std::get<0>(fs_.parseIndex[cur_id])] == false) {
						res.emplace_back(cur_id);
						seen[std::get<0>(fs_.parseIndex[cur_id])] = true;
						auto shift = (fs_.points[res[low_pos]]- fs_.points[cur_id]).round();
						if(shift.r() > 0.1) 
							fs_.points[cur_id] += shift;
					}
				}
				low_pos++;
			}
			return res;
		}
	};

	class HashedSpace {
	public:
		using AtomIndex = currents::AtomIndex;
		using PointType = currents::PointType;
		using CellType = currents::CellType;
		using BondList = ::std::vector<currents::BondType>;
		using FloatingPointType = typename CellType::value_type;
		using DimentionType = unsigned char;

		using SupListType = ::std::list<AtomIndex>;
		using SupType = ::std::array<::std::array<::std::array<SupListType, 3>, 3>, 3>;
		using SupPoint = geometry::Point<size_t>;

		constexpr HashedSpace(const CellType& cell, FloatingPointType maxbond) noexcept :
			cell_(cell) {
			calculateSep(maxbond);
		}
		constexpr bool is_effective() const noexcept {
			return sep_[0] > 3 || sep_[1] > 3 || sep_[2] > 3;
		}
		BondList create_hash_bonds(const ::std::vector<PointType>& points) const {
			BondList ret;
			size_t estimated_size = points.size() * points.size() * sizemod();
			ret.reserve(estimated_size);
			SupType supply_table;

			// Fill supply_table
			for (AtomIndex i = 0; i < points.size(); i++)
			{
				auto c = coordintate_of_point(points[i]);
				supply_table[c.get(0)][c.get(1)][c.get(2)].emplace_back(i);
			}

			// Create all bonds
			for (size_t i = 0; i < sep_[0]; i++) {
				for (size_t j = 0; j < sep_[1]; j++) {
					for (size_t k = 0; k < sep_[2]; k++) {
						box_working(ret, supply_table, i, j, k);
					}
				}
			}
			return ret;
		}

	private:
		void box_working(BondList& ret, const SupType& supply_table, size_t i, size_t j, size_t k) const {
			create_bonds_in_box(ret, supply_table[i][j][k]);

			bool is_x = sep_[0] != 1;
			bool is_y = sep_[1] != 1;
			bool is_z = sep_[2] != 1;

			auto dx = (i + 1 == sep_[0]) ? 0 : i + 1;
			auto dy = (j + 1 == sep_[1]) ? 0 : j + 1;
			auto dz = (k + 1 == sep_[2]) ? 0 : k + 1;

			// dx
			if (is_x) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][j][k]);
			}
			// dy
			if (is_y) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][dy][k]);
			}
			// dz
			if (is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][j][dz]);
			}

			// dxdy
			if (is_x && is_y) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][dy][k]);
			}
			// dxdz
			if (is_x && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][j][dz]);
			}
			// dydz
			if (is_y && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][dy][dz]);
			}

			// dxdydz
			if (is_x && is_y && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][dy][dz]);
			}
		}
		void create_bonds_in_box(BondList& ret, const SupListType& l) const {
			for (auto iter1 = l.begin(); iter1 != l.end(); iter1++)
			{
				auto iter2 = iter1;
				iter2++;
				for (; iter2 != l.end(); iter2++)
				{
					ret.emplace_back(*iter1, *iter2);
				}
			}
		}		
		constexpr void create_bonds_between_boxes(BondList& ret, const SupListType& l1, const SupListType& l2) const {
			for (auto v1 : l1)
			{
				for (auto v2 : l2)
				{
					ret.emplace_back(v1, v2);
					ret.back().validate();
				}
			}
		}

		static constexpr FloatingPointType modifier_ = 1.05;

		constexpr SupPoint coordintate_of_point(const PointType& p) const noexcept {
			return {
				static_cast<size_t>(std::floor(p.get(0) * sep_[0])),
				static_cast<size_t>(std::floor(p.get(1) * sep_[1])),
				static_cast<size_t>(std::floor(p.get(2) * sep_[2]))
			};
		}
		constexpr FloatingPointType sizemod() const noexcept {
			FloatingPointType ret = 1;
			for (DimentionType i = 0; i < 3; i++)
			{
				if (sep_[i] > 3) {
					ret *= 3;
					ret /= sep_[i];
					ret *= modifier_;
				}
			}
			return ret;
		}
		constexpr void calculateSep(FloatingPointType maxbond) {
			maxbond *= modifier_;
			for (DimentionType i = 0; i < 3; i++) {
				sep_[i] = static_cast<size_t>(floor(cell_.lat_dir(i) / maxbond));
				if (sep_[i] == 0) sep_[i] = 1;
			}
		}

	private:
		// Data
		const CellType& cell_;
		::std::array<size_t, 3> sep_;
	};
}