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
#include <cassert>
#include <string>
#include <vector>
#include <functional>
#include <numeric>

namespace cpplib {
	struct FAM_Struct {
		// Definitions
		using AtomType = currents::AtomTypeData;
		using size_type = currents::size_type;
		using PointType = currents::PointType;
		using FloatingPointType = PointType::value_type;
		using NodeType = Node<AtomType>;
		using DistancesType = Distances;
		using BondType = Bond;
		static_assert(::std::is_same<typename Bond::AtomIndex, typename NodeType::AtomIndex>::value, "Bond::AtomIndex and NodeType::AtomIndex should be the same");
		using BondExType = BondEx;
		using AtomIndex = NodeType::AtomIndex;
		using AtomContainerType = ::std::vector<AtomType>;
		using ParseIndexType = ::std::vector<size_type>;
		using PointConteinerType = ::std::vector<PointType>;
		using DistanceFunction = ::std::function<currents::FloatingPointType(const PointType& p1, const PointType& p2)>;
		// Data
		AtomContainerType types;
		PointConteinerType points;
		ParseIndexType parseIndex;

		size_type sizeUnique = 0;
		size_type sizePoints = 0;

		// Constructors
		inline FAM_Struct() noexcept = default;
		inline FAM_Struct(AtomContainerType&& t, PointConteinerType&& p) noexcept : types(std::move(t)), points(std::move(p)) {
			sizeUnique = types.size();
			sizePoints = points.size();
			parseIndex.resize(sizeUnique);
			::std::iota(parseIndex.begin(), parseIndex.end(), 0); // Fill with 0, 1...
		};

		// Methods
		auto findBonds(const DistancesType& distances, std::string& errorMSG, DistanceFunction distance_f) const {
			::std::vector<BondType> res;
			::std::vector<AtomIndex> invalidAtoms;
			for (size_type i = 0; i < sizePoints; i++) {
				const auto indexI = parseIndex[i];
				const auto type_i = types[indexI];
				for (size_type j = i + 1; j < sizePoints; j++) {
					const auto indexJ = parseIndex[j];
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
		auto findBondsEx(const DistancesType& distances, std::string& errorMSG, DistanceFunction distance_f) const {
			std::vector<BondExType> res;
			std::vector<AtomIndex> invalidAtoms;
			for (size_type i = 0; i < sizePoints; i++) {
				const auto indexI = parseIndex[i];
				const auto type_i = types[indexI];
				for (size_type j = i + 1; j < sizePoints; j++) {
					const auto indexJ = parseIndex[j];
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
			return std::fma(ret, static_cast<FloatingPointType>(2), static_cast<FloatingPointType>(0.0001));
		}

		// Debug Method
		void writeXYZ(const std::string& name) {
			std::ofstream out;
			const auto size = std::min(types.size(), points.size());
			for (size_t i = 0; i < size; i++)
			{
				out << types[i] << ' ' << points[i].get(0) << ' ' << points[i].get(1) << ' ' << points[i].get(2) << '\n';
			}
			out.close();
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

		explicit FAM_Cell(base&& cell) : base(::std::move(cell)) {}
		void GenerateSymm(FAM_Struct& fs, const std::vector<SymmType>& symm, const bool intoCell = true) const {
			const size_t p_s = fs.points.size();
			const size_t s_s = symm.size();

			// Find all translated atoms in "unique" atoms
			std::vector<bool> unique = this->FindUnique(fs, p_s, intoCell);

			for (size_t p = 0; p < p_s; p++) {
				if (unique[p] == false) // skip non unique atoms
					continue;
				const size_t p_start = fs.points.size();
				for (size_t s = 0; s < s_s; s++) {
					PointType newpoint = (intoCell ? symm[s].GenSymmNorm(fs.points[p]) : symm[s].GenSymm(fs.points[p]));
					if (isTheSame_(fs.points, newpoint, intoCell)) {
						continue;
					}
					const size_t au = isAnotherUnique(fs.points, newpoint, p, fs.sizeUnique, intoCell);
					if (au != static_cast<size_t>(-1)) {
						unique[au] = false;
					}
					if (intoCell) fs.points.push_back(newpoint.MoveToCell());
					else fs.points.push_back(newpoint);
					fs.parseIndex.push_back(p);
				}
			}
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
				fs.parseIndex[i - ds] = shift[fs.parseIndex[i]];
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
			auto dp = (p1 - p2).MoveToCell();
			FloatingPointType ret = 0;
			for (DimmentionType i = 0; i < static_cast<DimmentionType>(3); i++) {
				FloatingPointType val = dp.get(i);
				if (val > 0.5)
					dp.set(i, val - 1);
			}
			return (base::fracToCart() * dp).r();
		}
	private:
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
				net.emplace_back(fs_.types[i], 0, i);
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
						errorMsg += mend[static_cast<char>(type)];
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
					oneMol.emplace_back(fs_.points[std::get<0>(molecules[i]).operator[](j)], fs_.parseIndex[realID]);
				}
				right.emplace_back(std::move(oneMol), std::get<1>(molecules[i]), std::move(std::get<2>(molecules[i])));
			}
			auto outputStr = output(molecules, net);
			auto res = cpplib::currents::SearchGraphType::DatabaseGraphType::ResortString(outputStr.c_str()).substr(2);
			return std::make_tuple(res, errorMsg, std::move(right));
		}
		PointConteinerType& compaq(const DistancesType& distances, std::vector<BondType>& bonds) {

			// 1. Find closest atoms
			deb_write("FM::compaq Phase 1. Find closest atoms");
			std::vector<AtomIndex> closest(fs_.sizeUnique, 0);
			std::iota(closest.begin(), closest.end(), 0);
			for (AtomIndex i = 0; i < fs_.sizePoints; i++)
			{
				constexpr PointType m(0.5, 0.5, 0.5);
				auto& cur = closest[fs_.parseIndex[i]];
				if ((fs_.points[i] - m).r() < (fs_.points[cur] - m).r()) {
					cur = i;
				}
			}

			// 2. Create Nodes in net
			deb_write("FM::compaq Phase 2. Create Nodes in net");
			std::vector<NodeType> net;
			net.reserve(fs_.sizePoints);
			for (size_type i = 0; i < fs_.sizePoints; i++) {
				net.emplace_back(fs_.types[i], 0, i);
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

				size_type singleTableSize = static_cast<size_type>(singleTable.size());

				// 4.1. Shift Center of Mass
				deb_write("FM::compaq Phase 4.1. Shift Center of Mass");
				deb_write("Center of Mass");
				PointType center(0, 0, 0);
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
					if(singleTable[j] != fs_.parseIndex[singleTable[j]])
						std::swap(fs_.points[singleTable[j]], fs_.points[fs_.parseIndex[singleTable[j]]]);
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
					res += std::to_string(static_cast<int>(static_cast<char>(net[std::get<0>(molecules[i])[j]].getType())));
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
						bonds.emplace_back(low_pos, res.size());
					}
					res.emplace_back(cur_id);
				}
				seen[fs_.parseIndex[res[low_pos]]] = true;
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
					if (is_m == (static_cast<AtomIndex>((size_t)(-1))) && seen[fs_.parseIndex[cur_id]] == false) {
						res.emplace_back(cur_id);
						seen[fs_.parseIndex[cur_id]] = true;
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
}