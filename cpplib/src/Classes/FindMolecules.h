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

template<class A, class AI, class T>
struct FAM_Struct {
	// Definitions
	using AtomTypeConteinerType = std::vector<A>;
	using PointType = geometry::Point<T>;
	using PointConteinerType = std::vector<PointType>;
	using size_type = AI;
	using BondType = Bond<AI>;
	using BondExType = BondEx<AI, T>;
	using DistancesType = Distances<A, T, size_type>;
	using ParseIndexType = std::vector<size_t>;
	// Data
	AtomTypeConteinerType types;
	PointConteinerType points;
	ParseIndexType parseIndex;
	
	size_type sizeUnique = 0;
	size_type sizePoints = 0;

	// Constructors
	constexpr FAM_Struct() noexcept = default;
	constexpr FAM_Struct(AtomTypeConteinerType&& t, PointConteinerType&& p) noexcept : types(std::move(t)), points(std::move(p)) {
		sizeUnique = types.size();
		sizePoints = points.size();
		parseIndex.resize(sizeUnique);
		std::iota(parseIndex.begin(), parseIndex.end(), 0); // Fill with 0, 1...
	};

	// Methods
	constexpr auto findBonds(const DistancesType& distances, std::string& errorMSG, std::function<T(const PointType& p1, const PointType& p2)> distance_f) const {
		std::vector<BondType> res;
		std::vector<AI> invalidAtoms;
		for (size_type i = 0; i < sizePoints; i++) {
			const auto indexI = parseIndex[i];
			const auto type_i = types[indexI];
			for (size_type j = i + 1; j < sizePoints; j++) {
				const auto indexJ = parseIndex[j];
				const auto type_j = types[indexJ];
				T dist = distance_f(points[i], points[j]);
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
	constexpr auto findBondsEx(const DistancesType& distances, std::string& errorMSG, std::function<T(const PointType& p1, const PointType& p2)> distance_f) const {
		std::vector<BondExType> res;
		std::vector<AI> invalidAtoms;
		for (size_type i = 0; i < sizePoints; i++) {
			const auto indexI = parseIndex[i];
			const auto type_i = types[indexI];
			for (size_type j = i + 1; j < sizePoints; j++) {
				const auto indexJ = parseIndex[j];
				const auto type_j = types[indexJ];
				T dist = distance_f(points[i], points[j]);
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
	inline T findCutoff(const DistancesType& distances) const noexcept {
		T ret = 0;
		for (size_type i = 0; i < sizeUnique; i++) {
			auto var = distances.maxDistance(types[i], types[i]);
			if (var > ret) ret = var;
		}
		return std::fma(ret, static_cast<T>(2), static_cast<T>(0.0001));
	}

	// Debug Method
	void writeXYZ(const std::string& name) {
		std::ofstream out;
		const auto size = std::min(types.size(), points.size());
		for (size_t i = 0; i < size; i++)
		{
			out << types[i] << ' '<< points[i].get(0) << ' ' << points[i].get(1) << ' ' << points[i].get(2) << '\n';
		}
		out.close();
	}
};

template<class T>
struct FAM_Cell : public geometry::Cell<T> {
	using PointType = geometry::Point<T>;
	using PointConteinerType = std::vector<PointType>;
	using base = geometry::Cell<T>;
	explicit FAM_Cell(base&& cell) : base(std::move(cell)) {}
	template <class A, class AI> void GenerateSymm(FAM_Struct<A,AI,T>& fs, const std::vector<geometry::Symm<T>>& symm, const bool intoCell = true) const {
		const size_t p_s = fs.points.size();
		const size_t s_s = symm.size();

		// Find all translated atoms in "unique" atoms
		std::vector<bool> unique = this->FindUnique(fs, p_s, intoCell);

		for (size_t p = 0; p < p_s; p++) {
			if (unique[p] == false) // skip non unique atoms
				continue;
			const size_t p_start = fs.points.size();
			for (size_t s = 0; s < s_s; s++) {
				PointType newpoint = intoCell ? symm[s].GenSymmNorm(fs.points[p]) : symm[s].GenSymm(fs.points[p]);
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
	void CreateSupercell(PointConteinerType& points, T cutoff, uint8_t minimum = 1) {

		auto super = base::template findOptimalSupercell<uint8_t>(cutoff, minimum);
		// super_ is ready. Next step is resizing of actual points
		uint8_t mult_super_cell = super.get(2) * super.get(1) * super.get(0);
		size_t sizePoints = points.size();
		points.reserve(static_cast<typename PointConteinerType::size_type>(sizePoints) * mult_super_cell);
		for (uint8_t i = 0; i < 3; i++) {
			uint8_t cur_mult_minus1 = (super.get(i) - uint8_t(1));
			if (cur_mult_minus1 == 0)
			{
				continue;
			}
			uint8_t div2 = cur_mult_minus1 >> uint8_t(1);
			sizePoints = points.size();
			for (uint8_t j = 1; j <= div2; j++) {
				for (size_t k = 0; k < sizePoints; k++)
				{
					points.emplace_back(points[k]);
					points.back().set(i, (points.back().get(i) + static_cast<T>(j)));
					points.emplace_back(points[k]);
					points.back().set(i, (points.back().get(i) - static_cast<T>(j)));
				}
			}

			if ((cur_mult_minus1 & 1) == 1) {
				for (size_t k = 0; k < sizePoints; k++)
				{
					points.emplace_back(points[k]);
					if (points[k].get(i) < 0.5)
					{
						points.back().set(i, (points.back().get(i) + static_cast<T>(div2 + 1)));
					}
					else {
						points.back().set(i, (points.back().get(i) - static_cast<T>(div2 + 1)));
					}
				}
			}
			sizePoints = points.size();
			// Shrink cell
			for (size_t j = 0; j < sizePoints; j++) {
				points[j].set(i, (points[j].get(i) - (static_cast<T>(0.5))) / static_cast<T>(super.get(i)) + static_cast<T>(0.5));
			}
			base::lat_dir(i) *= super.get(i);
		}
		// Update base cell
		base::create(base::lat_dir(0), base::lat_dir(1), base::lat_dir(2), base::getAngleGrad(0), base::getAngleGrad(1), base::getAngleGrad(2), true);
	}
	T distanceInCell(const PointType& p1, const PointType& p2) const noexcept {
		auto dp = (p1 - p2).MoveToCell();
		T ret = 0;
		for (uint8_t i = 0; i < 3; i++) {
			T val = dp.get(i);
			if (val > 0.5)
				dp.set(i, val - 1);
		}
		return (base::fracToCart() * dp).r();
	}
private:
	template <class A, class AI> inline std::vector<bool> FindUnique(const FAM_Struct<A,AI,T>& fs, const size_t p_s, const bool intoCell) const {
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
		constexpr T distanceCutoff = 0.1;
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

template<class A, class H, class AI, class T> class FindMolecules {
public:
	using FAMSType = FAM_Struct<A, AI, T>;
	using AtomTypeConteinerType = std::vector<A>;
	using PointType = geometry::Point<T>;
	using PointConteinerType = std::vector<PointType>;
	using size_type = AI;
	using BondType = Bond<AI>;
	using DistancesType = Distances<A, T, size_type>;
	using NodeType = Node<A, H, AI>;
	using HashType = Hash<A, H, AI>;
private:
	FAMSType fs_;

public:
	constexpr FindMolecules() noexcept = delete;
	explicit constexpr FindMolecules(FAMSType && fs) : fs_(std::move(fs)) {}
	std::pair<std::string, std::vector<std::pair<std::vector<std::tuple<PointType, AI>>, int>>> findMolecules(const DistancesType& distances, std::vector<BondType>& bonds, const std::vector<AI>& invalids, std::string& errorMsg) {
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
			const A& type = net[i].getType();
			if ((type <= A(10) && contacts > 8) ||
				(type > A(10) && contacts > 14)) {
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
		std::vector<std::pair<std::vector<AI>, int>> molecules;
		std::vector<HashType> hashs;
		for (AI i = 0; i < fs_.sizeUnique; i++) {
			if (net[i].getType() == A(1) || seen[i] == true)
				continue;
			seen[i] = true;
			auto singleTable = findNextMolecule(i, net, seen);

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
					molecules[j].second++;
					break;
				}
			}
			if (!exist) {
				molecules.emplace_back(std::move(singleTable), 1);
				hashs.emplace_back(hash);
			}
		}
		auto molecules_size = molecules.size();
		std::vector<std::pair<std::vector<std::tuple<PointType, AI>>, int>> right(molecules_size);
		for (size_t i = 0; i < molecules_size; i++)
		{
			const auto mis = molecules[i].first.size();
			std::vector<std::tuple<PointType, AI>> oneMol;
			for (size_t j = 0; j < mis; j++)
			{
				auto realID = molecules[i].first[j];
				oneMol.emplace_back(fs_.points[molecules[i].first[j]], fs_.parseIndex[realID]);
			}
			right.emplace_back(std::move(oneMol), molecules[i].second);
		}
		auto res = output(molecules, net);
		if (errorMsg.empty())
			return std::make_pair(output(molecules, net), std::move(right));
		else
			return std::make_pair(output(molecules, net) + ";" + errorMsg, std::move(right));
	}


private:
	std::string output(const std::vector<std::pair<std::vector<AI>, int>>& molecules, const std::vector<NodeType>& net) const {
		size_type n = 1;
		size_type b = 0;
		std::string res;
		std::string atoms;
		std::string bonds;
		size_type s_net = net.size();
		std::vector<AI> revers(s_net, AI(-1));
		size_type m1s = molecules.size();

		for (size_type i = 0; i < m1s; i++) {
			size_type m2s = molecules[i].first.size();
			for (size_type j = 0; j < m2s; j++) {
				if (net[molecules[i].first[j]].getType() == A(1) && net[molecules[i].first[j]].neighboursSize() == 1) {
					continue;
				}

				revers[molecules[i].first[j]] = n;
				n++;
				res += std::to_string(static_cast<int>(static_cast<char>(net[molecules[i].first[j]].getType())));
				res += " ";
				res += std::to_string(calculateHAtoms(net[molecules[i].first[j]]));
				res += " ";
			}
		}
		n--;
		for (size_type i = 0; i < s_net; i++) {
			if (revers[i] == AI(-1))
				continue;
			auto nei_size = net[i].neighboursSize();
			for (size_type j = 0; j < nei_size; j++) {
				auto neiID = net[i].getNeighbour(j)->getID();
				if (neiID < i || revers[neiID] == AI(-1))
					continue;
				res += std::to_string(revers[i]);
				res += " ";
				res += std::to_string(revers[neiID]);
				res += " ";
				b++;
			}
		}

		return std::to_string(n) + " " + std::to_string(b) + " " + res;
	}
	H calculateHAtoms(const NodeType& node) const {
		H res = 0;
		auto ns = node.neighboursSize();
		for (AI i = 0; i < ns; i++)
		{
			if (node.getNeighbour(i)->getType() == A(1) && node.getNeighbour(i)->neighboursSize() == 1) {
				res++;
			}
		}
		return res;
	}
	
	constexpr std::vector<AI> findNextMolecule(const AI start, const std::vector<NodeType>& net, std::vector<bool>& seen) const {
		// 1. Create vector with single element and define limitators
		std::vector<AI> res(1, start);
		AI low_pos = 0;

		// 2. Breadth first
		while (low_pos < res.size()) {
			AI nei_size = net[res[low_pos]].neighboursSize();
			for (AI j = 0; j < nei_size; j++) {
				const auto cur_id = net[res[low_pos]].getNeighbour(j)->getID();
				if (is_member(cur_id, res) != -1) continue;
				res.emplace_back(cur_id);
			}
			seen[fs_.parseIndex[res[low_pos]]] = true;
			low_pos++;
		}
		return res;
	}
};
