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
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <unordered_set>
#include <utility>
#include <vector>

#include "BaseTypes.h"
#include "DSU.h"
#include "Distances.h"
#include "Engine.h"
#include "Geometry.h"

namespace cpplib {

	class Cluster {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;
		using ShiftType = geometry::Point<int8_t>;
		using Plane = geometry::Plane<FloatingPointType>;
		using BoxSet = ::std::unordered_set<ShiftType>;
		using SymmIndex = int;
		using CellType = geometry::Cell<FloatingPointType>;
		using Matrix = typename CellType::matrix_type;
		using AtomIndex = basic_types::AtomIndex;
		using SymmType = geometry::Symm<FloatingPointType>;
		using AtomTypeBase = basic_types::AtomTypeBase;
		using DistancesType = Distances;

		struct BondWithShift : public Bond {
			ShiftType shift{0, 0, 0};
			BondWithShift(AtomIndex a, AtomIndex b) : Bond(a, b) {}
			BondWithShift(Bond a, ShiftType b) : Bond(a), shift(b) {}
		};
		using BondList = ::std::vector<BondWithShift>;
		struct AnchorType {
			PointType point;
			FloatingPointType radius;
		};
		struct ClusterAtom {
			AtomIndex index;
			AtomTypeBase type;
			PointType point;
			SymmIndex symm;
			ShiftType shift;
		};
		struct TranslatedAtom {
			AtomIndex id;
			ShiftType shift;

			constexpr bool operator==(const TranslatedAtom& other) const noexcept = default;

			struct Hash {
				size_t operator()(const TranslatedAtom& ta) const noexcept {
					using ST = typename ShiftType::value_type;
					return std::hash<AtomIndex>()(ta.id) ^
						(std::hash<ST>()(ta.shift[0]) << 1) ^
						(std::hash<ST>()(ta.shift[1]) << 2) ^
						(std::hash<ST>()(ta.shift[2]) << 3);
				}
			};
		};
		struct Molecule {
			::std::vector<TranslatedAtom> nodes{};
			bool is_polymer = false;
		};

	private:
		CellType& cell;
		::std::vector<SymmType>& symm;
		::std::vector<AnchorType> anchors_frac;
		::std::vector<ClusterAtom> asymmetric_unit;
		FloatingPointType polymer_cutoff_radius;
		::std::array<Plane, 3> plane;
		bool has_poly = false;

	public:
		Cluster() = delete;
		Cluster(CellType& unit_cell,
				::std::vector<SymmType>& symms,
				::std::vector<AnchorType>&& anchors_fractal,
				const ::std::vector<PointType>& points,
				const ::std::vector<AtomTypeBase>& types,
				FloatingPointType polymer_cutoff)
			: cell(unit_cell), symm(symms), anchors_frac(::std::move(anchors_fractal)), polymer_cutoff_radius(polymer_cutoff)
		{
			auto s = points.size();
			assert(types.size() == s);
			asymmetric_unit.reserve(s);
			for (size_t i = 0; i < s; i++) {
				asymmetric_unit.emplace_back(AtomIndex(i), types[i], points[i], SymmIndex(0), ShiftType(0, 0, 0));
			}

			constexpr PointType zeroPoint(0, 0, 0);
			::std::array<PointType, 3> e = {
				cell.fracToCart() * ShiftType(1,0,0),
				cell.fracToCart() * ShiftType(0,1,0),
				cell.fracToCart() * ShiftType(0,0,1)
			};
			plane = { Plane(zeroPoint, e[1], e[2]), Plane(zeroPoint, e[0], e[2]), Plane(zeroPoint, e[0], e[1]) };


		}
		::std::vector<ClusterAtom> execute(const DistancesType& distances)
		{
			auto unit_01 = construct_unit_01();
			auto molecules01 = constructMoleculesInUnit01(unit_01, distances);

			auto molecule_pass = analyseMolecules(molecules01);
			auto unit01_molecule_indexes = molecule_indexes_create(molecules01, unit_01.size());
			// make Boxes
			BoxSet boxes = create_boxes();

			// Find nessesary molecules
			auto moleculeBoxes = create_molecule_boxes_nonpoly(molecules01,
															   unit_01,
															   molecule_pass,
															   boxes);

			// Grow polymers
			std::unordered_set<TranslatedAtom, TranslatedAtom::Hash> atoms;
			for (auto& molecule : molecules01) {
				if (molecule.is_polymer == true) {
					atoms.merge(grow_polymer(molecule, unit_01));
				}
			}

			// add nonpolymer molecules
			for (AtomIndex i = 0; i < moleculeBoxes.size(); i++) {
				for (const auto& shift : moleculeBoxes[i]) {
					for (const auto& atom : molecules01[i].nodes) {
						auto id = atom.id;
						auto sumshift = shift + atom.shift;
						atoms.emplace(id, sumshift);
					}
				}
			}

			// Create output vector
			::std::vector<ClusterAtom> ret;
			ret.reserve(atoms.size());

			for (const auto& atom : atoms) {
				ret.emplace_back(unit_01[atom.id].index,
								 unit_01[atom.id].type,
								 unit_01[atom.id].point + atom.shift,
								 unit_01[atom.id].symm,
								 unit_01[atom.id].shift + atom.shift);
			}
			return ret;
		}

	private:
		::std::vector<Molecule> constructMoleculesInUnit01(const ::std::vector<ClusterAtom>& unit_01, const DistancesType& distances) {
			// based on union-find
			geometry::HashedSpace<FloatingPointType, AtomIndex> hashed_space(cell, 4.0);
			::std::vector<PointType> points;
			::std::vector<AtomTypeBase> types;
			points.reserve(unit_01.size());
			types.reserve(unit_01.size());
			for (const auto& atom : unit_01) {
				points.emplace_back(atom.point);
				types.emplace_back(atom.type);
			}
			auto bonds = hashed_space.create_hash_bonds<BondWithShift>(points);

			distances.filter_bond_list(bonds, types, points, [this](const PointType& a, const PointType& b) {return cell.distance_in_01(a, b); });

			// add shift to bonds
			for (auto& bond : bonds) {
				PointType floatshift = (unit_01[bond.first].point - unit_01[bond.second].point).round();
				bond.shift = ShiftType(floatshift[0], floatshift[1], floatshift[2]);
			}

			// Construct molecules from unit cell and bonds
			auto molecules = create_molecules_near_unit01(unit_01, bonds);

			return molecules;
		}

		::std::vector<ClusterAtom> construct_unit_01() const {
			auto au_s = asymmetric_unit.size();
			size_t symm_s = symm.size();

			::std::vector<ClusterAtom> unit;
			unit.reserve(au_s * symm_s);

			for (SymmIndex i = 0; i < symm_s; i++) {
				for (AtomIndex j = 0; j < au_s; j++) {
					auto temp_point = symm[i].GenSymm(asymmetric_unit[j].point);
					auto floating_shift = -temp_point.floor();
					temp_point.MoveToCell();
					if (isPointInVector(unit, temp_point)) continue;

					ShiftType shift(static_cast<ShiftType::value_type>(floating_shift[0]),
									static_cast<ShiftType::value_type>(floating_shift[1]),
									static_cast<ShiftType::value_type>(floating_shift[2]));
					unit.emplace_back(j, asymmetric_unit[j].type, temp_point, i, shift);
				}
			}
			return unit;
		}

		inline bool isPointInVector(const ::std::vector<ClusterAtom>& vec, const PointType& point) const noexcept {
			return std::ranges::any_of(vec,
									   [&](const ClusterAtom& i) {
										   return i.point == point;
									   });
		}

		::std::vector<Molecule> create_molecules_near_unit01(const std::vector<ClusterAtom>& unit01,
															 const BondList& bonds) {
			auto atomsize = unit01.size();
			constexpr ShiftType zeroShift{ 0, 0, 0 };
			DSU dsu(atomsize);

			BondList filteredbonds;
			filteredbonds.reserve(bonds.size());

			for (const auto& bond : bonds) {

				// check if bond is correct and is in unit cell
				if (bond.first == bond.second)
					continue;

				if (bond.shift != zeroShift) {
					filteredbonds.emplace_back(bond);
					continue;
				}

				dsu.unite(bond.first, bond.second);
			}

			::std::vector<Molecule> molecules;
			molecules.reserve(dsu.get_count_components());
			const auto& components = dsu.get_components_ref();
			::std::vector<AtomIndex> unit_ref_to_molecule(atomsize);

			// Create Molecules based on components
			for (const auto& component : components) {
				if (component.empty())
					continue;

				// Temporary vector to store nodes of the molecule
				decltype(Molecule::nodes) nodes;
				nodes.reserve(component.size());
				for (const auto& i : component) {
					// Create nodes in molecules
					nodes.emplace_back(i, unit01[i].shift);

					// Add reference to molecule for each atom
					unit_ref_to_molecule[i] = molecules.size();
				}
				// Add molecule to the vector
				molecules.emplace_back(::std::move(nodes), false);
			}

			// Check bonds outside unit cell
			for (const auto& bond : filteredbonds) {
				auto a_ref = unit_ref_to_molecule[bond.first];
				auto b_ref = unit_ref_to_molecule[bond.second];

				if (a_ref == b_ref) { // is already in the same molecule

					// so it is polymer
					molecules[a_ref].is_polymer = true;
					has_poly = true;
					continue;
				}

				if (check_intersection(molecules[a_ref], molecules[b_ref])) {
					// is already combined
					continue;
				}

				if (molecules[a_ref].is_polymer || molecules[b_ref].is_polymer)
				{
					molecules[a_ref].is_polymer = true;
					molecules[b_ref].is_polymer = true;
				}

				// combine molecules
				AtomIndex current_size_a = molecules[a_ref].nodes.size();
				AtomIndex current_size_b = molecules[b_ref].nodes.size();

				molecules[a_ref].nodes.insert(molecules[a_ref].nodes.end(),
											  molecules[b_ref].nodes.begin(),
											  molecules[b_ref].nodes.end());

				molecules[b_ref].nodes.insert(molecules[b_ref].nodes.end(),
											  molecules[a_ref].nodes.begin(),
											  molecules[a_ref].nodes.begin() + current_size_a);

				// TODO: Check correctness.
				//change_shift(molecules[a_ref], current_size_a, bond.shift);
				//change_shift(molecules[b_ref], current_size_b, -bond.shift);
			}

			return molecules;
		}
		bool check_intersection(const Molecule& mol1, const Molecule& mol2) const {
			auto s1 = mol1.nodes.size();
			auto s2 = mol2.nodes.size();
			if (s2 > s1)
				return check_intersection(mol2, mol1); // Call method with swapped arguments

			// TODO: Maybe useless check
			if (s2 == 0)
				return false;

			AtomIndex id = mol1.nodes[0].id;
			return std::ranges::any_of(mol2.nodes, [id](const TranslatedAtom& node) {return id == node.id; });
		}
		void change_shift(Molecule& mol, AtomIndex startIndex, ShiftType shift) const {
			AtomIndex s = mol.nodes.size();
			for (AtomIndex i = startIndex; i < s; i++) {
				mol.nodes[i].shift += shift;
			}
		}

		BoxSet create_boxes() const {
			BoxSet boxes;
			auto rp111 = cell.fracToCart() * ShiftType(1, 1, 1);

			std::array<FloatingPointType, 3> dp = { plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111) };

			for (const auto& anchor : anchors_frac) {

				FloatingPointType cutoff = anchor.radius;
				constructBox(anchor, cutoff, dp, boxes);
			}
			return boxes;
		}




		std::vector<TranslatedAtom> analyseMolecules(const ::std::vector<Molecule>& molecules) const {
			std::vector<TranslatedAtom> ret;
			ret.resize(molecules.size(), { AtomIndex(0), ShiftType(0, 0, 0) });
			for (size_t i = 0; i < molecules.size(); i++)
			{
				ret[i].id = i;
			}

			// TODO: Need check (logic)

			for (AtomIndex i = 1; i < molecules.size(); i++) {
				auto nodeid = molecules[i].nodes[0].id;
				for (AtomIndex j = 0; j < i; j++) {
					auto iter = std::ranges::find_if(molecules[j].nodes, [nodeid](const TranslatedAtom& atom) {return atom.id == nodeid; });
					if (iter != molecules[j].nodes.end())
					{
						// TODO: Need check (logic)
						ret[i].shift = iter->shift - molecules[i].nodes[0].shift;
						ret[i].id = ret[j].id;
						break;
					}
				}
			}
			return ret;
		}
		void constructBox(const AnchorType& anchor,
						  FloatingPointType cutoff,
						  const std::array<FloatingPointType, 3>& dp,
						  BoxSet& box) const {

			ShiftType b(anchor.point.floor());
			std::array < FloatingPointType, 3> low{
				plane[0].distance(cell.fracToCart() * anchor.point) - b[0] * dp[0],
				plane[1].distance(cell.fracToCart() * anchor.point) - b[1] * dp[1],
				plane[2].distance(cell.fracToCart() * anchor.point) - b[2] * dp[2] };

			std::array < FloatingPointType, 3> high{
				dp[0] - low[0],
				dp[1] - low[1],
				dp[2] - low[2] };

			// [ -x, +x, -y, +y, -z, +z ]
			const std::array<ShiftType::value_type, 6> maxr{
				b[0] - static_cast<ShiftType::value_type>(ceil((cutoff - low[0]) / dp[0])),
				b[0] + static_cast<ShiftType::value_type>(ceil((cutoff - high[0]) / dp[0])),
				b[1] - static_cast<ShiftType::value_type>(ceil((cutoff - low[1]) / dp[1])),
				b[1] + static_cast<ShiftType::value_type>(ceil((cutoff - high[1]) / dp[1])),
				b[2] - static_cast<ShiftType::value_type>(ceil((cutoff - low[2]) / dp[2])),
				b[2] + static_cast<ShiftType::value_type>(ceil((cutoff - high[2]) / dp[2])), };


			for (char i = maxr[0]; i <= maxr[1]; i++) {
				for (char j = maxr[2]; j <= maxr[3]; j++) {
					for (char k = maxr[4]; k <= maxr[5]; k++) {
						box.emplace(i, j, k);
					}
				}
			}
		}

		::std::vector<AtomIndex> molecule_indexes_create(const std::vector<Molecule>& molecules, AtomIndex size) const {
			::std::vector<AtomIndex> ret(size, 0);
			constexpr ShiftType zeroshift(0, 0, 0);

			for (AtomIndex i = molecules.size() - 1; i > 0; i--)
			{
				for (AtomIndex j = 0; j < molecules[i].nodes.size(); j++)
				{
					if (molecules[i].nodes[j].shift == zeroshift) {
						ret[j] = i;
					}
				}
			}
			return ret;
		}
		::std::vector<BoxSet> create_molecule_boxes_nonpoly(const ::std::vector<Molecule>& molecules,
															const ::std::vector<ClusterAtom>& unit01,
															const ::std::vector<TranslatedAtom>& molecule_pass,
															const BoxSet& boxes) const
		{
			// Construct phantom of each possible molecule
			::std::vector<BoxSet> ret(molecule_pass.size());
			for (AtomIndex i = 0; i < molecule_pass.size(); i++) {
				if (molecules[i].is_polymer == true)
					continue;
				if (molecule_pass[i].id == i) {

					ret[i] = boxes;
				}
				else {
					for (const auto& box : boxes) {
						ret[molecule_pass[i].id].emplace(box + molecule_pass[i].shift);
					}
				}
			}

			// Check each combination of molecule, anchor and box
			for (AtomIndex i = 0; i < molecule_pass.size(); i++)
			{
				std::erase_if(ret[i], [this, i, &molecules, &unit01](const ShiftType& shift) {
					return !(this->check_molecule(shift, molecules[i], unit01)); });
			}

			return ret;
		}

		bool check_molecule(const ShiftType& shift, const Molecule& mol, const ::std::vector<ClusterAtom>& unit01) const {
			for (const auto& anchor : anchors_frac) {
				for (const auto& node : mol.nodes) {
					// Calculate distance to anchor
					PointType vec = unit01[node.id].point + node.shift + shift - anchor.point;

					if (vec.r() < anchor.radius) {
						return true;
					}
				}
			}
			return false;
		}


		std::unordered_set<TranslatedAtom, TranslatedAtom::Hash>
			grow_polymer(const Molecule& molecule,
						 const ::std::vector<ClusterAtom>& unit01) const
		{
			std::unordered_set<TranslatedAtom, TranslatedAtom::Hash> ret;

			BoxSet boxes;
			auto rp111 = cell.fracToCart() * ShiftType(1, 1, 1);

			std::array<FloatingPointType, 3> dp = { plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111) };

			for (const auto& anchor : anchors_frac) {
				constructBox(anchor, polymer_cutoff_radius, dp, boxes);
			}

			for (const auto& shift : boxes) {
				for (const auto& node : molecule.nodes) {
					for (const auto& anchor : anchors_frac) {
						PointType vec = unit01[node.id].point + shift - anchor.point;

						if (vec.r() < polymer_cutoff_radius) {
							ret.emplace(node.id, shift);
						}
					}
				}
			}

			return ret;
		}
	};

}