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
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/DSU.h"
#include "../Classes/Distances.h"
#include "../Classes/Engine.h"
#include "../Classes/Geometry.h"

#include <iostream>

namespace cpplib {
	class Cluster;
}

namespace cpplib::cluster_detail {

	using FloatingPointType = cpplib::basic_types::FloatingPointType;
	using PointType = cpplib::geometry::Point<FloatingPointType>;
	using ShiftType = cpplib::geometry::Point<int8_t>;
	using Plane = cpplib::geometry::Plane<FloatingPointType>;
	using SymmIndex = int;
	using CellType = cpplib::geometry::Cell<FloatingPointType>;
	using Matrix = typename CellType::matrix_type;
	using AtomIndex = cpplib::basic_types::AtomIndex;
	using SymmType = cpplib::geometry::Symm<FloatingPointType>;
	using AtomTypeBase = cpplib::basic_types::AtomTypeBase;

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
	struct TranslatedItem {
		AtomIndex id = 0;
		ShiftType shift;

		constexpr bool operator==(const TranslatedItem& other) const noexcept = default;

		struct Hash {
			size_t operator()(const TranslatedItem& ta) const noexcept {
				using ST = typename ShiftType::value_type;
				return std::hash<AtomIndex>()(ta.id) ^
					(std::hash<ST>()(ta.shift[0]) << 1) ^
					(std::hash<ST>()(ta.shift[1]) << 2) ^
					(std::hash<ST>()(ta.shift[2]) << 3);
			}
		};
	};
	using TranslatedAtom = TranslatedItem;
	using TranslatedMolecule = TranslatedItem;

	using BoxSet = ::std::unordered_set<ShiftType, ShiftType::Hash>;
	using TranslatedMoleculeSet = ::std::unordered_set<TranslatedMolecule, TranslatedMolecule::Hash>;

	struct Molecule {
		::std::vector<TranslatedAtom> nodes{};
		bool is_polymer = false;
	};

	static constexpr ShiftType zeroShift(0, 0, 0);

	class UnitCellBuilder {
	public:

		// Input configuration
		struct Config {
			bool remove_duplicates = true;
			FloatingPointType duplicate_tolerance = cpplib::geometry::crystallography_eq_position_eps_fractalspace;
		};

		// Extended result
		struct BuildResult {
			std::vector<ClusterAtom> atoms;
			size_t original_asymmetric_count = 0;
			size_t generated_count = 0;
			size_t duplicate_count = 0;

			bool success() const {
				return generated_count > 0;
			}
		};

		explicit UnitCellBuilder(const ::std::vector<SymmType>& symmetries,
								 Config config = Config{true, 
								                        cpplib::geometry::crystallography_eq_position_eps_fractalspace})
			: symmetries_(symmetries), config_(config) {
		}

		BuildResult build(const std::vector<PointType>& points,
						 const std::vector<AtomTypeBase>& types) const {

			BuildResult result;
			result.original_asymmetric_count = points.size();

			// 1. Create assymetric unit
			auto asymmetric_unit = create_asymmetric_unit(points, types);

			// 2. Apply symmetry\ies
			std::vector< ClusterAtom> all_atoms;
			all_atoms.reserve(asymmetric_unit.size() * symmetries_.size());

			for (int symm_idx = 0; symm_idx < symmetries_.size(); ++symm_idx) {
				for (const auto& atom : asymmetric_unit) {
					auto transformed = apply_symmetry(atom, symmetries_[symm_idx], symm_idx);

					if (config_.remove_duplicates &&
						is_duplicate(transformed, all_atoms, config_.duplicate_tolerance)) {
						result.duplicate_count++;
						continue;
					}

					all_atoms.push_back(std::move(transformed));
					result.generated_count++;
				}
			}

			result.atoms = std::move(all_atoms);
			return result;
		}

		// ==== TEST METHODS ====

		static std::vector<ClusterAtom> create_asymmetric_unit(
			const std::vector<geometry::Point<FloatingPointType>>& points,
			const std::vector<basic_types::AtomTypeBase>& types) {

			std::vector< ClusterAtom> atoms;
			atoms.reserve(points.size());

			for (size_t i = 0; i < points.size(); ++i) {
				atoms.emplace_back(
					basic_types::AtomIndex(i),
					types[i],
					points[i],
					0,  // symm_index
					ShiftType(0, 0, 0)
				);
			}
			return atoms;
		}

		ClusterAtom apply_symmetry(
			const ClusterAtom& atom,
			const geometry::Symm<FloatingPointType>& symmetry,
			int symm_index) const {

			auto temp_point = symmetry.GenSymm(atom.point);
			auto floating_shift = -temp_point.floor();
			temp_point.MoveToCell();

			 ShiftType shift(
				static_cast< ShiftType::value_type>(floating_shift[0]),
				static_cast< ShiftType::value_type>(floating_shift[1]),
				static_cast< ShiftType::value_type>(floating_shift[2])
			);

			return {
				atom.index,
				atom.type,
				temp_point,
				symm_index,
				shift
			};
		}

	private:
		bool is_duplicate(const ClusterAtom& atom,
						 const std::vector<ClusterAtom>& existing_atoms,
						 FloatingPointType tolerance) const {
			return std::any_of(existing_atoms.begin(), existing_atoms.end(),
				[&](const  ClusterAtom& existing) {
						return geometry::Point<FloatingPointType>::distance(
							atom.point, existing.point) < tolerance;
				});
		}

		const std::vector<geometry::Symm<FloatingPointType>>& symmetries_;
		Config config_;
	}; // class UnitCellBuilder

	class ConstructMolecules {
	public:
		using MoleculeIndex = AtomIndex;
		// Input configuration
		struct Config {
			bool generate_molecules_outside = true;
		};
		struct ResultType {
			::std::vector<Molecule> molecules;
			::std::vector<TranslatedMolecule> translated_molecules;
			::std::vector<MoleculeIndex> atom_to_trmol_id;
			size_t molecules_per_unit_cell = 0;
		};

		ConstructMolecules(const CellType& unit_cell,
						   const Distances& distances) noexcept
			: cell(unit_cell), dist(distances) {
		}

		ResultType execute(const ::std::vector<ClusterAtom>& atoms_01) const {
			ResultType result;
			
			auto atom_s = atoms_01.size();


			std::function< const PointType& (const ClusterAtom&)> point_iterator_lambda = 
				[](const ClusterAtom& a)->const PointType& { return a.point; };

			geometry::HashedSpace<FloatingPointType, AtomIndex> hashed_space(cell, 4.0);
			auto bonds = hashed_space.create_hash_bonds<BondWithShift>(atoms_01, 
																	   point_iterator_lambda);
			dist.filter_bond_list(bonds, atoms_01, [this](const PointType& a, const PointType& b) {return cell.distance_in_01(a, b); });

			// Based on union-find separation
			result.molecules.resize(atom_s);
			for (int i = 0; i < atom_s; i++)
			{
				result.molecules[i].nodes.emplace_back(i, zeroShift);
			}
			result.atom_to_trmol_id.resize(atom_s);
			std::iota(result.atom_to_trmol_id.begin(), result.atom_to_trmol_id.end(), 0);

			for (auto& bond : bonds) {
				if (bond.first == bond.second)
					continue;

				// Add shift to bonds
				PointType floatshift = (atoms_01[bond.first].point - atoms_01[bond.second].point).round();
				bond.shift = ShiftType(floatshift[0], floatshift[1], floatshift[2]);
				
				unite(bond, result.molecules, result.atom_to_trmol_id);

			}
            // Remove empty molecules
			decltype(result.molecules) tempmol;
			tempmol.reserve(atom_s);
			::std::vector<MoleculeIndex> mol_update_id(atom_s, 0);
			for (size_t i = 0; i < atom_s; i++)
			{
				if (result.molecules[i].nodes.empty() == true)
					continue;

				result.molecules_per_unit_cell++;
				MoleculeIndex new_index = tempmol.size();
				tempmol.push_back(std::move(result.molecules[i]));
				mol_update_id[i] = new_index;
			}
			result.molecules = std::move(tempmol);
			
			for (size_t j = 0; j < atom_s; j++)
			{
				result.atom_to_trmol_id[j] = mol_update_id[result.atom_to_trmol_id[j]];
			}

			for (size_t i = 0; i < atom_s; i++)
			{
				auto shift = find(i, result.molecules[result.atom_to_trmol_id[i]]).shift;
				
				result.atom_to_trmol_id[i] = find_or_push_mol({result.atom_to_trmol_id[i],-shift}, result.translated_molecules);
			}

			return result;
		}

		void unite(const BondWithShift& bond, std::vector<Molecule>& m, std::vector<AtomIndex>& a_to_m) const {

			Molecule& mol_a = m[a_to_m[bond.first]];
			Molecule& mol_b = m[a_to_m[bond.second]];

			// Check mol_a and mol_b are the same molecules
			if (&mol_a == &mol_b) {
				// Is it polymer?
				auto totalshift = find(bond.second, mol_a).shift - find(bond.first, mol_a).shift - bond.shift;
				if (totalshift != zeroShift) {
					// Yes! It is polymer!
					mol_a.is_polymer = true;
				}
				return;
			}

			auto mol_a_s = mol_a.nodes.size();
			auto mol_b_s = mol_b.nodes.size();
			auto sum_of_sizes = mol_a_s + mol_b_s;
			if (mol_a.nodes.capacity() < sum_of_sizes)
				mol_a.nodes.reserve(sum_of_sizes << 1);
			mol_a.nodes.insert(mol_a.nodes.end(), mol_b.nodes.begin(), mol_b.nodes.end());
			mol_b.nodes.clear();
			mol_b.nodes.reserve(0);

			mol_a.is_polymer = mol_a.is_polymer || mol_b.is_polymer;

			// change a_to_m
			for (size_t i = mol_a_s; i < sum_of_sizes; i++)
			{
				a_to_m[mol_a.nodes[i].id] = a_to_m[bond.first];
			}

			// If bond has shift: move atoms
			const auto fshift = find(bond.first, mol_a).shift - find(bond.second, mol_a).shift;
			const auto dshift = fshift + bond.shift;
			for (size_t i = mol_a_s; i < sum_of_sizes; i++)
			{
				mol_a.nodes[i].shift += dshift;
			}

		}

	private:
		const TranslatedAtom& find(AtomIndex id, const Molecule& m) const {
			auto s = m.nodes.size();
			for (size_t i = 0; i < s; i++)
			{
				if (m.nodes[i].id == id) {
					return m.nodes[i];
				}
			}
			throw m;
		}
		MoleculeIndex find_or_push_mol(TranslatedMolecule&& item, ::std::vector<TranslatedMolecule>& mols) const {
			auto s = mols.size();
			for (MoleculeIndex i = 0; i < s; i++) {
				if (mols[i] == item) {
					return i;
				}
			}
			mols.push_back(std::move(item));
			return s;
		}
		// Data
		const CellType& cell;
		const Distances& dist;
	};




} // namespace cpplib::cluster_detail




namespace cpplib {

	class Cluster {
	public:
		using FloatingPointType = cluster_detail::FloatingPointType;
		using PointType = cluster_detail::PointType;
		using ShiftType = cluster_detail::ShiftType;
		using Plane = cluster_detail::Plane;
		using BoxSet = cluster_detail::BoxSet;
		using SymmIndex = cluster_detail::SymmIndex;
		using CellType = cluster_detail::CellType;
		using Matrix = cluster_detail::Matrix;
		using AtomIndex = cluster_detail::AtomIndex;
		using SymmType = cluster_detail::SymmType;
		using AtomTypeBase = cluster_detail::AtomTypeBase;

		using AnchorType = cluster_detail::AnchorType;
		using ClusterAtom = cluster_detail::ClusterAtom;
		using BondWithShift = cluster_detail::BondWithShift;
		using BondList = cluster_detail::BondList;
		using TranslatedAtom = cluster_detail::TranslatedAtom;
		using TranslatedMolecule = cluster_detail::TranslatedMolecule;
		using Molecule = cluster_detail::Molecule;

	public:    
		Cluster(CellType& unit_cell,
				::std::vector<SymmType>& symms,
				::std::vector<AnchorType>&& anchors_fractal,
				::std::vector<PointType>&& points,
				::std::vector<AtomTypeBase>&& types,
				FloatingPointType polymer_cutoff)
			: cell(unit_cell),
			symm(symms),
			anchors_frac(std::move(anchors_fractal)),
			asymmetric_types(std::move(types)),
			asymmetric_points(std::move(points)),
			polymer_cutoff_radius(polymer_cutoff)
		
		{
			assert(asymmetric_types.size() == asymmetric_points.size());
		}
		

	private:
		CellType& cell;
		std::vector<SymmType>& symm;
		std::vector<AnchorType> anchors_frac;
		std::vector<AtomTypeBase> asymmetric_types;
		std::vector<PointType> asymmetric_points;
		FloatingPointType polymer_cutoff_radius;
	public:
		::std::vector<ClusterAtom> execute(const Distances& distances) const
		{
			constexpr PointType zeroPoint(0, 0, 0);
			::std::array<PointType, 3> e = {
				cell.fracToCart() * ShiftType(1,0,0),
				cell.fracToCart() * ShiftType(0,1,0),
				cell.fracToCart() * ShiftType(0,0,1)
			};
			const std::array<Plane, 3> plane = {Plane(zeroPoint, e[1], e[2]), Plane(zeroPoint, e[0], e[2]), Plane(zeroPoint, e[0], e[1])};

			// 1. Fill Utit cell [0,1) with atoms
			cluster_detail::UnitCellBuilder ucb(symm);
			auto unit_cell = ucb.build(asymmetric_points, asymmetric_types);
			if(unit_cell.success() == false)
				return {}; // unsuccessful generation of unit cell.
			const auto atom_size = unit_cell.atoms.size();


			// 2. Make molecules from atoms from unit cell.
			cluster_detail::ConstructMolecules cm(cell, distances);
			auto mols = cm.execute(unit_cell.atoms);

			// 3. Grow polymers
			::std::unordered_set<TranslatedAtom, TranslatedAtom::Hash> atoms;
			for (auto& molecule : mols.molecules) {
				if (molecule.is_polymer == true) {
					atoms.merge(grow_polymer(molecule, unit_cell.atoms, plane));
				}
			}

			// 4. Find boxes for molecules
			BoxSet boxes = create_boxes(plane);
			

			// 5. Add nonpolymer molecules
			::std::unordered_set<TranslatedMolecule, TranslatedMolecule::Hash> moleculeset;

			for (size_t i = 0; i < atom_size; i++) {
				const auto& tr_mol = mols.translated_molecules[mols.atom_to_trmol_id[i]];
				const auto& molref = mols.molecules[tr_mol.id];
				// Skip polymers
				if (molref.is_polymer == true) {
					continue;
				}
				for (const auto& shift : boxes) {
					auto totalMolShift = tr_mol.shift + shift;

					if (moleculeset.contains(TranslatedMolecule{tr_mol.id, totalMolShift})) {
						continue;
					}
					if (check_molecule(totalMolShift, molref, unit_cell.atoms) == false) {
						continue;
					}

					// Let's add new molecule to the set
					moleculeset.emplace(tr_mol.id, totalMolShift);

					// Add atoms to set
					for (AtomIndex j = 0; j < molref.nodes.size(); j++)
					{
						atoms.insert(TranslatedAtom{molref.nodes[j].id, molref.nodes[j].shift + totalMolShift});
					}
				}
			}

			// 6. Create output vector
			::std::vector<ClusterAtom> ret;
			ret.reserve(atoms.size());

			for (const auto& atom : atoms) {
				ret.emplace_back(unit_cell.atoms[atom.id].index,
								 unit_cell.atoms[atom.id].type,
								 unit_cell.atoms[atom.id].point + atom.shift,
								 unit_cell.atoms[atom.id].symm,
								 unit_cell.atoms[atom.id].shift + atom.shift);
			}
			return ret;
		}

	private:
		BoxSet create_boxes(const std::array<Plane, 3>& plane) const {
			BoxSet boxes;
			auto rp111 = cell.fracToCart() * ShiftType(1, 1, 1);



			std::array<FloatingPointType, 3> dp = { plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111) };

			for (const auto& anchor : anchors_frac) {

				FloatingPointType cutoff = anchor.radius;
				constructBox(anchor, cutoff, dp, boxes, plane);
			}
			return boxes;
		}




		void constructBox(const AnchorType& anchor,
						  FloatingPointType cutoff,
						  const std::array<FloatingPointType, 3>& dp,
						  BoxSet& box, 
						  const std::array<Plane, 3>& plane) const {

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
				b[0] - std::ceil((cutoff - low[0]) / dp[0]),
				b[0] + std::ceil((cutoff - high[0]) / dp[0]),
				b[1] - std::ceil((cutoff - low[1]) / dp[1]),
				b[1] + std::ceil((cutoff - high[1]) / dp[1]),
				b[2] - std::ceil((cutoff - low[2]) / dp[2]),
				b[2] + std::ceil((cutoff - high[2]) / dp[2]), };


			for (ShiftType::value_type i = maxr[0]; i <= maxr[1]; i++) {
				for (ShiftType::value_type j = maxr[2]; j <= maxr[3]; j++) {
					for (ShiftType::value_type k = maxr[4]; k <= maxr[5]; k++) {
						box.emplace(i, j, k);
					}
				}
			}
		}



		bool check_molecule(const ShiftType& shift, 
							const Molecule& mol, 
							const ::std::vector<ClusterAtom>& unit01) const {
			for (const auto& anchor : anchors_frac) {
				for (const auto& node : mol.nodes) {
					// Calculate distance to anchor
					PointType vec = cell.fracToCart()*(unit01[node.id].point + node.shift + shift - anchor.point);

					if (vec.r() < anchor.radius) {
						return true;
					}
				}
			}
			return false;
		}


		std::unordered_set<TranslatedAtom, TranslatedAtom::Hash>
			grow_polymer(const Molecule& molecule,
						 const ::std::vector<ClusterAtom>& unit01, 
						 const std::array<Plane, 3>& plane) const
		{
			std::unordered_set<TranslatedAtom, TranslatedAtom::Hash> ret;

			BoxSet boxes;
			auto rp111 = cell.fracToCart() * ShiftType(1, 1, 1);

			std::array<FloatingPointType, 3> dp = { plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111) };

			for (const auto& anchor : anchors_frac) {
				constructBox(anchor, polymer_cutoff_radius, dp, boxes, plane);
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