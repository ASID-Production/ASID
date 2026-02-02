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
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/Geometry.h"
#include "../Classes/Distances.h"

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
	using BondWithShift = typename geometry::SpatialGrid<FloatingPointType>::BondWithShift;
	using BondList = ::std::vector<BondWithShift>;


	class ClusterData {
	public:
		using PointType = geometry::Point<FloatingPointType>;

		std::vector<AtomIndex> indices;
		std::vector<AtomTypeBase> types;
		std::vector<PointType> points;
		std::vector<SymmIndex> symm_indices;
		std::vector<ShiftType> shifts;

		size_t size() const {
			//assert(consistency_check());
			return indices.size();
		}

		bool empty() const {
			return indices.empty();
		}

		bool consistency_check() const {
			const size_t s = indices.size();
			return s == types.size() &&
				s == points.size() &&
				s == symm_indices.size() &&
				s == shifts.size();
		}

		void reserve(size_t capacity) {
			indices.reserve(capacity);
			types.reserve(capacity);
			points.reserve(capacity);
			symm_indices.reserve(capacity);
			shifts.reserve(capacity);
		}

		void push_back(AtomIndex idx, AtomTypeBase type, const PointType& point,
					   SymmIndex symm, const ShiftType& shift) {
			indices.push_back(idx);
			types.push_back(type);
			points.push_back(point);
			symm_indices.push_back(symm);
			shifts.push_back(shift);
		}

		void clear() {
			indices.clear();
			types.clear();
			points.clear();
			symm_indices.clear();
			shifts.clear();
		}

		struct AtomView {
			AtomIndex index;
			AtomTypeBase& type;
			PointType& point;
			SymmIndex& symm;
			ShiftType& shift;
		};

		AtomView operator[](size_t i) {
			return {indices[i], types[i], points[i], symm_indices[i], shifts[i]};
		}

		struct ConstAtomView {
			AtomIndex index;
			const AtomTypeBase& type;
			const PointType& point;
			SymmIndex symm;
			const ShiftType& shift;
		};

		ConstAtomView operator[](size_t i) const {
			return {indices[i], types[i], points[i], symm_indices[i], shifts[i]};
		}
	};

	struct AnchorType {
		PointType point{};
		FloatingPointType radius = 0.0;
		AnchorType(const PointType& point, FloatingPointType radius) : point(point), radius(radius) {
			constexpr auto smin = std::numeric_limits<ShiftType::value_type>::min();
			constexpr auto smax = std::numeric_limits<ShiftType::value_type>::max();

			if (radius < 0.0) {
				throw std::runtime_error("Anchor radius must be positive");
			}


			if (point[0] < smin ||
				point[1] < smin ||
				point[2] < smin ||
				point[0] > smax ||
				point[1] > smax ||
				point[2] > smax)
			{
				throw std::runtime_error("Anchor point is out of bounds");
			}
		}
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
			ClusterData atoms;
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
			: symmetries_(symmetries), config_(config) {}

		BuildResult build(const std::vector<PointType>& points,
						 const std::vector<AtomTypeBase>& types) const {

			BuildResult result;
			result.original_asymmetric_count = points.size();

			// 1. Create assymetric unit
			const ClusterData asymmetric_unit = create_asymmetric_unit(points, types);

			// 2. Apply symmetry
			ClusterData all_atoms;
			all_atoms.reserve(asymmetric_unit.size() * symmetries_.size());

			for (int symm_idx = 0; symm_idx < (int)symmetries_.size(); ++symm_idx) {
				for (size_t i = 0; i < asymmetric_unit.size(); ++i) {
					auto transformed = apply_symmetry(asymmetric_unit[i], symmetries_[symm_idx], symm_idx);

					if (config_.remove_duplicates &&
						is_duplicate(transformed, all_atoms, config_.duplicate_tolerance)) {
						result.duplicate_count++;
						continue;
					}

					all_atoms.push_back(
						transformed.index,
						transformed.type,
						transformed.point,
						transformed.symm,
						transformed.shift
					);
					result.generated_count++;
				}
			}

			result.atoms = std::move(all_atoms);
			return result;
		}

		static ClusterData create_asymmetric_unit(
			const std::vector<geometry::Point<FloatingPointType>>& points,
			const std::vector<basic_types::AtomTypeBase>& types) {

			assert(points.size() == types.size());

			ClusterData atoms;
			atoms.indices.resize(points.size(), 0);
			atoms.points = points;
			atoms.types = types;
			atoms.symm_indices.assign(points.size(), 0);
			atoms.shifts.assign(points.size(), zeroShift);

			std::iota(atoms.indices.begin(), atoms.indices.end(), 0);

			return atoms;
		}
	private:

		struct TempAtom {
			AtomIndex index;
			AtomTypeBase type;
			PointType point;
			SymmIndex symm;
			ShiftType shift;
		};

		TempAtom apply_symmetry(
			const ClusterData::ConstAtomView& atom,
			const geometry::Symm<FloatingPointType>& symmetry,
			int symm_index) const {

			auto temp_point = symmetry.GenSymm(atom.point);

			PointType shifted_point;
			ShiftType shift;
			for (int i = 0; i < 3; ++i) {
				FloatingPointType coord = temp_point[i];
				auto s = static_cast<int8_t>(std::floor(coord));
				shifted_point[i] = coord - s;
				shift[i] = s;
			}

			return {atom.index, atom.type, shifted_point, symm_index, shift};
		}

		bool is_duplicate(const TempAtom& candidate, const ClusterData& existing, FloatingPointType tolerance) const {
			for (size_t i = 0; i < existing.size(); ++i) {
				if (geometry::Point<FloatingPointType>::distance(candidate.point, existing[i].point) < tolerance) {
					return true;
				}
			}
			return false;
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
			: cell(unit_cell), dist(distances) {}

		ResultType execute(const ClusterData& atoms_01) const {
			ResultType result;

			auto atom_s = atoms_01.size();


			geometry::SpatialGrid<FloatingPointType> sg;
			sg.build(atoms_01.points, cell, 4.0);
			auto bonds = sg.get_bonds(false);
			dist.filter_bond_list(bonds,
								  atoms_01.types,
								  atoms_01.points,
								  [this](const PointType& a, const PointType& b)
								  {
									  return (cell.fracToCart() * (a - b)).r();
								  });

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
			const auto bondshift = geometry::SpatialGrid<FloatingPointType>::decompress_shift(bond.shiftcode);

			// Check mol_a and mol_b are the same molecules
			if (&mol_a == &mol_b) {
				// Is it polymer?
				auto totalshift = find(bond.second, mol_a).shift - find(bond.first, mol_a).shift - bondshift;
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
			const auto dshift = fshift + bondshift;
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
		using ClusterData = cluster_detail::ClusterData;
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
		ClusterData execute(const Distances& distances, bool& hasPolymer) const
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
			if (unit_cell.success() == false)
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
					hasPolymer = true;
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
			ClusterData ret;
			ret.reserve(atoms.size());

			for (const auto& atom : atoms) {
				ret.push_back(unit_cell.atoms[atom.id].index,
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



			std::array<FloatingPointType, 3> dp = {plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111)};

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
				plane[2].distance(cell.fracToCart() * anchor.point) - b[2] * dp[2]};

			std::array < FloatingPointType, 3> high{
				dp[0] - low[0],
				dp[1] - low[1],
				dp[2] - low[2]};

			// [ -x, +x, -y, +y, -z, +z ]
			const std::array<ShiftType::value_type, 6> maxr{
				b[0] - static_cast<ShiftType::value_type>(std::ceil((cutoff - low[0]) / dp[0])),
				b[0] + static_cast<ShiftType::value_type>(std::ceil((cutoff - high[0]) / dp[0])),
				b[1] - static_cast<ShiftType::value_type>(std::ceil((cutoff - low[1]) / dp[1])),
				b[1] + static_cast<ShiftType::value_type>(std::ceil((cutoff - high[1]) / dp[1])),
				b[2] - static_cast<ShiftType::value_type>(std::ceil((cutoff - low[2]) / dp[2])),
				b[2] + static_cast<ShiftType::value_type>(std::ceil((cutoff - high[2]) / dp[2]))};


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
							const ClusterData& unit01) const {
			for (const auto& anchor : anchors_frac) {
				for (const auto& node : mol.nodes) {
					// Calculate distance to anchor
					PointType vec = cell.fracToCart() * (unit01.points[node.id] + node.shift + shift - anchor.point);

					if (vec.r() < anchor.radius) {
						return true;
					}
				}
			}
			return false;
		}


		std::unordered_set<TranslatedAtom, TranslatedAtom::Hash>
			grow_polymer(const Molecule& molecule,
						 const ClusterData& unit01,
						 const std::array<Plane, 3>& plane) const
		{
			std::unordered_set<TranslatedAtom, TranslatedAtom::Hash> ret;

			BoxSet boxes;
			auto rp111 = cell.fracToCart() * ShiftType(1, 1, 1);

			std::array<FloatingPointType, 3> dp = {plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111)};

			for (const auto& anchor : anchors_frac) {
				constructBox(anchor, polymer_cutoff_radius, dp, boxes, plane);
			}

			for (const auto& shift : boxes) {
				for (const auto& node : molecule.nodes) {
					for (const auto& anchor : anchors_frac) {
						PointType vec = unit01.points[node.id] + shift - anchor.point;

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