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
#include <vector>
#include <tuple>
#include <algorithm>

#include "../BaseHeaders/Currents.h"
#include "Geometry.h"
#include "DSU.h"

namespace cpplib {

	class Cluster {
	private:
		using FAMS = FAM_Struct;
		using FAMC = FAM_Cell;
	public:
		using ShiftType = FAMS::ShiftType;
		using BoxValueType = std::pair<ShiftType, std::vector<bool>>;
		using BoxType = std::map<ShiftType, std::vector<bool>>;
		using BoxValueTypeConst = BoxType::value_type;
		using PointType = FAMS::PointType;
		using FloatingPointType = PointType::value_type;
		using FracToCartMatrixType = cpplib::geometry::Matrix<FloatingPointType>;
		using Plane = cpplib::geometry::Plane<FloatingPointType>;
		using AnchorType = std::pair<PointType, FloatingPointType>;

	private:
		FAMS& fs;
		const FAMC& fc;
		const std::vector<AnchorType>& anchors;
		std::array<PointType, 3> e;
		std::array<Plane, 3> plane;
		std::vector<AnchorType::first_type> r_anchors;

	public:
		Cluster() = delete;
		Cluster(FAMS& f_str, const FAMC& f_cell, const std::vector<AnchorType>& anch) noexcept : fs(f_str), fc(f_cell), anchors(anch) {
			e = { fc.fracToCart() * ShiftType(1,0,0), fc.fracToCart() * ShiftType(0,1,0), fc.fracToCart() * ShiftType(0,0,1) };
			constexpr PointType zeroPoint(0, 0, 0);
			plane = { Plane(zeroPoint, e[1], e[2]), Plane(zeroPoint, e[0], e[2]), Plane(zeroPoint, e[0], e[1]) };
			r_anchors.reserve(anchors.size());
			for (const auto& a : anchors) {
				r_anchors.emplace_back(fc.fracToCart() * a.first);
			}
		};

		void CreateBox(BoxType& box, FloatingPointType over_radius = 0.0) {
			auto rp111 = fc.fracToCart() * ShiftType(1, 1, 1);
			std::array<FloatingPointType, 3> dp = { plane[0].distance(rp111),
													plane[1].distance(rp111),
													plane[2].distance(rp111) };

			for (const auto& pair : anchors)
			{
				ShiftType b(pair.first.floor());
				std::array < FloatingPointType, 3> low{
					plane[0].distance(fc.fracToCart() * pair.first) - b[0] * dp[0],
					plane[1].distance(fc.fracToCart() * pair.first) - b[1] * dp[1],
					plane[2].distance(fc.fracToCart() * pair.first) - b[2] * dp[2] };

				std::array < FloatingPointType, 3> high{
					dp[0] - low[0],
					dp[1] - low[1],
					dp[2] - low[2] };

				auto radius = over_radius < 0.01 ? pair.second : over_radius;

				// [ -x, +x, -y, +y, -z, +z ]
				const std::array<ShiftType::value_type, 6> maxr{
					b[0] - static_cast<ShiftType::value_type>(ceil((radius - low[0]) / dp[0])),
					b[0] + static_cast<ShiftType::value_type>(ceil((radius - high[0]) / dp[0])),
					b[1] - static_cast<ShiftType::value_type>(ceil((radius - low[1]) / dp[1])),
					b[1] + static_cast<ShiftType::value_type>(ceil((radius - high[1]) / dp[1])),
					b[2] - static_cast<ShiftType::value_type>(ceil((radius - low[2]) / dp[2])),
					b[2] + static_cast<ShiftType::value_type>(ceil((radius - high[2]) / dp[2])), };
				constructBox(maxr, box);
			}
		}


		void Grow(BoxType& basebox,
						 BoxType& newbox,
						 const std::vector<std::vector<ShiftType>>& molecule,
						 FAMS::AtomIndex ai,
						 const ShiftType& curshift) {
			auto baseshift = curshift - molecule[ai][0];
			for (FAMS::AtomIndex i = 0; i < molecule.size(); i++)
			{
				if (molecule[i].empty())
					continue;
				auto newshift = baseshift + molecule[i][0];
				auto it = basebox.find(newshift);
				if (it == basebox.end()) { // did not find
					auto newit = newbox.insert(std::make_pair(newshift, std::vector<bool>(molecule.size()))).first;
					newit->second[i] = true;
				}
				else { // found
					it->second[i] = true;
				}
			}
		}

		void GrowPoly(BoxType& box,
							 const std::vector<std::pair<std::vector<std::vector<ShiftType>>, bool>>& molecules,
							 const std::vector<bool>& polyflags,
							 const FloatingPointType over_radius) {

			if (std::none_of(polyflags.begin(), polyflags.end(), [](bool a) {return a; })) return;

			CreateBox(box, over_radius);

			std::vector<bool> polyatoms(fs.sizePoints, false);
			for (size_t i = 1; i < polyflags.size(); i++)
			{
				if (polyflags[i] == false) continue;
				for (FAMS::AtomIndex j = 0; j < fs.sizePoints; j++)
				{
					polyatoms[j] = polyatoms[j] || (!(molecules[i].first[j].empty()));
				}
			}

			for (FAMS::AtomIndex i = 0; i < fs.sizePoints; i++) {
				if (polyatoms[i] == false)
					continue;
				polyFill(box, i, over_radius);
			}
		}
	private:
		inline void polyFill(BoxType& box, const FAMS::AtomIndex i, const FloatingPointType over_radius) {
			for (auto& pair : box) {
				auto r_point = fc.fracToCart() * (pair.first + fs.points[i]);
				for (const auto& ranch : r_anchors)
				{
					if ((ranch - r_point).r() <= over_radius) {
						pair.second[i] = true;
						break;
					}
				}
			}
		}
		void constructBox(const std::array<ShiftType::value_type, 6>& maxr, BoxType& box) const {
			for (char i = maxr[0]; i <= maxr[1]; i++) {
				for (char j = maxr[2]; j <= maxr[3]; j++) {
					for (char k = maxr[4]; k <= maxr[5]; k++) {
						box.emplace(ShiftType(i, j, k), std::vector<bool>(fs.sizePoints, false));
					}
				}
			}
		}
	};


	class Cluster_New {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;
		using ShiftType = geometry::Point<int8_t>;
		using Box = ::std::vector<bool>;
		using BoxMap = ::std::map<ShiftType, Box>;
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
			AtomTypeBase type;
			PointType point;
			SymmIndex symm;
			ShiftType shift;
		};
		struct TranslatedAtom{
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
			::std::vector<ClusterAtom> nodes{};
			//::std::vector<bool> atoms_in_molecule{};
			bool is_polymer = false;
		};

	private:
		CellType cell;
		::std::vector<AnchorType> anchors_frac;
		::std::vector<ClusterAtom> asymmetric_unit;
		::std::vector<SymmType> symm;
		FloatingPointType polymer_cutoff_radius;

	public:
		Cluster_New() = delete;

		void create(const DistancesType& distances) {
			auto unit_01 = construct_unit_01();
			auto molecules = constructMoleculesInUnit01(unit_01, distances);

		}

	private:
		::std::vector<Molecule> constructMoleculesInUnit01(const ::std::vector<ClusterAtom>& unit_01, const DistancesType& distances) const {
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

			// reconstruct molecules

			// Construct molecules from unit cell and bonds
			// auto moleculles = reconstructMolecules(unit_01, bonds);

			// move molecules to [0,1]

			// make Boxes

			// fill Boxes



			return{};
		}

		::std::vector<ClusterAtom> construct_unit_01() const {
			auto au_s = asymmetric_unit.size();
			size_t symm_s = symm.size();

			::std::vector<ClusterAtom> unit;
			unit.reserve(au_s * symm_s);

			for (SymmIndex i = 0; i < symm_s; i++) {
				for (size_t j = 0; j < au_s; j++) {
					auto temp_point = symm[i].GenSymm(asymmetric_unit[j].point);
					auto floating_shift = -temp_point.floor();
					temp_point.MoveToCell();
					if (isPointInVector(unit, temp_point)) continue;

					ShiftType shift(static_cast<ShiftType::value_type>(floating_shift[0]),
									static_cast<ShiftType::value_type>(floating_shift[1]),
									static_cast<ShiftType::value_type>(floating_shift[2]));
					unit.emplace_back(asymmetric_unit[j].type, temp_point, i, shift);
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
													 const BondList& bonds) const {
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

			::std::vector<Molecule> molecules(dsu.get_count_components());


			for (const auto& bond : filteredbonds) {
				
				
			}





			// TODO




			return {};
		}

		void mergeMolecules(std::vector<Molecule>& molecules_unique,
							std::vector<AtomIndex>& atom_ref_to_molecule,
							const AtomIndex a_ref,
							const AtomIndex b_ref) const {



			atom_ref_to_molecule.insert(atom_ref_to_molecule.end(),
										std::make_move_iterator(atom_ref_to_molecule.begin()),
										std::make_move_iterator(atom_ref_to_molecule.end()));
			// TODO

			return;
		}
	};
}