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
#include "../BaseHeaders/Currents.h"
#include "FindMolecules.h"
#include "Geometry.h"

namespace cpplib {
	class Cluster {
	private:
		using FAMS = currents::FAMStructType;
		using FAMC = currents::FAMCellType;
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
				auto b = FAMC::toShift(pair.first.floor());
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
		void constructBox(const std::array<ShiftType::value_type, 6>& maxr, BoxType& box) {
			for (char i = maxr[0]; i <= maxr[1]; i++) {
				for (char j = maxr[2]; j <= maxr[3]; j++) {
					for (char k = maxr[4]; k <= maxr[5]; k++) {
						box.emplace(ShiftType(i, j, k), std::vector<bool>(fs.sizePoints, false));
					}
				}
			}
		}
	};
}