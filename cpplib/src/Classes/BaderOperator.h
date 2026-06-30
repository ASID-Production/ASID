// Copyright 2026 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
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
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <utility>
#include <vector>

#include <string>
#include <format>


#include "../Classes/Reader.h"
#include "../Classes/Splines.h"

namespace cpplib {

	// Global vector of RadialSplines of Electron Density
	inline const auto ElectronDensitySplines = DensityParser::parse_file(".\\ED.txt");

	template <typename T, size_t N = 2>
	struct PointsSoA {

		// Format: [min_x, min_y, min_z, max_x, max_y, max_z]
		template <typename T2>
		using BoundsArray = std::array<T2, 6>;

		std::vector<T> x;
		std::vector<T> y;
		std::vector<T> z;
		std::vector<uint32_t> voron_ids;
		std::vector<uint32_t> spline_ids;

		std::vector<uint32_t> ix;
		std::vector<uint32_t> iy;
		std::vector<uint32_t> iz;
		std::vector<uint32_t> ids_sh;
		std::vector<uint32_t> offset_sh;

		BoundsArray<T> cartesians;
		std::array<uint32_t, 3> grid_dim;
		T one_over_period = T(1);


		using PointType = geometry::Point<T>;

		static constexpr size_t NEAR = N;
		static constexpr size_t DIM = NEAR * 2 + 1;
		static constexpr size_t TOTAL_SHIFTS = DIM * DIM * DIM;

		template<typename I>
		static consteval auto unrollPositions() {
			using LocalPointType = typename geometry::Point<I>;
			std::array<LocalPointType, TOTAL_SHIFTS> result;
			constexpr size_t DIM2 = DIM * DIM;
			for (int i = 0; i < DIM; i++) {
				const size_t idim2 = i * DIM2;
				for (int j = 0; j < DIM; j++) {
					const size_t jdim = j * DIM;
					for (int k = 0; k < DIM; k++) {
						result[idim2 + jdim + k] =
							LocalPointType(i - static_cast<int>(NEAR),
										   j - static_cast<int>(NEAR),
										   k - static_cast<int>(NEAR));
					}
				}
			}
			return result;
		}


		static constexpr auto p_near = unrollPositions<char>();
		std::array<int32_t, TOTAL_SHIFTS> flat_shifts;


		void computeRuntimeFlatShifts() {

			const int32_t s_x = 1;
			const int32_t s_y = static_cast<int32_t>(grid_dim[0]);
			const int32_t s_z = s_y * static_cast<int32_t>(grid_dim[1]);

			for (size_t i = 0; i < TOTAL_SHIFTS; ++i) {
				int32_t dx = static_cast<int32_t>(p_near[i][0]);
				int32_t dy = static_cast<int32_t>(p_near[i][1]);
				int32_t dz = static_cast<int32_t>(p_near[i][2]);

				flat_shifts[i] = dz * s_z + dy * s_y + dx * s_x;
			}
		}


		PointsSoA() = default;
		void reserve(size_t capacity) {
			x.reserve(capacity);
			y.reserve(capacity);
			z.reserve(capacity);
			voron_ids.reserve(capacity);
			spline_ids.reserve(capacity);

			ix.reserve(capacity);
			iy.reserve(capacity);
			iz.reserve(capacity);
			ids_sh.reserve(capacity);
		}
		void addPoint(PointType point, uint32_t voron_id, uint32_t spine_id) {
			x.push_back(point[0]);
			y.push_back(point[1]);
			z.push_back(point[2]);

			voron_ids.push_back(voron_id);
			spline_ids.push_back(spine_id);
		}
		void calculateSpatialIndexes(T one_over_period_in) {
			const uint32_t size = static_cast<uint32_t>(x.size());
			one_over_period = one_over_period_in;

			// 1. Cleanup and reservation
			ix.clear(); ix.reserve(size);
			iy.clear(); iy.reserve(size);
			iz.clear(); iz.reserve(size);

			const int32_t min_fx = static_cast<int32_t>(std::floor(cartesians[0] * one_over_period));
			const int32_t min_fy = static_cast<int32_t>(std::floor(cartesians[1] * one_over_period));
			const int32_t min_fz = static_cast<int32_t>(std::floor(cartesians[2] * one_over_period));

			// Fill dim sizes
			grid_dim[0] = static_cast<uint32_t>(static_cast<int32_t>(std::ceil(cartesians[3] * one_over_period)) - min_fx + 1);
			grid_dim[1] = static_cast<uint32_t>(static_cast<int32_t>(std::ceil(cartesians[4] * one_over_period)) - min_fy + 1);
			grid_dim[2] = static_cast<uint32_t>(static_cast<int32_t>(std::ceil(cartesians[5] * one_over_period)) - min_fz + 1);

			const uint32_t grid_x = grid_dim[0];
			const uint32_t grid_xy = grid_x * grid_dim[1];

			uint32_t total_cells = grid_dim[0] * grid_dim[1] * grid_dim[2];

			// Fill offset_sh[] with zeros. Size = total_cells + 1
			offset_sh.assign(total_cells + 1, 0);

			std::vector<uint32_t> particle_cell_ids;
			particle_cell_ids.reserve(size);


			// 2. Counting particles in boxes
			for (uint32_t i = 0; i < size; ++i) {
				int32_t fx = static_cast<int32_t>(std::floor(x[i] * one_over_period));
				int32_t fy = static_cast<int32_t>(std::floor(y[i] * one_over_period));
				int32_t fz = static_cast<int32_t>(std::floor(z[i] * one_over_period));

				uint32_t local_ix = static_cast<uint32_t>(std::clamp(fx - min_fx, 0, static_cast<int32_t>(grid_dim[0]) - 1));
				uint32_t local_iy = static_cast<uint32_t>(std::clamp(fy - min_fy, 0, static_cast<int32_t>(grid_dim[1]) - 1));
				uint32_t local_iz = static_cast<uint32_t>(std::clamp(fz - min_fz, 0, static_cast<int32_t>(grid_dim[2]) - 1));

				ix.push_back(local_ix);
				iy.push_back(local_iy);
				iz.push_back(local_iz);

				// Calculate flat index (3D -> 1D)
				uint32_t cell_id = local_iz * grid_xy + local_iy * grid_x + local_ix;
				particle_cell_ids.push_back(cell_id);

				offset_sh[cell_id + 1]++;
			}

			// 3. Finalising offset array and creating a copy
			for (uint32_t c = 1; c <= total_cells; ++c) {
				offset_sh[c] += offset_sh[c - 1];
			}
			std::vector<uint32_t> current_offsets = offset_sh;

			// 4. Filling ids_sh
			ids_sh.resize(size);
			for (uint32_t i = 0; i < size; ++i) {
				uint32_t cell_id = particle_cell_ids[i];
				uint32_t dest_idx = current_offsets[cell_id];
				ids_sh[dest_idx] = i;
				current_offsets[cell_id]++;
			}
		}

		template<typename AtomType>
		void initializeWithPeriodicImages(const std::vector<geometry::Point<T>>& base_fractional_points,
										  const std::vector<AtomType>& base_spline_ids,
										  const geometry::Matrix<T>& mat_FracToCart,
										  const geometry::Matrix<T>& mat_CartToFrac,
										  T cutoff) {
			BoundsArray<int32_t> frac_bounds = getBoundsFrac(mat_CartToFrac, cutoff);
			getBoundsCart(mat_FracToCart, frac_bounds);

			int32_t min_i = frac_bounds[0];
			int32_t min_j = frac_bounds[1];
			int32_t min_k = frac_bounds[2];
			int32_t max_i = frac_bounds[3];
			int32_t max_j = frac_bounds[4];
			int32_t max_k = frac_bounds[5];

			size_t num_images = (max_i - min_i) * (max_j - min_j) * (max_k - min_k);
			size_t base_points_count = base_fractional_points.size();
			size_t total_estimated_points = base_points_count * num_images;

			this->reserve(total_estimated_points);
			x.clear(); y.clear(); z.clear();
			voron_ids.clear(); spline_ids.clear();

			for (int32_t i = min_i; i < max_i; ++i) {
				for (int32_t j = min_j; j < max_j; ++j) {
					for (int32_t k = min_k; k < max_k; ++k) {

						T shift_cart_x = mat_FracToCart.El(0, 0) * i + mat_FracToCart.El(0, 1) * j + mat_FracToCart.El(0, 2) * k;
						T shift_cart_y = mat_FracToCart.El(1, 0) * i + mat_FracToCart.El(1, 1) * j + mat_FracToCart.El(1, 2) * k;
						T shift_cart_z = mat_FracToCart.El(2, 0) * i + mat_FracToCart.El(2, 1) * j + mat_FracToCart.El(2, 2) * k;

						for (size_t pt_idx = 0; pt_idx < base_points_count; ++pt_idx) {
							const auto& f_point = base_fractional_points[pt_idx];

							T base_cart_x = mat_FracToCart.El(0, 0) * f_point[0] + mat_FracToCart.El(0, 1) * f_point[1] + mat_FracToCart.El(0, 2) * f_point[2];
							T base_cart_y = mat_FracToCart.El(1, 0) * f_point[0] + mat_FracToCart.El(1, 1) * f_point[1] + mat_FracToCart.El(1, 2) * f_point[2];
							T base_cart_z = mat_FracToCart.El(2, 0) * f_point[0] + mat_FracToCart.El(2, 1) * f_point[1] + mat_FracToCart.El(2, 2) * f_point[2];

							x.push_back(base_cart_x + shift_cart_x);
							y.push_back(base_cart_y + shift_cart_y);
							z.push_back(base_cart_z + shift_cart_z);

							spline_ids.push_back(static_cast<uint32_t>(base_spline_ids[pt_idx]));
							voron_ids.push_back(pt_idx);
						}
					}
				}
			}

			calculateSpatialIndexes(one_over_period);
			computeRuntimeFlatShifts();
		}


		constexpr static BoundsArray<int32_t> getBoundsFrac(const geometry::Matrix<T> mat_CartToFrac, T cutoff) {
			BoundsArray<int32_t> indexes;
			for (int i = 0; i < 3; ++i) {
				auto comp_x = mat_CartToFrac.El(i, 0);
				auto comp_y = mat_CartToFrac.El(i, 1);
				auto comp_z = mat_CartToFrac.El(i, 2);

				auto row_length = std::sqrt(comp_x * comp_x + comp_y * comp_y + comp_z * comp_z);
				T delta = cutoff * row_length;

				indexes[i] = static_cast<int32_t>(std::floor(-delta));
				indexes[i + 3] = static_cast<int32_t>(std::ceil(static_cast<T>(1.0) + delta));
			}
			return indexes;
		}

		constexpr void getBoundsCart(const geometry::Matrix<T> mat_FracToCart, const BoundsArray<int32_t> indexes) {
			for (int i = 0; i < 3; ++i) {
				T cart_min = 0;
				T cart_max = 0;

				for (int j = 0; j < 3; ++j) {
					T elem = mat_FracToCart.El(i, j);

					T v1 = elem * indexes[j];
					T v2 = elem * indexes[j + 3];

					cart_min += std::min(v1, v2);
					cart_max += std::max(v1, v2);
				}

				cartesians[i] = cart_min;
				cartesians[i + 3] = cart_max;
			}
		}

		std::array<int32_t, TOTAL_SHIFTS> generateFlatOffsets() const {

			std::array<int32_t, TOTAL_SHIFTS>  flatOffsets;

			const int32_t s_y = static_cast<int32_t>(grid_dim[0]);
			const int32_t s_z = s_y * static_cast<int32_t>(grid_dim[1]);

			size_t i = 0;

			constexpr int32_t pos_n = static_cast<int32_t>(NEAR);
			constexpr int32_t neg_n = -static_cast<int32_t>(NEAR);

			for (int32_t dz = neg_n; dz <= pos_n; ++dz) {
				int32_t offset_z = dz * s_z;

				for (int32_t dy = neg_n; dy <= pos_n; ++dy) {
					int32_t offset_zy = offset_z + (dy * s_y);

					for (int32_t dx = neg_n; dx <= pos_n; ++dx) {
						int32_t final_offset = offset_zy + dx;
						flatOffsets[i] = final_offset;
						i++;
					}
				}
			}

			return flatOffsets;
		}

		uint32_t getCellIndexForPoint(T point_x, T point_y, T point_z) const {
			const int32_t min_fx = static_cast<int32_t>(std::floor(cartesians[0] * one_over_period));
			const int32_t min_fy = static_cast<int32_t>(std::floor(cartesians[1] * one_over_period));
			const int32_t min_fz = static_cast<int32_t>(std::floor(cartesians[2] * one_over_period));

			int32_t fx = static_cast<int32_t>(std::floor(point_x * one_over_period));
			int32_t fy = static_cast<int32_t>(std::floor(point_y * one_over_period));
			int32_t fz = static_cast<int32_t>(std::floor(point_z * one_over_period));

			uint32_t local_ix = static_cast<uint32_t>(std::clamp(fx - min_fx, 0, static_cast<int32_t>(grid_dim[0]) - 1));
			uint32_t local_iy = static_cast<uint32_t>(std::clamp(fy - min_fy, 0, static_cast<int32_t>(grid_dim[1]) - 1));
			uint32_t local_iz = static_cast<uint32_t>(std::clamp(fz - min_fz, 0, static_cast<int32_t>(grid_dim[2]) - 1));

			const uint32_t grid_x = grid_dim[0];
			const uint32_t grid_xy = grid_x * grid_dim[1];

			return local_iz * grid_xy + local_iy * grid_x + local_ix;
		}
		inline uint32_t getCellIndexForPoint(const PointType& point) const {
			return getCellIndexForPoint(point[0], point[1], point[2]);
		}


	};

	class CriticalPoint {
	public:
		using PointType = typename RadialSpline::PointType;
		using value_type = typename PointType::value_type;
		static constexpr value_type MAX_STEP = 0.1;
		static constexpr value_type MAX_STEP_SQ = MAX_STEP * MAX_STEP;
		static constexpr value_type EPS = 1e-12;

		enum class TYPE {
			N = 0,
			B = 1,
			R = 2,
			C = 3
		};
		CriticalPoint() = default;
		CriticalPoint(TYPE type, const PointType& pos) : type_(type), pos_(pos) {
		}
		inline PointType FindNextPosition() const {
			switch (type_) {
				using enum TYPE;
				case N:
					return PartitionedRFO<0>();
				case B:
					return PartitionedRFO<1>();
				case R:
					return PartitionedRFO<2>();
				case C:
					return PartitionedRFO<3>();
			}
			return pos_;
		}
		TripleDouble CalculateEDinPoint(const PointsSoA<value_type>& psoa,
										const std::vector<RadialSpline>& splines) {
			TripleDouble result;

			const uint32_t home_cell = psoa.getCellIndexForPoint(pos_);
			const uint32_t total_cells = static_cast<uint32_t>(psoa.offset_sh.size() - 1);

			for (int32_t shift:psoa.flat_shifts) {
				int64_t cell_idx = static_cast<int64_t>(home_cell) + shift;
				if (cell_idx < 0 || cell_idx >= static_cast<int64_t>(total_cells))
					continue;

				uint32_t start = psoa.offset_sh[cell_idx];
				uint32_t end = psoa.offset_sh[cell_idx + 1];

				for (uint32_t idx = start; idx < end; ++idx) {
					uint32_t atom_i = psoa.ids_sh[idx];

					PointType delta(pos_[0] - psoa.x[atom_i],
									pos_[1] - psoa.y[atom_i],
									pos_[2] - psoa.z[atom_i]);

					uint32_t spline_id = psoa.spline_ids[atom_i];
					result += splines[spline_id].evaluate(delta);
				}
			}

			return result;
		}



		constexpr void UpdatePos(const PointType& pos) noexcept {
			pos_ = pos;
		}
		constexpr void UpdateVal(const TripleDouble& value) noexcept {
			value_ = value;
		}
	private:
		PointType EigenVectorFollowing() const {
			const auto& g = value_.grad;
			const auto& h = value_.hess;

			auto eig = h.EigenvaluesAndVectors();

			int idx = 0;
			if (type_ == TYPE::B) {
				value_type max_val = eig.values[0];
				if (eig.values[1] > max_val) {
					max_val = eig.values[1];
					idx = 1;
				}
				if (eig.values[2] > max_val) {
					max_val = eig.values[2];
					idx = 2;
				}
			} else {
				value_type min_val = eig.values[0];
				if (eig.values[1] < min_val) {
					min_val = eig.values[1];
					idx = 1;
				}
				if (eig.values[2] < min_val) {
					min_val = eig.values[2];
					idx = 2;
				}
			}

			const auto& v = eig.vectors[idx];
			value_type lambda = eig.values[idx];

			value_type grad_dot_v = g[0] * v[0] + g[1] * v[1] + g[2] * v[2];

			value_type alpha = 0.0;

			if (std::abs(lambda) > EPS) {
				alpha = -grad_dot_v / lambda;
			} else {
				alpha = -grad_dot_v * 0.01;
			}

			alpha = std::clamp(alpha, -MAX_STEP, MAX_STEP);

			return pos_ + PointType(v[0] * alpha, v[1] * alpha, v[2] * alpha);
		}


		PointType NewtonRaphsonPredict() const {
			const auto& g = value_.grad;
			const auto& h = value_.hess;

			value_type det = h.Det();
			if (std::abs(det) < EPS) {
				value_type step = 0.01;
				return pos_ + PointType(-g[0] * step, -g[1] * step, -g[2] * step);
			}

			const auto H_inv = h.Invert();
			value_type delta_x = -(H_inv.El(0, 0) * g[0] + H_inv.El(0, 1) * g[1] + H_inv.El(0, 2) * g[2]);
			value_type delta_y = -(H_inv.El(1, 0) * g[0] + H_inv.El(1, 1) * g[1] + H_inv.El(1, 2) * g[2]);
			value_type delta_z = -(H_inv.El(2, 0) * g[0] + H_inv.El(2, 1) * g[1] + H_inv.El(2, 2) * g[2]);

			const value_type alpha_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
			if (alpha_sq > MAX_STEP_SQ) {
				const value_type step_over_alpha = std::sqrt(MAX_STEP_SQ / alpha_sq);
				delta_x *= step_over_alpha;
				delta_y *= step_over_alpha;
				delta_z *= step_over_alpha;
			}

			return pos_ + PointType(delta_x, delta_y, delta_z);
		}

		template<int Target>
		PointType PartitionedRFO() const {
			static_assert(Target >= 0 && Target <= 3,
						  "Target index must be 0 (min), 1 (saddle), 2 (ring), or 3 (max)");

			const auto& g = value_.grad;
			auto A = value_.hess; // a copy
			auto eig = A.EigenvaluesAndVectors();  // values[0] >= values[1] >= values[2]

			const std::array<value_type, 3> lam = {eig.values[2], eig.values[1], eig.values[0]};

			value_type eta;
			if constexpr (Target == 0) {
				eta = lam[0] - std::max(1e-4, std::abs(lam[0]) * 1e-3);
			} else if constexpr (Target == 3) {
				eta = lam[2] + std::max(1e-4, std::abs(lam[2]) * 1e-3);
			} else {
				// Saddle Target (1 or 2)
				value_type left = lam[Target - 1];
				value_type right = lam[Target];
				if (right - left < 1e-8) [[unlikely]] {
					right = left + 1e-6;
					left = left - 1e-6;
				}
				eta = 0.5 * (left + right);
			}

			// A = H - eta * I
			A.El(0, 0) -= eta;
			A.El(1, 1) -= eta;
			A.El(2, 2) -= eta;

			const value_type rhs[3] = {-g[0], -g[1], -g[2]};

			const value_type e00 = A.El(0, 0);
			const value_type e01 = A.El(0, 1);
			const value_type e02 = A.El(0, 2);
			const value_type e10 = A.El(1, 0);
			const value_type e11 = A.El(1, 1);
			const value_type e12 = A.El(1, 2);
			const value_type e20 = A.El(2, 0);
			const value_type e21 = A.El(2, 1);
			const value_type e22 = A.El(2, 2);

			const value_type c00 = e11 * e22 - e12 * e12;
			const value_type c11 = e00 * e22 - e02 * e02;
			const value_type c22 = e00 * e11 - e01 * e01;
			const value_type c01 = e12 * e02 - e01 * e22;
			const value_type c02 = e01 * e12 - e02 * e11;
			const value_type c12 = e02 * e12 - e00 * e12;

			const value_type detA = e00 * (e11 * e22 - e12 * e21)
				                  - e01 * (e10 * e22 - e12 * e20)
				                  + e02 * (e10 * e21 - e11 * e20);

			if (std::abs(detA) < EPS) [[unlikely]] {
				const value_type step = 0.1;
				return pos_ + PointType(-g[0] * step, -g[1] * step, -g[2] * step);
			}

			const value_type invDet = 1.0 / detA;
			value_type dx[3];
			dx[0] = (c00 * rhs[0] + c01 * rhs[1] + c02 * rhs[2]) * invDet;
			dx[1] = (c01 * rhs[0] + c11 * rhs[1] + c12 * rhs[2]) * invDet;
			dx[2] = (c02 * rhs[0] + c12 * rhs[1] + c22 * rhs[2]) * invDet;

			const value_type len_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
			if (len_sq > MAX_STEP_SQ) [[unlikely]] {
				const value_type scale = std::sqrt(MAX_STEP_SQ / len_sq);
				dx[0] *= scale;
				dx[1] *= scale;
				dx[2] *= scale;
			}

			return pos_ + PointType(dx[0], dx[1], dx[2]);
		}

	public:
		TYPE type_ = TYPE::N;
		PointType pos_{};
		TripleDouble value_{};
	};

	class BaderOperator {
	public:
		using PointType = typename geometry::Point<TripleDouble::value_type>;
		using value_type = CubicSpline<>::value_type;
		static constexpr value_type EPS = voronoi::EPSILON;

		static value_type GetOptimalRadius(value_type radius,
										   value_type vertical_eps,
										   const std::vector<RadialSpline>& splines,
										   const std::vector<basic_types::AtomTypeBase>& unfiltered_types) {
			value_type max_x = 0;

			std::vector<bool> visited(splines.size(), false);

			for (const auto& active_index : unfiltered_types) {
				const size_t idx = static_cast<size_t>(active_index);

				if (visited[idx]) {
					continue;
				}
				visited[idx] = true;

				value_type phi, dphi, ddphi;
				auto& active_spline = splines[idx].spline;

				active_spline.eval(radius, phi, dphi, ddphi);
				value_type x = active_spline.find_value(radius, phi * vertical_eps);

				max_x = std::max(x, max_x);
			}

			return max_x;
		}

		static std::vector<CriticalPoint> GenerateInitialCriticalPoints(const voronoi::VoronoiFused& vf, const geometry::Matrix<value_type>& fracToCart) {

			std::vector<CriticalPoint> result;
			result.reserve(vf.vertices.size() + vf.edges.size() + vf.polygons.size());

			auto add_candidate = [&](const PointType& p, CriticalPoint::TYPE t) {
				value_type x = p.a[0] - std::floor(p.a[0]);
				value_type y = p.a[1] - std::floor(p.a[1]);
				value_type z = p.a[2] - std::floor(p.a[2]);

				if (x >= 1.0 - EPS) x = 0.0;
				if (y >= 1.0 - EPS) y = 0.0;
				if (z >= 1.0 - EPS) z = 0.0;
				result.emplace_back(t, fracToCart * CriticalPoint::PointType(x, y, z));
				};

			// 1. Type C (Cage)
			for (const PointType& p : vf.vertices) {
				add_candidate(p, CriticalPoint::TYPE::C);
			}

			// 2. Type R (Ring)
			for (const auto& edge : vf.edges) {
				const auto& v1 = vf.vertices[edge.vert_ids[0]];
				const auto& v2 = vf.vertices[edge.vert_ids[1]];
				PointType mid((v1[0] + v2[0]) * 0.5,
							  (v1[1] + v2[1]) * 0.5,
							  (v1[2] + v2[2]) * 0.5);
				add_candidate(mid, CriticalPoint::TYPE::R);
			}

			// 3. Type B (Bond)
			for (const auto& poly : vf.polygons) {
				PointType center(0.0, 0.0, 0.0);
				for (uint32_t vid:poly.vert_ids) {
					center += vf.vertices[vid];
				}
				center /= static_cast<value_type>(poly.vert_ids.size());
				add_candidate(center, CriticalPoint::TYPE::B);
			}



			return result;
		}

		static void OptimizePoint(CriticalPoint& cp, const PointsSoA<value_type>& psoa, value_type end_sq) {

			TripleDouble cur_vals = cp.CalculateEDinPoint(psoa, ElectronDensitySplines);
			cp.UpdateVal(cur_vals);
			PointType next_point = cp.FindNextPosition();
			value_type drs = (next_point - cp.pos_).rSq();
			value_type ddrs = 10.0;
			value_type min_ddrs = end_sq * end_sq;
			int counter = 0;
			value_type adaptive = 0.001;
			while ((drs > end_sq && ddrs > min_ddrs) || drs >= CriticalPoint::MAX_STEP_SQ*adaptive) {
				cp.UpdatePos(next_point);
				cur_vals = cp.CalculateEDinPoint(psoa, ElectronDensitySplines);
				cp.UpdateVal(cur_vals);
				value_type old_drs = drs;
				next_point = cp.FindNextPosition();
				drs = (next_point - cp.pos_).rSq();
				ddrs = old_drs - drs;
				counter++;
				adaptive *= 1.2;
			}
			cp.UpdatePos(next_point);
			cur_vals = cp.CalculateEDinPoint(psoa, ElectronDensitySplines);
			cp.UpdateVal(cur_vals);
		}

		static std::string criticalPointsToPDB_fast(const std::vector<CriticalPoint>&cps) {
			static constexpr std::array<const char*, 4> typeToElem = {
				"C", // N
				"N", // B
				"O", // R
				"F"  // C
			};
			static constexpr size_t LINE_LENGTH = 80;

			static constexpr std::string_view HEADER = "REMARK Generated by critical point converter\n";
			static constexpr std::string_view FOOTER = "END\n";

			std::string result;
			result.reserve(HEADER.size() + cps.size() * LINE_LENGTH + FOOTER.size());

			result += HEADER;

			for (size_t i = 0; i < cps.size(); ++i) {
				const auto& cp = cps[i];
				const int serial = static_cast<int>(i + 1);

				const char* elem = typeToElem[static_cast<int>(cp.type_)];

				const double x = cp.pos_[0];
				const double y = cp.pos_[1];
				const double z = cp.pos_[2];

				std::format_to(
					std::back_inserter(result),
					"ATOM  {:>5}  {:<4} CP  {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00{:>12}\n",
					serial, elem, serial,
					x, y, z,
					elem
				);
			}

			result += FOOTER;
			return result;
		}
		static void writeCriticalPointsToFile(const std::vector<CriticalPoint>& cps, const char* filename) {
			const std::string data = criticalPointsToPDB_fast(cps);

			std::ofstream file(filename, std::ios::binary);

			//constexpr size_t bufSize = 1 << 20; 
			//std::unique_ptr<char[]> buf(new char[bufSize]);
			//file.rdbuf()->pubsetbuf(buf.get(), bufSize);

			file.write(data.data(), data.size());
		}
	};
}
