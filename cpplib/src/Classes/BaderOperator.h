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
#include <iostream>


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

			calculateSpatialIndexes(NEAR/cutoff);
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
		static TripleDouble CalculateEDinPoint(PointType point, const PointsSoA<value_type>& psoa,
										const std::vector<RadialSpline>& splines) {
			TripleDouble result;

			const uint32_t home_cell = psoa.getCellIndexForPoint(point);
			const uint32_t total_cells = static_cast<uint32_t>(psoa.offset_sh.size() - 1);

			int total_splines = 0;
			int zero_splines = 0;

			for (int32_t shift:psoa.flat_shifts) {
				int64_t cell_idx = static_cast<int64_t>(home_cell) + shift;
				if (cell_idx < 0 || cell_idx >= static_cast<int64_t>(total_cells))
					continue;

				uint32_t start = psoa.offset_sh[cell_idx];
				uint32_t end = psoa.offset_sh[cell_idx + 1];

				for (uint32_t idx = start; idx < end; ++idx) {
					uint32_t atom_i = psoa.ids_sh[idx];

					PointType delta(point[0] - psoa.x[atom_i],
									point[1] - psoa.y[atom_i],
									point[2] - psoa.z[atom_i]);

					uint32_t spline_id = psoa.spline_ids[atom_i];
					auto TD = splines[spline_id].evaluate(delta);
					if (TD.val < 1.0e-10)
						zero_splines++;
					result += TD;
					total_splines++;
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
	public:
		PointType GradientPathInitVector(value_type step_length) const {
			auto eigen = value_.hess.EigenvaluesAndVectors();
			return eigen.vectors[0] * step_length;
		}
		template<value_type step_length>
		std::pair<std::vector<PointType>, std::vector<PointType>> CalculatePaths(
			const PointsSoA<value_type>& psoa,
			const std::vector<RadialSpline>& splines) const 
		{
			auto shift_pos = GradientPathInitVector(step_length);
			PointType next_pos = pos_ + shift_pos;
			PointType next_neg = pos_ - shift_pos;
			CriticalPoint cp_pos(TYPE::B, next_pos);
			CriticalPoint cp_neg(TYPE::B, next_neg);
			cp_pos.value_ = CalculateEDinPoint(next_pos, psoa, splines);
			cp_neg.value_ = CalculateEDinPoint(next_neg, psoa, splines);
			return std::make_pair(cp_pos.SinglePath<step_length>(psoa, splines),
								  cp_neg.SinglePath<step_length>(psoa, splines));
		}

		template<value_type step_length>
		std::vector<PointType> SinglePath(const PointsSoA<value_type>& psoa,
										  const std::vector<RadialSpline>& splines) {
			static_assert(step_length <= 5.0);
			constexpr size_t max_step = 5 / step_length - 1;

			std::vector<PointType> ret_val;
			ret_val.reserve(max_step + 1);
			ret_val.push_back(pos_);

			for (size_t i = 0; i < max_step; i++)
			{
				auto r = value_.grad.r();
				auto next_point = pos_ + value_.grad * (1 / r) * step_length;
				auto next_value = CalculateEDinPoint(next_point, psoa, splines);
				if (PointType::Scalar(value_.grad, next_value.grad) <= 0) {
					break;
				}
				pos_ = next_point;
				value_ = next_value;
				ret_val.push_back(pos_);
			}
			return ret_val;
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
		value_type solve_rfo_shift_trust(
			const std::vector<int>& indices,
			const std::array<value_type, 3>& lam,
			const std::array<value_type, 3>& f,
			bool maximize,
			value_type R) const
		{
			if (indices.empty()) return value_type(0);

			value_type lambda_bound = maximize?lam[indices[0]]:lam[indices[0]];
			for (int idx:indices) {
				if (maximize) lambda_bound = std::max(lambda_bound, lam[idx]);
				else          lambda_bound = std::min(lambda_bound, lam[idx]);
			}
			constexpr value_type tiny = 1e-12;
			value_type lo, hi;
			value_type f_max_abs = 0;
			for (int idx:indices) f_max_abs = std::max(f_max_abs, std::abs(f[idx]));

			if (maximize) {
				lo = lambda_bound + tiny;
				hi = lambda_bound + std::max(f_max_abs / R, value_type(1.0));
			} else {
				hi = lambda_bound - tiny;
				lo = lambda_bound - std::max(f_max_abs / R, value_type(1.0));
			}

			// Лямбда и f теперь правильно соотносятся через оригинальный idx
			auto Phi = [&](value_type eta) {
				value_type sum = -R * R;
				for (int idx:indices) {
					value_type diff = lam[idx] - eta;
					sum += (f[idx] * f[idx]) / (diff * diff);
				}
				return sum;
				};

			const int max_iter = 40;
			for (int iter = 0; iter < max_iter; ++iter) {
				value_type mid = (lo + hi) * 0.5;
				value_type phi_mid = Phi(mid);
				if (std::abs(phi_mid) < 1e-10 * R * R) return mid;
				if (maximize) {
					if (phi_mid > 0) lo = mid; else hi = mid;
				} else {
					if (phi_mid > 0) hi = mid; else lo = mid;
				}
				if (hi - lo < 1e-14) return mid;
			}
			return (lo + hi) * 0.5;
		}

		template<int Target>
		PointType PartitionedRFO() const {
			static_assert(Target >= 0 && Target <= 3,
						  "Target must be 0 (min), 1 (saddle1), 2 (saddle2), or 3 (max)");

			const auto& g = value_.grad;
			auto A = value_.hess;
			auto eig = A.EigenvaluesAndVectors();
			auto& lam = eig.values;
			auto& vecs = eig.vectors;

			// Проекции градиента на моды
			std::array<value_type, 3> f;
			for (int i = 0; i < 3; ++i) {
				f[i] = vecs[i][0] * g[0] + vecs[i][1] * g[1] + vecs[i][2] * g[2];
			}

			// Разделение на подпространства
			// Для седла 1 порядка (Target=1): первые 1 моды (минимальные) максимизируем, остальные минимизируем
			int max_modes = Target;
			std::vector<int> max_indices, min_indices;

			// Mode following (ведение моды)
			if (max_modes > 0) {
				if (!mode_following_init_) {
					for (int i = 0; i < max_modes; ++i) {
						max_indices.push_back(i);
						prev_max_vecs_[i] = vecs[i];
					}
					mode_following_init_ = true;
				} else {
					std::vector<bool> used_new(3, false);
					for (int j_old = 0; j_old < max_modes; ++j_old) {
						value_type best_overlap = -1.0;
						int best_i = -1;
						value_type best_sign = 1.0;

						for (int i_new = 0; i_new < 3; ++i_new) {
							if (used_new[i_new]) continue;
							value_type dot = vecs[i_new][0] * prev_max_vecs_[j_old][0]
								+ vecs[i_new][1] * prev_max_vecs_[j_old][1]
								+ vecs[i_new][2] * prev_max_vecs_[j_old][2];
							value_type ov = std::abs(dot);
							if (ov > best_overlap) {
								best_overlap = ov;
								best_i = i_new;
								best_sign = (dot < 0.0)?-1.0:1.0;
							}
						}
						max_indices.push_back(best_i);
						used_new[best_i] = true;

						// Стабилизируем знак вектора в истории, чтобы избежать осцилляций флипа
						prev_max_vecs_[j_old] = PointType(vecs[best_i][0] * best_sign,
														  vecs[best_i][1] * best_sign,
														  vecs[best_i][2] * best_sign);
					}
				}
			}

			for (int i = 0; i < 3; ++i) {
				if (std::find(max_indices.begin(), max_indices.end(), i) == max_indices.end())
					min_indices.push_back(i);
			}

			// Каноническое решение RFO уравнений (вместо QA-сдвига)
			// Функция ищет корень уравнения: Подпространство_RFO(eta) = 0
			auto solve_p_rfo = [&](const std::vector<int>& indices, bool is_max) -> value_type {
				if (indices.empty()) return 0.0;

				// Поиск границ для метода Ньютона/Бисекции
				value_type max_lam = -1e20, min_lam = 1e20;
				for (int idx:indices) {
					if (lam[idx] > max_lam) max_lam = lam[idx];
					if (lam[idx] < min_lam) min_lam = lam[idx];
				}

				// Функция RFO: \sum { f_i^2 / (eta - lam_i) } - eta = 0
				auto RFO_Eq = [&](value_type eta) {
					value_type sum = -eta;
					for (int idx:indices) {
						sum += (f[idx] * f[idx]) / (eta - lam[idx]);
					}
					return sum;
					};

				value_type lo, hi;
				if (is_max) {
					hi = min_lam - 1e-6;
					lo = min_lam - std::max(1.0, std::abs(min_lam) * 2.0);
				} else {
					lo = max_lam + 1e-6;
					hi = max_lam + std::max(1.0, std::abs(max_lam) * 2.0);
				}

				// Простой и надежный Bisection
				for (int iter = 0; iter < 60; ++iter) {
					value_type mid = 0.5 * (lo + hi);
					value_type val = RFO_Eq(mid);
					if (std::abs(val) < 1e-11 || (hi - lo) < 1e-12) return mid;

					if (val > 0) lo = mid;
					else         hi = mid;
				}
				return 0.5 * (lo + hi);
				};

			value_type eta_max = solve_p_rfo(max_indices, true);
			value_type eta_min = solve_p_rfo(min_indices, false);

			PointType dx(0, 0, 0);
			for (int i = 0; i < 3; ++i) {
				bool is_max = (std::find(max_indices.begin(), max_indices.end(), i) != max_indices.end());
				value_type eta = is_max?eta_max:eta_min;

				value_type denom = eta - lam[i];
				if (std::abs(denom) < 1e-10) denom = (denom >= 0)?1e-10:-1e-10;

				value_type c = f[i] / denom; 
				dx[0] += c * vecs[i][0];
				dx[1] += c * vecs[i][1];
				dx[2] += c * vecs[i][2];
			}

			const value_type len_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
			if (len_sq > MAX_STEP_SQ) {
				const value_type scale = std::sqrt(MAX_STEP_SQ / len_sq);
				dx[0] *= scale; dx[1] *= scale; dx[2] *= scale;
			}

			return pos_ + dx;
		}

	public:
		TYPE type_ = TYPE::N;
		PointType pos_{};
		TripleDouble value_{};
	private:
		mutable std::array<PointType, 3> prev_max_vecs_;   
		mutable bool mode_following_init_ = false;         
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

			// 1. Type B (Bond)
			for (const auto& poly : vf.polygons) {
				PointType center = (vf.polyhedra[poly.atom_ids[0]].center + vf.polyhedra[poly.atom_ids[1]].center + poly.second_shift) * 0.5;
				add_candidate(center, CriticalPoint::TYPE::B);
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

			// 3. Type C (Cage)
			for (const PointType& p : vf.vertices) {
				add_candidate(p, CriticalPoint::TYPE::C);
			}
			return result;
		}

		static inline CriticalPoint::TYPE calculateType(CriticalPoint& cp) {
			auto eigen = cp.value_.hess.EigenvaluesAndVectors();
			return static_cast<CriticalPoint::TYPE>(static_cast<int>(eigen.values[0] > 0) +
													static_cast<int>(eigen.values[1] > 0) +
													static_cast<int>(eigen.values[2] > 0));
		}

		static void OptimizePoint(CriticalPoint& cp, const PointsSoA<value_type>& psoa, value_type end_sq) {

			constexpr value_type SHRINK_FACTOR = 0.7071067811865475;
			constexpr value_type GROW_FACTOR = 1.15;


			value_type current_trust_radius = CriticalPoint::MAX_STEP;
			PointType prev_dx(0, 0, 0);
			bool has_prev = false;

			const int max_iterations = 64;
			int iter = 0;

			while (iter++ < max_iterations) {
				PointType next_point = cp.FindNextPosition();

				auto dx = next_point - cp.pos_;
				value_type drs = dx.rSq();

				if (drs < end_sq) {
					break;
				}
				if (has_prev) {
					value_type dot_dx = dx[0] * prev_dx[0] + dx[1] * prev_dx[1] + dx[2] * prev_dx[2];

					if (dot_dx < 0.0) {
						current_trust_radius = std::max(CriticalPoint::MAX_STEP * value_type(1e-5), current_trust_radius * SHRINK_FACTOR);
						next_point = cp.pos_ + dx * SHRINK_FACTOR;
						dx = next_point - cp.pos_;
					} else {
						current_trust_radius = std::min(CriticalPoint::MAX_STEP, current_trust_radius * GROW_FACTOR);
					}
				}

				prev_dx = dx;
				has_prev = true;

				cp.UpdatePos(next_point);
				cp.UpdateVal(CriticalPoint::CalculateEDinPoint(next_point, psoa, ElectronDensitySplines));
			}
			while (iter-- > 0 && cp.value_.getGradSq() > 1e-12) {
				PointType next_point = cp.NewtonRaphsonPredict();
				cp.UpdatePos(next_point);
				cp.UpdateVal(CriticalPoint::CalculateEDinPoint(next_point, psoa, ElectronDensitySplines));
			}
			auto rank = calculateType(cp);
			if (rank != cp.type_) {
				cp.type_ = rank;
			}
			return;
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

			file.write(data.data(), data.size());
		}
		static void deleteAllDublicates(std::vector<CriticalPoint>& cps, value_type eps) {
			value_type eps_sq = eps * eps;
			std::sort(cps.begin(), cps.end(), 
					  [](const CriticalPoint& a, const CriticalPoint& b) {
						  return a.pos_[0] < b.pos_[0];
					  });
			size_t cps_s = cps.size();
			std::vector<bool> is_active(cps_s, true);
			size_t j = 0;
			auto Upgrade_j = [&](size_t i) {
				j = std::max(i, j);
				auto x = cps[i].pos_[0];
				for (; j < cps_s - 1; j++)
				{
					if (std::abs(x - cps[j + 1].pos_[0]) > eps) {
						break;
					}
				}
			};
			for (size_t i = 0; i < cps_s; i++)
			{
				if (is_active[i] == false) continue;
				Upgrade_j(i);
				for (size_t k = i + 1; k <= j; k++)
				{
					if (is_active[k] == false) continue;
					if (PointType::distanceSq(cps[i].pos_, cps[k].pos_) < eps_sq)
						is_active[k] = false;
				}
			}
			size_t write_idx = 0;
			for (size_t read_idx = 0; read_idx < cps_s; ++read_idx) {
				if (is_active[read_idx]) {
					if (write_idx != read_idx) {
						cps[write_idx] = std::move(cps[read_idx]);
					}
					write_idx++;
				}
			}
			cps.resize(write_idx);
		}
	};
}
