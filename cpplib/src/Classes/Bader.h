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

#include "../Classes/Geometry.h"
#include "../Classes/Voronoi.h"

namespace cpplib {

	struct TripleDouble {
		using value_type = double;
		using grad_type = std::array<value_type, 3>;
		using hess_type = geometry::Matrix<value_type>;

		value_type val = 0.0;
		grad_type grad = {};
		hess_type hess = {};

		TripleDouble() = default;
		explicit TripleDouble(value_type v) : val(v) {
		}

		TripleDouble& operator+=(const TripleDouble& other) {
			val += other.val;
			for (size_t i = 0; i < 3; ++i) {
				grad[i] += other.grad[i];
			}
			for (size_t i = 0; i < 9; ++i) {
				hess[i] += other.hess[i];
			}
			return *this;
		}

		friend TripleDouble operator+(const TripleDouble& a, const TripleDouble& b) {
			TripleDouble result;
			result.val = a.val + b.val;
			for (size_t i = 0; i < 3; ++i) {
				result.grad[i] = a.grad[i] + b.grad[i];
			}
			for (size_t i = 0; i < 9; ++i) {
				result.hess[i] = a.hess[i] + b.hess[i];
			}
			return result;
		}
	};

	template <size_t SIZE = 576>
	struct alignas(64) CubicSpline {
	public:
		static_assert(SIZE % 64 == 0, "SIZE must be a multiple of 64 bytes");

		using value_type = TripleDouble::value_type;

		CubicSpline(value_type r_min, value_type r_max, size_t n_knots, const std::vector<value_type>& values) {
			// Validate that the incoming active nodes match the provided intervals count
			assert(values.size() == n_knots + 1);
			assert(n_knots <= SIZE); // Ensure it fits into our maximum fixed capacity
			assert(values.back() == 0.0);
			assert(r_min > 0.0 && r_max > r_min);

			// Initialize grid transformation parameters based on the ACTIVE file data boundaries
			ln_r0_ = std::log(r_min);
			r_max_ = r_max;
			m_ = n_knots; // Store active intervals count to eliminate runtime log-calculations

			// inv_delta_ is the scaling factor
			inv_delta_ = static_cast<double>(n_knots) / (std::log(r_max) - ln_r0_);
			const value_type step_delta = 1.0 / inv_delta_;

			// Fixed compile-time dimensions for stack allocations
			constexpr size_t n_max = SIZE + 1;
			constexpr size_t m_max = SIZE;

			// Runtime boundaries for the active dataset region
			const size_t n = n_knots + 1; // Active matrix size
			const size_t m = n_knots;     // Active intervals count

			// 1. Generate active nodes directly into the class array
			for (size_t i = 0; i < m; ++i) {
				knots_[i] = std::exp(std::fma(static_cast<double>(i), step_delta, ln_r0_));
			}
			// Zero-out the unused capacity tail of the knots array
			std::fill(knots_.begin() + m, knots_.end(), 0.0);

			// 2. Calculate active grid steps (intervals) using fixed-size stack memory
			std::array<value_type, m_max> h; h.fill(0.0);
			for (size_t i = 0; i < m - 1; ++i) {
				h[i] = knots_[i + 1] - knots_[i];
				assert(h[i] > 0);
			}
			h[m - 1] = r_max_ - knots_[m - 1];
			assert(h[m - 1] > 0);

			// Build the tridiagonal system for the active second derivatives on the stack
			std::array<value_type, n_max> sub;  sub.fill(0.0);
			std::array<value_type, n_max> diag; diag.fill(0.0);
			std::array<value_type, n_max> sup;  sup.fill(0.0);
			std::array<value_type, n_max> rhs;  rhs.fill(0.0);

			diag[0] = 1.0; // Natural spline boundary condition (S''(x0) = 0)

			for (size_t i = 1; i < m; ++i) {
				sub[i] = h[i - 1];
				diag[i] = 2.0 * (h[i - 1] + h[i]);
				sup[i] = h[i];
				rhs[i] = 6.0 * ((values[i + 1] - values[i]) / h[i] -
								(values[i] - values[i - 1]) / h[i - 1]);
			}

			// Right boundary condition: First derivative is zero (S'(x_{n-1}) = 0)
			const size_t last = n - 1;
			sub[last] = h[m - 1];
			diag[last] = 2.0 * h[m - 1];
			sup[last] = 0.0;
			rhs[last] = 6.0 * values[last - 1] / h[m - 1];

			// Solve the system using Thomas algorithm (TDMA) strictly on stack arrays
			std::array<value_type, n_max> gamma; gamma.fill(0.0);
			std::array<value_type, n_max> beta;  beta.fill(0.0);

			gamma[0] = sup[0] / diag[0];
			beta[0] = rhs[0] / diag[0];

			for (size_t i = 1; i < n; ++i) {
				value_type denom = diag[i] - sub[i] * gamma[i - 1];
				gamma[i] = sup[i] / denom;
				beta[i] = (rhs[i] - sub[i] * beta[i - 1]) / denom;
			}

			// Back-substitution to find second derivatives
			std::array<value_type, n_max> second_deriv; second_deriv.fill(0.0);
			second_deriv[last] = beta[last];
			for (size_t i = last; i-- > 0; ) {
				second_deriv[i] = beta[i] - gamma[i] * second_deriv[i + 1];
			}

			// Calculate and store polynomial coefficients into the active part of class arrays
			for (size_t i = 0; i < m; ++i) {
				value_type hi = h[i];
				value_type yi = values[i];
				value_type yip1 = values[i + 1];
				value_type mi = second_deriv[i];
				value_type mip1 = second_deriv[i + 1];

				a_[i] = (mip1 - mi) / (6.0 * hi);
				b_[i] = mi / 2.0;
				c_[i] = (yip1 - yi) / hi - hi * (2.0 * mi + mip1) / 6.0;
				d_[i] = yi;
			}

			// 3. Safely zero-out the remaining unused capacity of the member arrays
			std::fill(a_.begin() + m, a_.end(), value_type(0.0));
			std::fill(b_.begin() + m, b_.end(), value_type(0.0));
			std::fill(c_.begin() + m, c_.end(), value_type(0.0));
			std::fill(d_.begin() + m, d_.end(), value_type(0.0));
		}


		/// @brief Evaluate spline value and its first two derivatives at point r
		/// @param r Point to evaluate
		/// @param phi returns value of spline at r
		/// @param dphi returns value of 1-st derivative of spline at r
		/// @param ddphi returns value of 2-nd derivative of spline at r
		void eval(value_type r, value_type& phi, value_type& dphi, value_type& ddphi) const {
			if (r < knots_[0] || r >= radius()) {
				phi = dphi = ddphi = 0.0;
				return;
			}
			auto idx = find_active_spline(r);

			// Safety clamp to ensure idx is within coeffs_ range
			if (idx >= SIZE) [[unlikely]] {
				phi = dphi = ddphi = 0.0;
				return;
			}

			const auto a = a_[idx];
			const auto b = b_[idx];
			const auto c = c_[idx];
			const auto d = d_[idx];
			value_type dr = r - knots_[idx];

			// Horner's method for stable evaluation
			phi = std::fma(std::fma(std::fma(dr, a, b), dr, c), dr, d);
			dphi = std::fma(std::fma(3.0 * dr, a, 2.0 * b), dr, c);
			ddphi = std::fma(6.0 * dr, a, 2.0 * b);
		}

		/// @brief Find the point r where the spline value matches target_val
		/// @param start Initial guess (guaranteed to be closer to 0 than the true root, meaning f(start) > target_val)
		/// @param target_val The function value to search for
		/// @return The point r, or radius() if the value is out of bounds or decayed
		value_type find_value(value_type start, value_type target_val) const {
			// 1. Fast boundary checks
			if (target_val >= d_[0]) return knots_[0];
			if (target_val <= 0.0)   return radius(); // Since values.back() == 0.0

			// 2. Hybrid O(1) jump followed by local correction
			size_t idx = find_active_spline_hybrid(start, target_val);

			// Safety clamp
			if (idx >= SIZE) [[unlikely]] return radius();

			// 3. Compute an extremely accurate initial guess inside the interval.
			value_type f_left = d_[idx];
			value_type f_right = (idx + 1 < SIZE)?d_[idx + 1]:0.0;

			// If the right node is zero (edge of the grid), protect against log(0)
			if (f_right <= 0.0) f_right = 1e-300;

			value_type t = (std::log(target_val) - std::log(f_left)) / (std::log(f_right) - std::log(f_left));

			// Starting r value
			value_type r_left = knots_[idx];
			value_type r_right = (idx + 1 < SIZE)?knots_[idx + 1]:r_max_;

			// Linear interpolation in log-space for the initial r guess
			value_type r = std::exp(std::log(r_left) + t * (std::log(r_right) - std::log(r_left)));

			// 4. Halley's method (Cubic convergence rate)
			const auto a = a_[idx];
			const auto b = b_[idx];
			const auto c = c_[idx];
			const auto d = d_[idx];

			for (int iter = 0; iter < 2; ++iter) {
				value_type f_val, df, d2f;
				eval(r, f_val, df, d2f); 
				
				value_type df_log = r * df;
				value_type d2f_log = std::fma(r * r, d2f, df_log);

				value_type numerator = 2.0 * f * df_log;
				value_type denominator = std::fma(2.0 * df_log, df_log, -(f * d2f_log));

				value_type du = numerator / denominator;
				r *= std::exp(-du);
			}

			return r;
		}

		constexpr value_type radius() const {
			return r_max_;
		}

	private:
		// O(1)-optimized finder of active spline
		size_t find_active_spline(value_type r) const {
			return static_cast<size_t>((std::log(r) - ln_r0_) * inv_delta_);
		}

		/// @brief Bounded O(1) parabolic log-jump to the target neighborhood with branchless-friendly correction
		size_t find_active_spline_hybrid(value_type start, value_type target_val) const {
			// Minimum index dictated by the 'start' guarantee (O(1) grid transformation)
			size_t idx_start = find_active_spline(start);
			if (idx_start >= m_ - 1) return m_ - 1;

			// Values at the boundaries of the REMAINDER of the grid
			value_type val_first = d_[idx_start];
			value_type val_last = d_[m_ - 1];

			// Strict monotonicity defense against numerical artifacts
			if (target_val >= val_first) return idx_start;

			// 1. Three-point logarithmic mapping to capture the changing decay rate.
			// We sample the start, the end, and the exact middle node of the remaining grid.
			size_t idx_mid = idx_start + ((m_ - 1) - idx_start) / 2;
			value_type val_mid = d_[idx_mid];

			// SAFETY: clamp values to avoid NaN / -inf
			value_type safe_target = (target_val > 1e-300)?target_val:1e-300;
			value_type safe_first = (val_first > 1e-300)?val_first:1e-300;
			value_type safe_mid = (val_mid > 1e-300)?val_mid:1e-300;
			value_type safe_last = (val_last > 1e-300)?val_last:1e-300;

			value_type ln_y = std::log(safe_target);
			value_type ln_y0 = std::log(safe_first);
			value_type ln_ym = std::log(safe_mid);
			value_type ln_yN = std::log(safe_last);

			// Map log-values to a normalized [0, 1] range relative to the midpoint
			value_type d1 = ln_ym - ln_y0;
			value_type d2 = ln_yN - ln_y0;
			value_type dy = ln_y - ln_y0;

			// Avoid division by zero if the remaining grid is too small
			value_type norm_pos = 0.0;
			if (std::abs(d2 * (d1 - 0.5 * d2)) > 1e-15) {
				// Quadratic fit (Parabolic inverse interpolation): x = a*y^2 + b*y
				// Calculates how far to jump through the transition region
				value_type cA = (d1 - 0.5 * d2) / (d1 * d2 * (0.5 * d1 - 0.5 * d2 + 1e-300)); // Quadratic coefficient
				value_type cB = (1.0 - cA * d2 * d2) / d2;                                    // Linear coefficient
				norm_pos = (cA * dy + cB) * dy;
			} else {
				// Fallback to linear log-interpolation if the segment is nearly linear
				norm_pos = dy / d2;
			}

			norm_pos = std::max(0.0, std::min(norm_pos, 1.0));

			// Convert normalized position back to grid index units
			size_t remaining_intervals = (m_ - 1) - idx_start;
			value_type guessed_offset = norm_pos * static_cast<value_type>(remaining_intervals);

			auto guessed_idx = static_cast<size_t>(idx_start + static_cast<size_t>(guessed_offset));
			guessed_idx = std::max(idx_start, std::min(guessed_idx, m_ - 1));

			// 2. Micro-tuning correction loops.
			// Since the quadratic fit maps the smooth transition curve with high fidelity,
			// these loops will evaluate instantly, often performing 0 iterations.
			while (guessed_idx < (m_ - 1) && d_[guessed_idx + 1] > target_val) {
				guessed_idx++;
			}
			while (guessed_idx > idx_start && d_[guessed_idx] < target_val) {
				guessed_idx--;
			}

			return guessed_idx;
		}

		std::array<value_type, SIZE> a_;
		std::array<value_type, SIZE> b_;
		std::array<value_type, SIZE> c_;
		std::array<value_type, SIZE> d_;
		std::array<value_type, SIZE> knots_; // knot_[SIZE] is r_max_

		value_type r_max_; // x = knot[last]
		value_type ln_r0_; // ln(x_0)
		value_type inv_delta_; // 1.0 / ln (x_1 / x_0)
		size_t m_; // Precalculated active intervals count to optimize hot path
	};

	class RadialSpline {
	public:
		using value_type = CubicSpline<>::value_type;
		using PointType = geometry::Point<value_type>;

		// Constructs the underlying CubicSpline directly in-place without any copying/moving
		RadialSpline(value_type r_min, value_type r_max, size_t n_knots, const std::vector<value_type>& values)
			: spline(r_min, r_max, n_knots, values) {
		}

		/// @brief Evaluates value, gradient, and Hessian for a radially symmetric spline
		/// @param x_minus_c is the relative vector from the center (r = x - c)
		/// @return TripleDouble value at point x_minus_c
		TripleDouble evaluate(const PointType& x_minus_c) const {
			value_type r2 = x_minus_c.rSq();
			value_type r = std::sqrt(r2);
			if (r >= spline.radius()) {
				return TripleDouble{}; // Return zero-initialized struct
			}
			value_type phi;
			value_type dphi;
			value_type ddphi;
			spline.eval(r, phi, dphi, ddphi);

			TripleDouble result;
			result.val = phi;

			if (r > voronoi::EPSILON) {
				value_type inv_r = 1.0 / r;
				value_type inv_r2 = inv_r * inv_r;
				value_type dphi_over_r = dphi * inv_r;
				value_type coeff_hess = ddphi - dphi_over_r;
				value_type coeff_hess_over_r2 = coeff_hess * inv_r2;

				for (uint8_t i = 0; i < 3; ++i) {
					result.grad[i] = dphi_over_r * x_minus_c[i];
					value_type coeff_i = coeff_hess_over_r2 * x_minus_c[i];
					for (uint8_t j = 0; j < 3; ++j) {
						result.hess[i * 3 + j] = coeff_i * x_minus_c[j];
					}
					result.hess[i * 3 + i] += dphi_over_r;
				}
			}

			return result;
		}

	public:
		CubicSpline<> spline;
	};

	template <typename T, size_t NEAR = 2>
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


		using PointType = geometry::Point<T>;

		static constexpr size_t TOTAL_SHIFTS = (NEAR * 2 + 1) * (NEAR * 2 + 1) * (NEAR * 2 + 1);

		template<typename I, size_t N>
		static consteval auto unrollPositions() {
			constexpr size_t dim = 2 * N + 1;
			constexpr size_t size = dim * dim * dim;
			using LocalPointType = typename geometry::Point<I>;
			std::array<LocalPointType, size> result;

			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					for (int k = 0; k < dim; k++) {
						result[i * dim * dim + j * dim + k] =
							LocalPointType(i - static_cast<int>(N),
										   j - static_cast<int>(N),
										   k - static_cast<int>(N));
					}
				}
			}
			return result;
		}

		static constexpr auto p_near = unrollPositions<char, NEAR>();

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
		void calculateSpatialIndexes(T one_over_period) {
			const uint32_t size = static_cast<uint32_t>(x.size());

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

		constexpr BoundsArray<int32_t> getBoundsFrac(const geometry::Matrix<T> mat_CartToFrac, T cutoff) const {
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

	};

	class CriticalPoint {
	public:
		using PointType = typename RadialSpline::PointType;
		static constexpr PointType::value_type MAX_STEP = 0.2;
		static constexpr PointType::value_type MAX_STEP_SQ = MAX_STEP * MAX_STEP;
		static constexpr double EPS = 1e-12;

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
					assert(false);
					return pos_;
				case B:
					return EigenVectorFollowing();
				case R:
					return EigenVectorFollowing();
				case C:
					return NewtonRaphsonPredict();
			}
			return pos_;
		}
		void UpdatePoint(const PointType& pos, const TripleDouble& value) {
			pos_ = pos;
			value_ = value;
		}
	private:
		PointType EigenVectorFollowing() const {
			const auto& g = value_.grad;
			const auto& h = value_.hess;

			auto eig = h.EigenvaluesAndVectors();

			
			int idx = 0;
			if (type_ == TYPE::B) {
				double max_val = eig.values[0];
				if (eig.values[1] > max_val) {
					max_val = eig.values[1]; 
					idx = 1;
				}
				if (eig.values[2] > max_val) {
					idx = 2;
				}
			} else {
				double min_val = eig.values[0];
				if (eig.values[1] < min_val) {
					min_val = eig.values[1]; 
					idx = 1;
				}
				if (eig.values[2] < min_val) {
					idx = 2;
				}
			}

			const auto& v = eig.vectors[idx];
			double lambda = eig.values[idx];

			double grad_dot_v = g[0] * v[0] + g[1] * v[1] + g[2] * v[2];

			double alpha = 0.0;

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

			double det = h.Det();
			if (std::abs(det) < EPS) {
				double step = 0.01;
				return pos_ + PointType(-g[0] * step, -g[1] * step, -g[2] * step);
			}

			const auto H_inv = h.Invert();
			double delta_x = -(H_inv.El(0, 0) * g[0] + H_inv.El(0, 1) * g[1] + H_inv.El(0, 2) * g[2]);
			double delta_y = -(H_inv.El(1, 0) * g[0] + H_inv.El(1, 1) * g[1] + H_inv.El(1, 2) * g[2]);
			double delta_z = -(H_inv.El(2, 0) * g[0] + H_inv.El(2, 1) * g[1] + H_inv.El(2, 2) * g[2]);

			const double alpha_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
			if (alpha_sq > MAX_STEP_SQ) { // TODO: Maybe add [[unlikely]]
				const double step_over_alpha = std::sqrt(MAX_STEP_SQ / alpha_sq);

				delta_x *= step_over_alpha;
				delta_y *= step_over_alpha;
				delta_z *= step_over_alpha;
			}

			return pos_ + PointType(delta_x, delta_y, delta_z);
		}

	private:
		TYPE type_ = TYPE::N;
		PointType pos_{};
		TripleDouble value_{};
	};

	class BaderOperator {
	public:
		using value_type = CubicSpline<>::value_type;
		static PointsSoA<value_type> CreatePointGrid(const voronoi::VoronoiFused& vf, value_type radius) {
			PointsSoA<value_type> grid;
			// TODO
			

			return grid;
		}

		static value_type GetOptimalRadius(value_type radius, value_type vertical_eps, const std::vector<RadialSpline>& splines, const std::vector<uint8_t>& active_elems) {
			value_type max_x = 0;
			const size_t active_elems_size = active_elems.size();
			for (size_t i = 0; i < active_elems_size; i++)
			{
				value_type phi;
				value_type dphi;
				value_type ddphi;
				uint8_t active_index = active_elems[i];
				auto& active_spline = splines[active_index].spline;
				active_spline.eval(radius, phi, dphi, ddphi);
				value_type x = active_spline.find_value(radius, phi * vertical_eps);
				max_x = std::max(x, max_x);
			}
			return max_x;
		}



		// TODO
	};



}
