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
		using grad_type = geometry::Point<value_type>;
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
		constexpr value_type getGradSq() const noexcept {
			return grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
		}
	};

	template <size_t N = 576>
	struct alignas(64) CubicSpline {
	public:
		static constexpr size_t SIZE = N;
		static_assert(SIZE % 64 == 0, "Template Argument must be a multiple of 64 bytes");

		using value_type = TripleDouble::value_type;

		constexpr CubicSpline() = default;
		CubicSpline(value_type r_min, value_type r_max, size_t n_knots, const std::vector<value_type>& values) {
			create(r_min, r_max, n_knots, values);
		}

		void create(value_type r_min, value_type r_max, size_t n_knots, const std::vector<value_type>& values) {
			// Validate that the incoming active nodes match the provided intervals count
			assert(values.size() == n_knots + 1);
			assert(n_knots <= SIZE); // Ensure it fits into our maximum fixed capacity
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
			if (target_val >= d_[0]) return knots_[0];
			if (target_val <= 0.0)   return radius();

			size_t idx = find_active_spline_hybrid(start, target_val);

			const value_type r_left = knots_[idx];
			const value_type r_right = knots_[idx + 1];
			const value_type h_node = r_right - r_left;

			const value_type f_left = d_[idx];
			const value_type f_right = d_[idx + 1];

			value_type t = (target_val - f_left) / (f_right - f_left + 1e-300);
			t = std::clamp(t, 0.0, 1.0);

			value_type dr = t * h_node;

			const auto a = a_[idx];
			const auto b = b_[idx];
			const auto c = c_[idx];
			const auto d = d_[idx];

			for (int iter = 0; iter < 2; ++iter) {
				value_type f_val = std::fma(std::fma(std::fma(dr, a, b), dr, c), dr, d) - target_val;
				value_type df = std::fma(std::fma(3.0 * dr, a, 2.0 * b), dr, c);
				if (std::abs(df) < 1e-15) [[unlikely]] {
					break;
				}

				dr -= f_val / df;
				dr = std::clamp(dr, 0.0, h_node);
			}

			return r_left + dr;
		}

		constexpr value_type radius() const {
			return r_max_;
		}
		constexpr value_type r0() const {
			return knots_[0];
		}

	private:
		// O(1)-optimized finder of active spline
		size_t find_active_spline(value_type r) const {
			return static_cast<size_t>((std::log(r) - ln_r0_) * inv_delta_);
		}

		size_t find_active_spline_hybrid(value_type start, value_type target_val) const {
			size_t idx_start = find_active_spline(start);
			if (idx_start >= m_ - 1) return m_ - 1;

			if (target_val >= d_[idx_start]) return idx_start;

			auto first = d_.begin() + idx_start;
			auto last = d_.begin() + m_; 

			auto it = std::upper_bound(first, last, target_val, std::greater<value_type>());

			size_t idx = std::distance(d_.begin(), it);

			if (idx > idx_start) {
				idx--;
			}

			return std::clamp(idx, idx_start, m_ - 2);
		}


		std::array<value_type, SIZE> a_{};
		std::array<value_type, SIZE> b_{};
		std::array<value_type, SIZE> c_{};
		std::array<value_type, SIZE> d_{};
		std::array<value_type, SIZE> knots_{}; // knot_[SIZE] is r_max_

		value_type r_max_ = 0; // x = knot[last]
		value_type ln_r0_ = 1; // ln(x_0)
		value_type inv_delta_ = 1; // 1.0 / ln (x_1 / x_0)
		size_t m_ = 0; // Precalculated active intervals count to optimize hot path
	};

	class RadialSpline {
	public:
		using value_type = CubicSpline<>::value_type;
		using PointType = geometry::Point<value_type>;
		static constexpr size_t SIZE = CubicSpline<>::SIZE;
		// Constructs the underlying CubicSpline directly in-place without any copying/moving
		RadialSpline(value_type r_min, value_type r_max, size_t n_knots, const std::vector<value_type>& values)
			: spline(r_min, r_max, n_knots, values) {}
		constexpr RadialSpline() = default;
		/// @brief Evaluates value, gradient, and Hessian for a radially symmetric spline
		/// @param x_minus_c is the relative vector from the center (r = x - c)
		/// @return TripleDouble value at point x_minus_c
		TripleDouble evaluate(const PointType& x_minus_c) const {
			const value_type r = x_minus_c.r();
			if ((r >= spline.radius()) || (r <= spline.r0())) {
				return TripleDouble{}; // Return zero-initialized struct
			}
			value_type phi;
			value_type dphi;
			value_type ddphi;
			spline.eval(r, phi, dphi, ddphi);

			TripleDouble result(phi);

			value_type inv_r = 1.0 / r;
			value_type inv_r2 = inv_r * inv_r;
			value_type dphi_over_r = dphi * inv_r;
			value_type coeff_hess = ddphi - dphi_over_r;
			value_type coeff_hess_over_r2 = coeff_hess * inv_r2;

			for (uint8_t i = 0; i < 3; ++i) {
				result.grad[i] = dphi_over_r * x_minus_c[i];
				value_type coeff_i = coeff_hess_over_r2 * x_minus_c[i];
				for (uint8_t j = 0; j < 3; ++j) {
					result.hess.El(i, j) = coeff_i * x_minus_c[j];
				}
				result.hess.El(i, i) += dphi_over_r;
			}
			return result;
		}

	public:
		CubicSpline<> spline;
	};
}
