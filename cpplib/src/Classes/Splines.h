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
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include "../Classes/Geometry.h"

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
            // Trusted HPC input: ED metadata/values are checked by DensityParser.
            std::vector<value_type> knots(n_knots + 1);
            const auto log_min = std::log(r_min);
            const auto delta = (std::log(r_max) - log_min) / n_knots;
            knots.front() = r_min;
            for (size_t i = 1; i < n_knots; ++i)
                knots[i] = std::exp(std::fma(static_cast<double>(i), delta, log_min));
            knots.back() = r_max;
            create(knots, values);
        }
        CubicSpline(const std::vector<value_type>& knots, const std::vector<value_type>& values) {
            create(knots, values);
        }
        void create(const std::vector<value_type>& knots, const std::vector<value_type>& values) {
            const size_t n_knots = knots.size() - 1;
            m_ = n_knots;
            r_max_ = knots.back();
            std::copy(knots.begin(), knots.end(), knots_.begin());
            std::fill(knots_.begin() + knots.size(), knots_.end(), 0.0);

			// Fixed compile-time dimensions for stack allocations
			constexpr size_t n_max = SIZE + 1;
			constexpr size_t m_max = SIZE;

			// Runtime boundaries for the active dataset region
			const size_t n = n_knots + 1; // Active matrix size
			const size_t m = n_knots;     // Active intervals count

			// 2. Calculate active grid steps (intervals) using fixed-size stack memory
			std::array<value_type, m_max> h; h.fill(0.0);
			for (size_t i = 0; i < m - 1; ++i) {
				h[i] = knots_[i + 1] - knots_[i];
			}
			h[m - 1] = r_max_ - knots_[m - 1];

			// Build the tridiagonal system for the active second derivatives on the stack
			std::array<value_type, n_max> sub;  sub.fill(0.0);
			std::array<value_type, n_max> diag; diag.fill(0.0);
			std::array<value_type, n_max> sup;  sup.fill(0.0);
			std::array<value_type, n_max> rhs;  rhs.fill(0.0);

			// At an explicit zero node impose radial smoothness S'(0)=0.
            // Positive-radius ED grids retain their natural left boundary.
            diag[0] = knots_[0] == 0 ? 2 * h[0] : 1;
            if (knots_[0] == 0) {
                sup[0] = h[0];
                rhs[0] = 6 * (values[1] - values[0]) / h[0];
            }

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
			rhs[last] = 6.0 * (values[last - 1] - values[last]) / h[m - 1];

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
                if (i == 0 && knots_[0] == 0) c_[i] = 0; // exact clamped boundary
				d_[i] = yi;
			}

			d_[m] = values[m];

			// 3. Safely zero-out the remaining unused capacity of the member arrays
			std::fill(a_.begin() + m, a_.end(), value_type(0.0));
			std::fill(b_.begin() + m, b_.end(), value_type(0.0));
			std::fill(c_.begin() + m, c_.end(), value_type(0.0));
			std::fill(d_.begin() + m + 1, d_.end(), value_type(0.0));
		}

		/// @brief Evaluate spline value and its first two derivatives at point r
		/// @param r Point to evaluate
		/// @param phi returns value of spline at r
		/// @param dphi returns value of 1-st derivative of spline at r
		/// @param ddphi returns value of 2-nd derivative of spline at r
		void eval(value_type r, value_type& phi, value_type& dphi, value_type& ddphi) const {
            if (m_ == 0 || r > radius()) {
                phi = dphi = ddphi = 0;
                return;
            }
            if (r < knots_[0]) {
                phi = d_[0];
                dphi = ddphi = 0;
                return;
            }
            const auto idx = find_active_spline(r);

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
            if (m_ == 0) return radius();
            if (target_val >= d_[0]) return knots_[0];
            if (target_val <= d_[m_]) return radius();
            if (start >= radius()) return radius();
            for (size_t i = find_active_spline(start); i < m_; ++i) {
                if (d_[i] < target_val || d_[i+1] > target_val) continue;
                value_type left = knots_[i], right = knots_[i+1];
                for (int iteration = 0; iteration <= std::numeric_limits<value_type>::digits; ++iteration) {
                    const auto mid = left + (right-left) / 2;
                    if (mid == left || mid == right) break;
                    value_type value, first, second;
                    eval(mid, value, first, second);
                    if (value > target_val) left = mid; else right = mid;
                }
                return left + (right-left) / 2;
            }
            return radius();
        }

		constexpr value_type radius() const {
			return r_max_;
		}
		constexpr value_type r0() const {
			return knots_[0];
		}

	private:
        size_t find_active_spline(value_type r) const {
            // Use the first interval at/below its node; use the last
            // polynomial at the endpoint. Interior queries use their actual knots.
            if (r <= knots_[0]) return 0;
            if (r >= r_max_) return m_ - 1;
            return static_cast<size_t>(std::upper_bound(knots_.begin(), knots_.begin()+m_+1, r)
                                       - knots_.begin() - 1);
        }

		std::array<value_type, SIZE> a_{};
		std::array<value_type, SIZE> b_{};
		std::array<value_type, SIZE> c_{};
		std::array<value_type, SIZE + 1> d_{};
		std::array<value_type, SIZE + 1> knots_{}; // knot_[SIZE] is r_max_

		value_type r_max_ = 0; // x = knot[last]
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
			const value_type r = std::hypot(x_minus_c[0], x_minus_c[1], x_minus_c[2]);
			if (r > spline.radius()) {
				return TripleDouble{}; // Return zero-initialized struct
			}
			value_type phi;
			value_type dphi;
			value_type ddphi;
			spline.eval(r, phi, dphi, ddphi);

			TripleDouble result(phi);

            // A constant core contributes zero gradient and Hessian. Other atoms'
            // contributions are still accumulated by the caller.
            if (r < spline.r0()) return result;
            if (r == 0) {
                for (size_t i = 0; i < 3; ++i) result.hess.El(i,i) = ddphi;
                return result;
            }
            // Unit-vector form avoids 1/r^2 overflow. There is no epsilon sphere:
            // every positive radius uses the radial derivatives of the same spline.
            const auto u = x_minus_c / r;
            const auto tangential = dphi / r;
            for (size_t i = 0; i < 3; ++i) {
                result.grad[i] = dphi * u[i];
                for (size_t j = 0; j < 3; ++j)
                    result.hess.El(i,j) = (ddphi - tangential) * u[i] * u[j];
                result.hess.El(i,i) += tangential;
            }
			return result;
		}

	public:
		CubicSpline<> spline;
	};
}
