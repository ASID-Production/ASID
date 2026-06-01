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
		using hess_type = std::array<value_type, 9>;

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

	struct CubicSpline {
	public:
		using value_type = TripleDouble::value_type;
		struct IntervalCoeffs {
			value_type a; // x^3 coefficient
			value_type b; // x^2 coefficient
			value_type c; // x^1 coefficient
			value_type d; // constant term
		};

		CubicSpline(std::vector<value_type>&& knots, const std::vector<value_type>& values) : knots_(std::move(knots)) {
			assert(knots_.size() == values.size());
			assert(knots_.size() >= 2);
			assert(values.back() == 0.0);   // Required by the problem statement

			ln_r0_ = std::log(knots_[0]);

			const size_t n = knots_.size();          // Number of nodes
			const size_t m = n - 1;                  // Number of intervals

			// Grid steps (intervals)
			std::vector<value_type> h(m);
			for (size_t i = 0; i < m; ++i) {
				h[i] = knots_[i + 1] - knots_[i];
				assert(h[i] > 0);
			}

			// Build tri-diagonal system for second derivatives S''(x_i)
			std::vector<value_type> sub(n, 0.0);  // sub-diagonal (a_i)
			std::vector<value_type> diag(n, 0.0); // main diagonal (b_i)
			std::vector<value_type> sup(n, 0.0);  // super-diagonal (c_i)
			std::vector<value_type> rhs(n, 0.0);  // right-hand side (d_i)

			// Left boundary condition: S''(x0) = 0 (Natural spline)
			diag[0] = 1.0;
			sup[0] = 0.0;
			rhs[0] = 0.0;

			// Internal nodes i = 1 .. n-2
			for (size_t i = 1; i < n - 1; ++i) {
				sub[i] = h[i - 1];
				diag[i] = 2.0 * (h[i - 1] + h[i]);
				sup[i] = h[i];
				rhs[i] = 6.0 * ((values[i + 1] - values[i]) / h[i] -
								(values[i] - values[i - 1]) / h[i - 1]);
			}

			// Right boundary condition: S'(x_{n-1}) = 0
			// Derived from: h_{m-1}/6 * M_{n-2} + h_{m-1}/3 * M_{n-1} = -(y_{n-1}-y_{n-2})/h_{m-1}
			const size_t last = n - 1;
			sub[last] = h[m - 1];                 // Coefficient for M_{n-2}
			diag[last] = 2.0 * h[m - 1];          // Coefficient for M_{n-1}
			sup[last] = 0.0;
			rhs[last] = 6.0 * values[last - 1] / h[m - 1]; // Equivalent to -6*(0 - y)/h since values.back() == 0

			// Solve the system using Thomas algorithm (TDMA)
			std::vector<value_type> gamma(n, 0.0), beta(n, 0.0);
			gamma[0] = sup[0] / diag[0];
			beta[0] = rhs[0] / diag[0];

			for (size_t i = 1; i < n; ++i) {
				value_type denom = diag[i] - sub[i] * gamma[i - 1];
				gamma[i] = sup[i] / denom;
				beta[i] = (rhs[i] - sub[i] * beta[i - 1]) / denom;
			}

			// Back-substitution
			std::vector<value_type> second_deriv(n, 0.0);
			second_deriv[last] = beta[last];
			for (size_t i = last; i-- > 0; ) {
				second_deriv[i] = beta[i] - gamma[i] * second_deriv[i + 1];
			}

			// Calculate cubic polynomial coefficients for each interval
			coeffs_.resize(m);
			for (size_t i = 0; i < m; ++i) {
				value_type hi = h[i];
				value_type yi = values[i];
				value_type yip1 = values[i + 1];
				value_type mi = second_deriv[i];
				value_type mip1 = second_deriv[i + 1];

				// Polynomial S_i(dr) = a*dr^3 + b*dr^2 + c*dr + d, where dr = x - x_i
				coeffs_[i].a = (mip1 - mi) / (6.0 * hi);
				coeffs_[i].b = mi / 2.0;
				coeffs_[i].c = (yip1 - yi) / hi - hi * (2.0 * mi + mip1) / 6.0;
				coeffs_[i].d = yi;
			}
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
			if (idx >= coeffs_.size()) idx = coeffs_.size() - 1;

			value_type dr = r - knots_[idx];
			const auto& c = coeffs_[idx];

			// Horner's method for stable evaluation
			phi = std::fma(std::fma(std::fma(dr, c.a, c.b), dr, c.c), dr, c.d);
			dphi = std::fma(std::fma(3.0 * dr, c.a, 2.0 * c.b), dr, c.c);
			ddphi = std::fma(6.0 * dr, c.a, 2.0 * c.b);
		}

		value_type radius() const {
			return knots_.back();
		}

	private:
		// O(1)-optimized finder of active spline
		size_t find_active_spline(value_type r) const {
			constexpr double INV_ALPHA = 40.98379665578782;
			return static_cast<size_t>((std::log(r) - ln_r0_) * INV_ALPHA);
		}

		std::vector<value_type> knots_;
		std::vector<IntervalCoeffs> coeffs_;
		value_type ln_r0_;
	};

	class RadialSpline {
	public:
		using value_type = CubicSpline::value_type;
		using PointType = geometry::Point<value_type>;

		explicit RadialSpline(CubicSpline&& spline) : spline_(std::move(spline)) {
		}

		/// @brief Evaluates value, gradient, and Hessian for a radially symmetric spline
		/// @param x_minus_c is the relative vector from the center (r = x - c)
		/// @return TripleDouble value at point x_minus_c
		TripleDouble evaluate(const PointType& x_minus_c) const {
			value_type r2 = x_minus_c.rSq();
			value_type r = std::sqrt(r2);
			if (r >= spline_.radius()) {
				return TripleDouble{}; // Return zero-initialized struct
			}
			value_type phi;
			value_type dphi;
			value_type ddphi;
			spline_.eval(r, phi, dphi, ddphi);

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

	private:
		CubicSpline spline_;
	};

	template <typename T>
	struct PointsSoA {
		std::vector<T> x;
		std::vector<T> y;
		std::vector<T> z;
		std::vector<uint32_t> ix;
		std::vector<uint32_t> iy;
		std::vector<uint32_t> iz;

		std::vector<uint32_t> ids;
		std::vector<uint32_t> spine_ids;

		size_t min_ix = 0;
		size_t max_ix = 0;

		size_t min_iy = 0;
		size_t max_iy = 0;

		size_t min_iz = 0;
		size_t max_iz = 0;

		using PointType = geometry::Point<T>;
		
		static constexpr size_t NEAR = 2;

		template<typename I, size_t N>
		static consteval auto unrollPositions() {
			constexpr size_t dim = 2 * N + 1;
			constexpr size_t size = dim * dim * dim;
			using LocalPointType = typename geometry::Point<I>;
			std::array<LocalPointType, size> result;

			for (int i = -static_cast<int>(N); i <= N; i++) {
				for (int j = -static_cast<int>(N); j <= N; j++) {
					for (int k = -static_cast<int>(N); k <= N; k++) {
						result[i * dim * dim + j * dim + k] = LocalPointType(i, j, k);
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
			ix.reserve(capacity);
			iy.reserve(capacity);
			iz.reserve(capacity);

			ids.reserve(capacity);
			spine_ids.reserve(capacity);
		}
		void addPoint(PointType point, uint32_t id, uint32_t spine_id) {
			x.push_back(point[0]);
			y.push_back(point[1]);
			z.push_back(point[2]);

			ids.push_back(id);
			spine_ids.push_back(spine_id);
		}
		void calculateSpartialIndexes(T one_over_period) {
			const size_t size = x.size();
			for (size_t i = 0; i < size; ++i) {
				ix.push_back(static_cast<uint32_t>(std::floor(x[i] * one_over_period)));
				iy.push_back(static_cast<uint32_t>(std::floor(y[i] * one_over_period)));
				iz.push_back(static_cast<uint32_t>(std::floor(z[i] * one_over_period)));
			}
		}
		static std::array<uint32_t, 6> calculateMinMaxIndexes(geometry::Matrix<T>& mat, T cutoff) {
			constexpr std::array<PointType, 8> p01 = {{
				{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
				{1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}
			}};

			std::array<PointType, 8> real;
			for (int i = 0; i < 8; ++i) {
				real[i] = mat * p01[i];
			}

			std::array<T, 6> result = {0, 1, 0, 1, 0, 1};

			// TODO: finish
		}
	};


	class CriticalPoint {
	public:
		using PointType = typename RadialSpline::PointType;
		enum class TYPE {
			N = 0,
			B = 1,
			R = 2,
			C = 3
		};
		CriticalPoint() = default;
		CriticalPoint(TYPE type, const PointType& pos) : type_(type), pos_(pos) {
		}
		PointType FindNextPosition() const {
			switch (type_) {
				using enum TYPE;
				case N:
					assert(false);
					return pos_;
					break;
				case B:
					return EigenVectorFollowing();
					break;
				case R:
					return EigenVectorFollowing();
					break;
				case C:
					return NewtonRaphsonPredict();
					break;
			}
			return pos_;
		}
		void UpdatePoint(const PointType& pos, const TripleDouble& value) {
			pos_ = pos;
			value_ = value;
		}
	private:
		PointType EigenVectorFollowing() const {
			// TODO
			return pos_;
		}

		PointType NewtonRaphsonPredict() const {
			// TODO
			return pos_;
		}

	private:
		TYPE type_ = TYPE::N;
		PointType pos_{};
		TripleDouble value_{};
	};

	class BaderOperator {
		// TODO
	};



}
