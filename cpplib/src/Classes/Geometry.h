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
#include <cstdint> // for int8_t etc.
#include <type_traits>
#include <array>
#include <vector>
#include <list> // for std::list
#include <cmath> // for std::sqrt and std::floor
#include <utility> // for std::move
#include <algorithm> // for std::sort
#include <ranges>
#include <numbers>
#include <optional>
#include <concepts>

#include "../BaseHeaders/Concepts.h"

namespace cpplib::geometry {
	template <class T> inline T GradtoRad(T a) { return a * static_cast<T>(0.0174532925199432957692); }
	template <class T> inline T RadtoGrad(T a) { return a * static_cast<T>(57.295779513082320877); }

	constexpr double crystallography_eq_position_eps_realspace = 0.01; // Angstrom
	constexpr double crystallography_eq_position_eps_fractalspace = crystallography_eq_position_eps_realspace / 100;

	template<class T>
	struct Point {
	private:
		static constexpr T eq_pos = static_cast<T>(crystallography_eq_position_eps_fractalspace);
	public:
		using value_type = T;
		using array_type = ::std::array<value_type, 3>;

		struct Hash {
			std::size_t operator()(const std::array<T, 3>& point) const {
				if constexpr (std::is_floating_point_v<T>) {
					std::size_t hx = std::hash<T>{}(point[0]);
					std::size_t hy = std::hash<T>{}(point[1]);
					std::size_t hz = std::hash<T>{}(point[2]);

					return hx ^ (hy << 1) ^ (hz << 2) ^ (hx >> 31);
				}
				else {
					return 
						(point[0] * std::size_t(73856093)) ^ 
						(point[1] * std::size_t(19349663)) ^ 
						(point[2] * std::size_t(83492791));
				}
			}
		};


	public:
		array_type a = { 0,0,0 };

	public:
		// Constructors
		constexpr Point() noexcept = default;
		constexpr Point(value_type x, value_type y, value_type z) noexcept : a{ x, y, z } {};
		explicit constexpr Point(const array_type& other) noexcept : a(other) {};
		explicit constexpr Point(array_type&& other) noexcept : a(::std::move(other)) {};

		template <typename T2> 
			requires ((::std::integral<T2> || ::std::floating_point<T2>) && ::std::is_convertible<T2,T>::value)
		explicit constexpr Point(const Point<T2>& other) noexcept {
			a[0] = static_cast<T>(other[0]);
			a[1] = static_cast<T>(other[1]);
			a[2] = static_cast<T>(other[2]);
		}

		constexpr value_type r() const noexcept {
			return sqrt(fma(a[0], a[0], fma(a[1], a[1], a[2] * a[2])));
		}
		constexpr Point& MoveToCell() noexcept {
			a[0] -= ::std::floor(a[0]);
			a[1] -= ::std::floor(a[1]);
			a[2] -= ::std::floor(a[2]);
			return *this;
		}

		// Static constexpr functions
		static constexpr value_type Scalar(const Point& left, const Point& right) noexcept {
			return (left.a[0] * right.a[0] + left.a[1] * right.a[1] + left.a[2] * right.a[2]);
		}
		static constexpr Point Vector(const Point& left, const Point& right) noexcept {
			return Point(left.a[1] * right.a[2] - left.a[2] * right.a[1], left.a[2] * right.a[0] - left.a[0] * right.a[2], left.a[0] * right.a[1] - left.a[1] * right.a[0]);
		}
		static constexpr value_type distance(const Point& a, const Point& b) noexcept {
			value_type d0 = a.a[0] - b.a[0];
			value_type d1 = a.a[1] - b.a[1];
			value_type d2 = a.a[2] - b.a[2];
			return sqrt(fma(d0, d0, fma(d1, d1, d2 * d2)));
		}
		static constexpr value_type angleRad(const Point& a, const Point& b, const Point& c) noexcept {
			auto ab = distance(a, b);
			auto ac = distance(a, c);
			auto bc = distance(b, c);
			return ::std::acos(fma(ab, ab, fma(bc, bc, -ac * ac)) / (ab * bc * 2));
		}
		static constexpr value_type angleGrad(const Point& a, const Point& b, const Point& c) noexcept {
			return RadtoGrad(angleRad(a, b, c));
		}

		static constexpr value_type torsionRad(const Point& a, const Point& b, const Point& c, const Point& d) noexcept {
			auto b1 = c - b;
			auto b0 = a - b;
			b1 = b1 / (b1.r());
			auto b2 = d - c;
			auto v = b0 - b1 * (b0.a[0] * b1.a[0] + b0.a[1] * b1.a[1] + b0.a[2] * b1.a[2]);
			auto w = b2 - b1 * (b2.a[0] * b1.a[0] + b2.a[1] * b1.a[1] + b2.a[2] * b1.a[2]);
			auto x = v.a[0] * w.a[0] + v.a[1] * w.a[1] + v.a[2] * w.a[2];
			auto y = (b1.a[1] * v.a[2] - b1.a[2] * v.a[1]) * w.a[0] +
				(b1.a[2] * v.a[0] - b1.a[0] * v.a[2]) * w.a[1]
				+ (b1.a[0] * v.a[1] - b1.a[1] * v.a[0]) * w.a[2];
			return ::std::atan2(y, x);
		}
		static constexpr value_type torsionGrad(const Point& a, const Point& b, const Point& c, const Point& d) noexcept {
			return RadtoGrad(torsionRad(a, b, c, d));
		}
		constexpr Point round() const {
			return Point(std::round(a[0]), std::round(a[1]), std::round(a[2]));
		}
		constexpr Point floor() const {
			return Point(std::floor(a[0]), std::floor(a[1]), std::floor(a[2]));
		}

		constexpr value_type operator[](const uint8_t i) const noexcept { return a[i]; }
		constexpr value_type& operator[](const uint8_t i) noexcept { return a[i]; }

		// Operators
		constexpr Point operator-() const noexcept {
			return Point(-a[0], -a[1], -a[2]);
		}
		template<class OT> friend constexpr auto operator+(const Point& left, const Point<OT>& right) {
			using RT = typename std::conditional_t<std::is_same_v<T, OT>, T, decltype(left[0] + right[0])>;
			return Point<RT>(left.a[0] + right.a[0],
							 left.a[1] + right.a[1],
							 left.a[2] + right.a[2]);
		}

		template<class OT> friend constexpr auto operator+(const Point<T>& left, const OT b) noexcept {
			using RT = typename std::conditional_t<std::is_same_v<T, OT>, T, decltype(left.a[0] + b)>;
			return Point<RT>(left.a[0] + b, left.a[1] + b, left.a[2] + b);
		}
		template<class OT> friend constexpr auto operator-(const Point<T>& left, const Point<OT>& right) noexcept {
			using RT = typename std::conditional_t<std::is_same_v<T, OT>, T, decltype(left.a[0] - right.a[0])>;
			return Point<RT>(left.a[0] - right.a[0], left.a[1] - right.a[1], left.a[2] - right.a[2]);
		}
		template<class OT> friend constexpr auto operator-(const Point<T>& left, const OT b) noexcept {
			using RT = typename std::conditional_t<std::is_same_v<T, OT>, T, decltype(left.a[0] - b)>;
			return Point<RT>(left.a[0] - b, left.a[1] - b, left.a[2] - b);
		}

		friend constexpr Point operator*(const Point& left, const Point& right) noexcept {
			return Point(left.a[0] * right.a[0], left.a[1] * right.a[1], left.a[2] * right.a[2]);
		}
		friend constexpr Point operator*(const Point& left, const value_type b) noexcept {
			return Point(left.a[0] * b, left.a[1] * b, left.a[2] * b);
		}
		friend constexpr Point operator/(const Point& left, const Point& right) noexcept {
			return Point(left.a[0] / right.a[0],
						 left.a[1] / right.a[1],
						 left.a[2] / right.a[2]);
		}

		friend constexpr Point operator/(const Point& left, const value_type b) noexcept {
			return Point(left.a[0] / b, left.a[1] / b, left.a[2] / b);
		}

		template <class OT> inline Point& operator+=(const Point<OT>& right) noexcept {
			a[0] += static_cast<T>(right.a[0]);
			a[1] += static_cast<T>(right.a[1]);
			a[2] += static_cast<T>(right.a[2]);
			return *this;
		}
		inline Point& operator+=(const value_type right) noexcept {
			a[0] += right;
			a[1] += right;
			a[2] += right;
			return *this;
		}
		template <class OT> inline Point<OT>& operator-=(const Point<OT>& right) noexcept {
			a[0] += static_cast<T>(right.a[0]);
			a[1] += static_cast<T>(right.a[1]);
			a[2] += static_cast<T>(right.a[2]);
			return *this;
		}
		inline Point& operator-=(const value_type right) noexcept {
			a[0] -= right;
			a[1] -= right;
			a[2] -= right;
			return *this;
		}
		inline Point& operator*=(const value_type right) noexcept {
			a[0] *= right;
			a[1] *= right;
			a[2] *= right;
			return *this;
		}
		inline Point& operator/=(const value_type right) noexcept {
			a[0] /= right;
			a[1] /= right;
			a[2] /= right;
			return *this;
		}
		constexpr bool operator==(const Point& other) const noexcept {
			if (std::is_floating_point_v<T>)
				return abs(a[0] - other.a[0]) < T(0.00001) && abs(a[1] - other.a[1]) < T(0.00001) && abs(a[2] - other.a[2]) < T(0.00001);
			else
				return a[0] == other.a[0] && a[1] == other.a[1] && a[2] == other.a[2];
		}

		constexpr auto operator<=>(const Point& other) const noexcept = default;
	};

	template<class T>
	class Matrix {
	public:
		using size_t = unsigned char;
		using value_type = T;
		using array_type = ::std::array< ::std::array<T, 3>, 3>;
		using const_array_type = const array_type;
	private:
		array_type A{ { {0,0,0},{0,0,0},{0,0,0} } };
		template<class T2> friend class Matrix; // for constructors
	public:
		constexpr Matrix() noexcept = default;

		template<class T2> explicit constexpr Matrix(const Matrix<T2>& r) noexcept {
			A[0][0] = static_cast<T>(r.A[0][0]);
			A[0][1] = static_cast<T>(r.A[0][1]);
			A[0][2] = static_cast<T>(r.A[0][2]);
			A[1][0] = static_cast<T>(r.A[1][0]);
			A[1][1] = static_cast<T>(r.A[1][1]);
			A[1][2] = static_cast<T>(r.A[1][2]);
			A[2][0] = static_cast<T>(r.A[2][0]);
			A[2][1] = static_cast<T>(r.A[2][1]);
			A[2][2] = static_cast<T>(r.A[2][2]);
		}
		template<class T2> explicit constexpr Matrix(Matrix<T2>&& r) noexcept {

			A[0][0] = static_cast<T&&>(r.A[0][0]);
			A[0][1] = static_cast<T&&>(r.A[0][1]);
			A[0][2] = static_cast<T&&>(r.A[0][2]);
			A[1][0] = static_cast<T&&>(r.A[1][0]);
			A[1][1] = static_cast<T&&>(r.A[1][1]);
			A[1][2] = static_cast<T&&>(r.A[1][2]);
			A[2][0] = static_cast<T&&>(r.A[2][0]);
			A[2][1] = static_cast<T&&>(r.A[2][1]);
			A[2][2] = static_cast<T&&>(r.A[2][2]);
		}
		explicit constexpr Matrix(const T v) noexcept : A{ { {v,0,0},{0,v,0},{0,0,v} } } {}
		explicit constexpr Matrix(const T** input_massive) noexcept {
			for (size_t i = 0; i < 3; i++) {
				for (size_t j = 0; j < 3; j++) {
					A[i][j] = input_massive[i][j];
				}
			}
		}
		explicit constexpr Matrix(const T* input_massive) noexcept {
			for (size_t i = 0, k = 0; i < 3; i++) {
				for (size_t j = 0; j < 3; j++, k++) {
					A[i][j] = input_massive[k];
				}
			}
		}
		explicit constexpr Matrix(const_array_type& in) noexcept : A(in) {}
		explicit constexpr Matrix(array_type&& in) noexcept : A(std::move(in)) {}
		[[nodiscard]] constexpr T& El(const size_t a, const size_t b) noexcept {
			return A[a][b];
		}
		[[nodiscard]] constexpr T El(const size_t a, const size_t b) const noexcept {
			return A[a][b];
		}
		template<class T2>
		[[nodiscard]] constexpr Matrix<decltype(T()* T2())> operator*(const Matrix<T2>& right) const noexcept {
			using resv = decltype(T()* T2());
			std::array<resv, 3> a1 = { A[0][0] * right.A[0][0] + A[0][1] * right.A[1][0] + A[0][2] * right.A[2][0],
										A[0][0] * right.A[0][1] + A[0][1] * right.A[1][1] + A[0][2] * right.A[2][1],
										A[0][0] * right.A[0][2] + A[0][1] * right.A[1][2] + A[0][2] * right.A[2][2] };
			std::array<resv, 3> a2 = { A[1][0] * right.A[0][0] + A[1][1] * right.A[1][0] + A[1][2] * right.A[2][0],
										A[1][0] * right.A[0][1] + A[1][1] * right.A[1][1] + A[1][2] * right.A[2][1],
										A[1][0] * right.A[0][2] + A[1][1] * right.A[1][2] + A[1][2] * right.A[2][2] };
			std::array<resv, 3> a3 = { A[2][0] * right.A[0][0] + A[2][1] * right.A[1][0] + A[2][2] * right.A[2][0],
										A[2][0] * right.A[0][1] + A[2][1] * right.A[1][1] + A[2][2] * right.A[2][1],
										A[2][0] * right.A[0][2] + A[2][1] * right.A[1][2] + A[2][2] * right.A[2][2] };
			std::array<std::array<resv, 3>, 3> b{ a1, a2, a3 };
			return Matrix<resv>(b);
		}
		template<class T2>
		[[nodiscard]] constexpr std::array<decltype(T()* T2()), 3> operator*(const std::array<T2, 3>& right) const {
			std::array<decltype(T()* T2()), 3> res;
			res[0] = A[0][0] * right[0] + A[0][1] * right[1] + A[0][2] * right[2];
			res[1] = A[1][0] * right[0] + A[1][1] * right[1] + A[1][2] * right[2];
			res[2] = A[2][0] * right[0] + A[2][1] * right[1] + A[2][2] * right[2];
			return res;
		}
		template<class T2>
		[[nodiscard]] constexpr Matrix<decltype(T() / T2())> operator/(const T2 right) const noexcept {
			std::array<T, 3> a1 = { A[0][0] / right, A[0][1] / right, A[0][2] / right };
			std::array<T, 3> a2 = { A[1][0] / right, A[1][1] / right, A[1][2] / right };
			std::array<T, 3> a3 = { A[2][0] / right, A[2][1] / right, A[2][2] / right };
			array_type b{ a1, a2, a3 };
			return Matrix(std::move(b));
		}
		[[nodiscard]] constexpr Matrix<T> Transponate() const noexcept {
			std::array<T, 3> a1 = { A[0][0],A[1][0],A[2][0] };
			std::array<T, 3> a2 = { A[0][1],A[1][1],A[2][1] };
			std::array<T, 3> a3 = { A[0][2],A[1][2],A[2][2] };
			array_type b{ a1, a2, a3 };
			return Matrix(b);
		}
		[[nodiscard]] constexpr Matrix<T> Invert() const {
			const T det = Det();
			std::array<T, 3> a1 = { (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det, (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det, (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det };
			std::array<T, 3> a2 = { (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det, (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det, (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / det };
			std::array<T, 3> a3 = { (A[1][0] * A[2][1] - A[2][0] * A[1][1]) / det, (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / det, (A[1][1] * A[0][0] - A[1][0] * A[0][1]) / det };
			array_type b = { a1,a2,a3 };
			return Matrix(b);
		}
		[[nodiscard]] constexpr Matrix<T> Modul() const noexcept {
			constexpr T zero = 0;
			std::array<T, 3> a1 = { A[0][0] < zero ? A[0][0] : -A[0][0], A[0][1] < zero ? A[0][1] : -A[0][1], A[0][2] < zero ? A[0][2] : -A[0][2] };
			std::array<T, 3> a2 = { A[1][0] < zero ? A[1][0] : -A[1][0], A[1][1] < zero ? A[1][1] : -A[1][1], A[1][2] < zero ? A[1][2] : -A[1][2] };
			std::array<T, 3> a3 = { A[2][0] < zero ? A[2][0] : -A[2][0], A[2][1] < zero ? A[2][1] : -A[2][1], A[2][2] < zero ? A[2][2] : -A[2][2] };
			array_type b{ a1, a2, a3 };
			return Matrix(std::move(b));
		}
		constexpr double Trace() const noexcept {
			return (A[0][0] + A[1][1] + A[2][2]) / 3.0;
		}
		constexpr T Det() const noexcept {
			return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1]
				- A[0][2] * A[1][1] * A[2][0] - A[0][1] * A[1][0] * A[2][2] - A[0][0] * A[1][2] * A[2][1];
		}
		template<class T2>
		constexpr void MultMatrixByArray(const std::array<T2, 3>& sup) noexcept {
			return MultMatrixByArray(sup[0], sup[1], sup[2]);
		}
		template<class T2>
		constexpr void MultMatrixByArray(const T2 x, const T2 y, const T2 z) noexcept {
			A[0][0] *= x;
			A[1][0] *= x;
			A[2][0] *= x;
			A[0][1] *= y;
			A[1][1] *= y;
			A[2][1] *= y;
			A[0][2] *= z;
			A[1][2] *= z;
			A[2][2] *= z;
		}
		template<class T2>
		friend constexpr Point<decltype(T()* T2())> operator*(const Matrix<T>& left, const Point<T2>& right) noexcept {
			Point<decltype(T()* T2())> res;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					res[i] = std::fma(left.El(i, j), right[j], res[i]);
				}
			}
			return res;
		}

	};


	template<class T>
	constexpr Matrix<T> EqualMatrix(static_cast<T>(1));

	template<class T>
	struct Plane {
		using value_type = T; // the same as T
		std::array<T, 4> a = { T(0), T(0), T(0), T(0) }; // [ A, B, C, D]

		constexpr Plane() noexcept = default;
		constexpr Plane(const Point<T>& a1, const Point<T>& a2, const Point<T>& a3) noexcept {
			a[0] = a1[1] * (a2[2] - a3[2]) + a2[1] * (a3[2] - a1[2]) + a3[1] * (a1[2] - a2[2]);
			a[1] = a1[2] * (a2[0] - a3[0]) + a2[2] * (a3[0] - a1[0]) + a3[2] * (a1[0] - a2[0]);
			a[2] = a1[0] * (a2[1] - a3[1]) + a2[0] * (a3[1] - a1[1]) + a3[0] * (a1[1] - a2[1]);
			a[3] = -(a1[0] * (a2[1] * a3[2] - a3[1] * a2[2]) +
					 a2[0] * (a3[1] * a1[2] - a1[1] * a3[2]) +
					 a3[0] * (a1[1] * a2[2] - a2[1] * a1[2]));
		}

		// Construct a plane from a point and a normal
		constexpr Plane(const Point<T>& point, const Point<T>& normal) noexcept {
			// Set plane
			a[0] = normal[0];
			a[1] = normal[1];
			a[2] = normal[2];
			a[3] = -(normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2]);
		}

		constexpr Plane(const Plane& p, const Point<T>& a1) noexcept {
			a[0] = p.a[0];
			a[1] = p.a[1];
			a[2] = p.a[2];
			a[3] = -(a[0] * a1[0] + a[1] * a1[1] + a[2] * a1[2]);
		}

		Point<T> make_projection(const Point<T>& p) const {
			auto pcos = (a[0] * p[0] + a[1] * p[1] + a[2] * p[2] + a[3]) / (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
			return Point<T>(p[0] - a[0] * pcos, p[1] - a[1] * pcos, p[2] - a[2] * pcos);
		}
		T distance(const Point<T>& p) const {
			return abs(a[0] * p[0] + a[1] * p[1] + a[2] * p[2] + a[3]) / sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		}

		constexpr Point<T> normal() const noexcept {
			Point<T> norm(a[0], a[1], a[2]);
			if (T norm_length = norm.r(); norm_length > 1e-10)
				return norm / norm_length;
			else {
				return Point<T>(0, 0, 1);
			}
		}

		T side(const Point<T>& p) const {
			return a[0] * p[0] + a[1] * p[1] + a[2] * p[2] + a[3];
		}
	};

	template<class T>
	class Polygon {
	public:
		using PointType = Point<T>;
		using PlaneType = Plane<T>;
		static constexpr T limit = ::std::numeric_limits<T>::epsilon();
	private:
		std::vector<PointType> vertices_;
		PlaneType plane_;
		bool is_valid_;

	public:
		// Default Constructor
		Polygon() = default;

		// Construct from vector of vertices
		explicit Polygon(const std::vector<PointType>& verts) : vertices_(verts) {
			if (verts.size() >= 3) {
				plane_ = PlaneType(verts[0], verts[1], verts[2]);
				is_valid_ = true;
			}
			else {
				is_valid_ = false;
			}
		}

		Polygon(const PlaneType& plane, const PointType& center, T side_length)
			: plane_(plane), is_valid_(true) {
			// Create normal vector by hand (need to check initial length)
			PointType normal(plane.a[0], plane.a[1], plane.a[2]);

			if (T normal_length = normal.r(); normal_length > 1e-10) {
				normal = normal / normal_length;
			}
			else {
				// if normal is zero, the plane is XY plane
				normal = PointType(0, 0, 1);
				plane_ = PlaneType(PointType(0, 0, 0), PointType(1, 0, 0), PointType(0, 1, 0));
			}

			// Find first vector in plane
			PointType u;
			if (std::abs(normal[0]) > std::abs(normal[1])) {
				u = PointType(-normal[2], 0, normal[0]);
			}
			else {
				u = PointType(0, normal[2], -normal[1]);
			}
			u = u / u.r();

			// Find second vector in plane
			PointType v = PointType::Vector(normal, u);
			v = v / v.r();

			// Resize to half side length
			T half_side = side_length / 2;
			u = u * half_side;
			v = v * half_side;

			// Create vertices of square
			vertices_.reserve(4);
			vertices_.push_back(center + u + v);
			vertices_.push_back(center - u + v);
			vertices_.push_back(center - u - v);
			vertices_.push_back(center + u - v);
		}

		Polygon(const std::vector<PointType>& verts, const PlaneType& pl)
			: vertices_(verts), plane_(pl), is_valid_(true) {
			bool t = isConvex();
			_ASSERT(isConvex());
		}

		// Method to clip the polygon by a plane
		std::optional<std::pair<PointType,PointType>> clipByPlane(const PlaneType& clipping_plane) {
			if (!is_valid_ || vertices_.empty()) 
				return std::optional<std::pair<PointType, PointType>>();

			auto intersect = [](const PointType& a, const PointType& b, const PlaneType& plane)->PointType {
				PointType d = b - a;
				T denominator = d[0] * plane.a[0] + d[1] * plane.a[1] + d[2] * plane.a[2];
				_ASSERT(abs(denominator) > limit);
				T t = -plane.side(a) / denominator;
				return a + d * t;
				};
			auto vs = vertices_.size();

			auto [e1, e2] = findIntersectionPoints(clipping_plane);

			// Check: all points are on the same side of the plane?
			if (e1 == vs) {
				if (clipping_plane.side(vertices_[0]) < 0) {
					vertices_.clear();
					is_valid_ = false;
				}
				return std::optional<std::pair<PointType, PointType>>();
			}
			else 
			{
				// erase from e1 (including) to e2 (excluding)
				// find intersection points
				PointType inter1 = intersect(vertices_[(e1 + vs - 1) % vs], vertices_[e1], clipping_plane);
				PointType inter2 = intersect(vertices_[(e2 + vs - 1) % vs], vertices_[e2], clipping_plane);
				if (e1 > e2) {
					vertices_.erase(vertices_.begin() + e1, vertices_.end());
					vertices_.erase(vertices_.begin(), vertices_.begin() + e2);
					vertices_.push_back(inter1);
					vertices_.push_back(inter2);
				}
				else if (e2 > e1) {
					vertices_[e1] = inter1;
					vertices_[e2] = inter2;
					if ((e2 - e1) > 2) {
						vertices_.erase(vertices_.begin() + e1 + 1, vertices_.begin() + e2 - 1);
					}
				}
				else if (e1 == e2) {
					// overwrite e1 and add one more point
					vertices_[e1] = inter1;
					vertices_.insert(vertices_.begin() + e1 + 1, inter2);
				}
				_ASSERT(isConvex());
				return std::make_optional(std::pair<PointType, PointType>(inter1, inter2));
			}


			return std::optional<std::pair<PointType, PointType>>();
		}

		constexpr const PointType& operator[](size_t i) const noexcept {
			return vertices_[i];
		}
		constexpr size_t size() const noexcept {
			return vertices_.size();
		}

		bool isConvex() const {
			if (vertices_.size() < 3) return false;
			auto normal = plane_.normal();
			for (size_t i = 0; i < vertices_.size(); i++) {
				PointType current = vertices_[i];
				PointType next = vertices_[(i + 1) % vertices_.size()];
				PointType nextNext = vertices_[(i + 2) % vertices_.size()];

				PointType edge1 = next - current;
				PointType edge2 = nextNext - next;

				PointType cross = PointType::Vector(edge1, edge2);
				if (PointType::Scalar(cross, normal) < -limit) {
					return false;
				}
			}
			return true;
		}
		void fixVertexOrder() {
			if (vertices_.size() < 3) return;

			// Calculate the center of the polygon 
			Point<T> center(0, 0, 0);
			for (const auto& vertex : vertices_) {
				center += vertex;
			}
			center = center / static_cast<T>(vertices_.size());

			// Calculate the normal vector of the plane
			Point<T> normal = plane_.normal();

			// Sort the vertices in counter-clockwise order around the center
			std::ranges::sort(vertices_,
							  [&](const Point<T>& a, const Point<T>& b) {
								  Point<T> vecA = a - center;
								  Point<T> vecB = b - center;

								  // ¬ычисл€ем угол через векторное произведение
								  Point<T> cross = Point<T>::Vector(vecA, vecB);
								  T dot_with_normal = Point<T>::Scalar(cross, normal);

								  return dot_with_normal > 0;
							  });
		}

		bool isValid() const {
			return is_valid_ && vertices_.size() >= 3;
		}

		const std::vector<PointType>& getVertices() const {
			return vertices_;
		}
	private:
		std::pair<size_t,size_t> findIntersectionPoints(const PlaneType& clipping_plane) const {
			size_t e1 = vertices_.size();
			size_t e2 = vertices_.size();
			size_t vs = vertices_.size();

			PointType inter1;
			PointType inter2;

			// Prepare dist to plane vector
			::std::vector<T> dists(vs, 0);
			for (size_t iter = 0; iter < vs; iter++)
			{
				dists[iter] = clipping_plane.side(vertices_[iter]);
			}

			if (dists.front() < 0 && dists.back() > 0)
			{
				e1 = 0;
			}

			if (dists.front() > 0 && dists.back() < 0)
			{
				e2 = 0;
			}

			for (size_t i = 1; i < vs; i++)
			{
				if (dists[i - 1] > 0 && dists[i] < 0)
				{
					e1 = i;
				}

				if (dists[i - 1] < 0 && dists[i] > 0)
				{
					e2 = i;
				}
			}

			return std::make_pair(e1, e2);
		}
	};

	template<class T>
	class VoronoiCell {
	public:
		using PointType = Point<T>;
		using Face = Polygon<T>;
		using PlaneType = typename Face::PlaneType;
		using FaceVector = ::std::vector<Face>;

	private:
		FaceVector faces_;
		PointType seed_ = PointType(0, 0, 0);

	public:
		// Default constructor creates cube around [0,0,0]
		constexpr VoronoiCell() noexcept {
			initiate_cube_faces_on_seed();
		}
		// Creates cube around seed
		constexpr explicit VoronoiCell(const PointType& seed, bool init = true) noexcept : seed_(seed) {
			if (init) { initiate_cube_faces_on_seed(); }
		}

		/// <summary>
		/// Cut both VoronoiCells by each other
		/// </summary>
		/// <returns> 0 - if correct even if any of the cells are empty, 
		///           1 - if Cells are too close</returns>
		static int interact(VoronoiCell& a, VoronoiCell& b) {
			// Define if faces are not empty
			bool empty_a = a.faces_.empty();
			bool empty_b = b.faces_.empty();

			if (empty_a && empty_b)
				return 0; // Both cells are empty

			PointType d = b.seed_ - a.seed_;

			// Move interval "d" to [-0.5; 0.5]
			for (int i = 0; i < 3; ++i) {
				if (d[i] > 0.5) d[i] -= 1.0;
				else if (d[i] < -0.5) d[i] += 1.0;
			}

			// If points too close - stop
			T d_length = d.r();
			if (d_length < 1e-6)
				return 1;

			PointType normal = d / d_length;
			PointType half_d = d * 0.5;

			// Create plane and clip a by it
			if (!empty_a) {
				PointType midpoint_a = a.seed_ + half_d;
				PlaneType plane(midpoint_a, -normal);
				a.clipByPlaneAndAddNewFace(plane);
			}

			// Using inverted plane for b
			if (!empty_b) {
				PointType midpoint_b = b.seed_ - half_d;
				PlaneType inverted_plane(midpoint_b, normal);
				b.clipByPlaneAndAddNewFace(inverted_plane);
			}

			return 0;
		}
		constexpr const PointType& getSeed() const noexcept { return seed_; }
		constexpr const FaceVector& getFaces() const noexcept { return faces_; }

	private:
		void clipByPlaneAndAddNewFace(const PlaneType& clipping_plane) {
			std::vector<PointType> intersection_points;
			std::vector<Face> new_faces;

			// Collect all intersection points and create new faces
			for (auto& face : faces_) {
				// Save original vertices for intersection detection
				std::vector<PointType> original_vertices;
				for (size_t i = 0; i < face.size(); ++i) {
					original_vertices.push_back(face[i]);
				}

				// Clip the face
				face.clipByPlane(clipping_plane);

				// If face remains valid, add it
				if (face.size() >= 3) {
					new_faces.push_back(face);

					// Collect intersection points with this face
					collectIntersectionPoints(original_vertices, clipping_plane, intersection_points);
				}
				// If face becomes invalid (less than 3 vertices), it is not added - thus removed
			}

			// Replace old faces with new ones
			faces_ = std::move(new_faces);

			// Create new face from intersection points if there are enough
			if (intersection_points.size() >= 3) {
				createNewFaceFromIntersections(intersection_points, clipping_plane);
			}
		}

		void collectIntersectionPoints(const std::vector<PointType>& vertices,
									   const PlaneType& plane,
									   std::vector<PointType>& intersection_points) {
			const size_t n = vertices.size();
			if (n < 3) return;

			PointType prev_vertex = vertices.back();
			T prev_dist = plane.side(prev_vertex);

			for (size_t i = 0; i < n; i++) {
				const PointType& current_vertex = vertices[i];
				const T current_dist = plane.side(current_vertex);

				// If edge intersects the plane, find intersection point
				if (prev_dist * current_dist < 0) {
					const T t = prev_dist / (prev_dist - current_dist);
					const PointType intersection = prev_vertex + (current_vertex - prev_vertex) * t;
					intersection_points.push_back(intersection);
				}

				prev_vertex = current_vertex;
				prev_dist = current_dist;
			}
		}

		void createNewFaceFromIntersections(std::vector<PointType>& points,
											const PlaneType& plane) {
			if (points.size() < 3) return;

			// Order points in correct order (counter-clockwise relative to normal)
			orderPointsOnPlane(points, plane);
			// Remove duplicates
			removeDuplicatePoints(points);

			// Create new face
			Face new_face(points, plane);
			if (new_face.isConvex()) {
				faces_.push_back(new_face);
			}
		}

		void removeDuplicatePoints(std::vector<PointType>& points) const {
			const auto s = points.size();

			// Remove consecutive duplicates
			for (size_t i = s - 1; i > 0; --i)
			{
				if (PointType::distance(points[i], points[i - 1]) < static_cast<T>(1e-6)) {
					points.erase(points.begin() + i);
				}
			}

			// Check first and last element
			if (PointType::distance(points.front(), points.back()) < static_cast<T>(1e-6)) {
				points.pop_back();
			}
		}

		void orderPointsOnPlane(std::vector<PointType>& points, const PlaneType& plane) {
			if (points.size() < 3) return;

			// Find center of mass of points
			PointType center(0, 0, 0);
			for (const auto& p : points) {
				center += p;
			}
			center = center / static_cast<T>(points.size());

			// Get plane normal
			PointType normal = plane.normal();

			// Choose arbitrary vector in plane (perpendicular to normal)
			PointType reference_vector;
			if (std::abs(normal[0]) > std::abs(normal[1])) {
				reference_vector = PointType(-normal[2], 0, normal[0]);
			}
			else {
				reference_vector = PointType(0, normal[2], -normal[1]);
			}
			reference_vector = reference_vector / reference_vector.r();

			// Sort points by angle relative to center
			std::ranges::sort(points,
							  [&](const PointType& a, const PointType& b)
							  {
								  PointType vecA = a - center;
								  PointType vecB = b - center;

								  // Project onto plane
								  PointType projA = vecA - normal * PointType::Scalar(vecA, normal);
								  PointType projB = vecB - normal * PointType::Scalar(vecB, normal);

								  if (projA.r() < 1e-10 || projB.r() < 1e-10) {
									  return false; // Points too close to center
								  }

								  // Normalize
								  projA = projA / projA.r();
								  projB = projB / projB.r();

								  // Calculate angles using scalar and vector products
								  T cosA = PointType::Scalar(reference_vector, projA);
								  T sinA = PointType::Scalar(PointType::Vector(reference_vector, projA), normal);
								  T angleA = std::atan2(sinA, cosA);

								  T cosB = PointType::Scalar(reference_vector, projB);
								  T sinB = PointType::Scalar(PointType::Vector(reference_vector, projB), normal);
								  T angleB = std::atan2(sinB, cosB);

								  return angleA < angleB;
							  });
		}


		inline void clipByPlane(const PlaneType& clipping_plane) {
			for (auto& face : faces_) {
				face.clipByPlane(clipping_plane);
			}
		}
		constexpr void initiate_cube_faces_on_seed() {
			std::array<PointType, 8> cube;
			for (int i = 0; i < 8; ++i) {
				cube[i] = base_vertices[i] + seed_;
			}
			faces_.reserve(6);
			for (int i = 0; i < 6; ++i) {
				faces_.emplace_back(Face({ cube[face_indices[i][0]], 
										   cube[face_indices[i][1]], 
										   cube[face_indices[i][2]], 
										   cube[face_indices[i][3]] }));
			}
		}

		// Base array of vertices for a cube
		static constexpr std::array<PointType, 8> base_vertices = { {
			PointType{-0.5, -0.5, -0.5}, // 0
			PointType{ 0.5, -0.5, -0.5}, // 1
			PointType{ 0.5,  0.5, -0.5}, // 2
			PointType{-0.5,  0.5, -0.5}, // 3
			PointType{-0.5, -0.5,  0.5}, // 4
			PointType{ 0.5, -0.5,  0.5}, // 5
			PointType{ 0.5,  0.5,  0.5}, // 6
			PointType{-0.5,  0.5,  0.5}  // 7
		} };

		// Indexes for each face of the cube (conter-clockwise from outside)
		static constexpr std::array<std::array<int, 4>, 6> face_indices = { {
			{4, 7, 6, 5}, // front face
			{0, 1, 2, 3}, // back face
			{0, 3, 7, 4}, // left face
			{1, 5, 6, 2}, // right face
			{0, 4, 5, 1}, // bottom face
			{3, 2, 6, 7}  // top face
		} };
	};

	template<class T, class AI> class HashedSpace;

	template<class T>
	class VoronoiDiagram {
	public:
		using PointType = Point<T>;
		using VoronCell = VoronoiCell<T>;

		template <class AI>
		using BondList = ::std::vector<::std::pair<AI, AI>>;
		using PointVector = ::std::vector<PointType>;
		using CellVector = ::std::vector<VoronCell>;
		using BoolVector = ::std::vector<bool>;

		enum class State : unsigned char {
			Uninitialized = 0,
			Cubic_cells = 1,
			Correct_cells = 2
		};
	private:
		//Data
		BoolVector flags_;
		CellVector cells_;
		State state = State::Uninitialized;
	public:
		constexpr VoronoiDiagram() noexcept = default;
		constexpr explicit VoronoiDiagram(const PointVector& points, const BoolVector& flags = BoolVector()) noexcept : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points.size(), true);
			}
			addPoints(points, flags_);
		}
		
		template<class AI>
		constexpr VoronoiDiagram(const PointVector& points, const BondList<AI>& bonds, const BoolVector& flags = BoolVector(true, points.size())) noexcept : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points.size(), true);
			}
			addPoints(points, flags_);
			calculateFaces<AI>(bonds);
		}

		void addPoints(const PointVector& points, const BoolVector& flags) noexcept {
			cells_.reserve(points.size());
			for (size_t i = 0; i < points.size(); i++)
			{
				cells_.emplace_back(points[i], flags[i]);
			}
			state = State::Cubic_cells;
		}

		template <class AI>
		int calculateFaces(const BondList<AI>& bonds) noexcept {
			if (state == State::Uninitialized) 
				return 1; // Error: VoronoiDiagram not initialized
			for (auto& bond : bonds) {
				auto interaction_result = VoronCell::interact(cells_[bond.first], cells_[bond.second]);
				if (interaction_result != 0) 
					return 2; // Error: Cells too close
			}
			state = State::Correct_cells;
			return 0;
		}

		constexpr T calculateLongestDiagonal(const Matrix<T>& mat) const noexcept {
			T ret = 0;
			for (auto& vcell : cells_) {
				const PointType seed = vcell.getSeed();
				for (const auto& face : vcell.getFaces()) {
					for (size_t i = 0; i < face.size(); i++)
					{
						T val = (mat*(face[i] - seed)).r();
						if (val > ret) ret = val;
					}
				}
			}
			return ret*2;
		}

		CellVector extractCells() noexcept {
			if (state == State::Uninitialized) {
				return {};
			}
			state = State::Uninitialized;
            return std::move(cells_);
		}
	};

	template<class T>
	struct Cell {
	public:
		using value_type = T;
		using PointType = Point<T>;
		using array_type = ::std::array<T, 3>;
		using matrix_type = geometry::Matrix<T>;
	private:
		array_type lattice_;
		array_type angleRad_;
		array_type angleGrad_;
		matrix_type fracToCart_;
		matrix_type cartToFrac_;
		constexpr void createFromFracToCart(const matrix_type& Mat) noexcept {
			fracToCart_ = Mat;
			cartToFrac_ = fracToCart_.Invert();
			takeCellFromFracToCart();
		}
		constexpr void takeCellFromFracToCart() noexcept {
			lattice_[0] = fracToCart_.El(0, 0);
			lattice_[1] = std::sqrt((fracToCart_.El(1, 1) * fracToCart_.El(1, 1)) + (fracToCart_.El(0, 1) * fracToCart_.El(0, 1)));
			value_type CosG = fracToCart_.El(0, 1) / lattice_[1]; //CosG
			value_type SinG = std::sqrt(1 - (CosG * CosG)); //SinG
			angleRad_[2] = std::acos(CosG);
			value_type tem2 = fracToCart_.El(0, 2);
			value_type temp = tem2 * CosG + fracToCart_.El(1, 2) * SinG;
			value_type tem1 = fracToCart_.El(2, 2) * SinG;
			lattice_[2] = std::sqrt((tem1 * tem1) + (temp * temp) - 2 * CosG * temp * tem2 + (tem2 * tem2)) / SinG;
			angleRad_[1] = std::acos(fracToCart_.El(0, 2) / lattice_[2]);
			angleRad_[0] = std::acos(temp / lattice_[2]);
			for (int i = 0; i < 3; i++)
				angleGrad_[i] = RadtoGrad(angleRad_[i]);
		}
		constexpr void createFromCartToFrac(const matrix_type& Mat) noexcept {
			cartToFrac_ = Mat;
			fracToCart_ = cartToFrac_.Invert();
			takeCellFromFracToCart();
		}
		constexpr void createMatrix() noexcept {
			value_type cos0 = std::cos(angleRad_[0]);
			value_type cos1 = std::cos(angleRad_[1]);
			value_type cos2 = std::cos(angleRad_[2]);
			value_type sin2 = std::sin(angleRad_[2]);
			fracToCart_.El(0, 0) = lattice_[0];
			fracToCart_.El(0, 1) = lattice_[1] * cos2;
			fracToCart_.El(1, 1) = lattice_[1] * sin2;
			fracToCart_.El(0, 2) = lattice_[2] * cos1;
			fracToCart_.El(1, 2) = lattice_[2] * (cos0 - cos1 * cos2) / sin2;
			fracToCart_.El(2, 2) = lattice_[2] * std::sqrt(sin2 * sin2 - cos0 * cos0 - cos1 * cos1 + 2 * cos0 * cos1 * cos2) / sin2;
			cartToFrac_ = fracToCart_.Invert();
		}
	public:
		explicit constexpr Cell(const matrix_type& Mat, const bool is_FracToCart = true) {
			create(Mat, is_FracToCart);
		}
		explicit constexpr Cell(const value_type a = 10, const value_type b = 10, const value_type c = 10, const value_type alpha = 90, const value_type beta = 90, const value_type gamma = 90, const bool is_grad = true) {
			create(a, b, c, alpha, beta, gamma, is_grad);
		}
		explicit constexpr Cell(const std::array<T, 6>& ar, const bool is_grad = true) {
			create(ar[0], ar[1], ar[2], ar[3], ar[4], ar[5], is_grad);
		}
		constexpr Cell(Cell&&) noexcept = default;
		constexpr Cell(const Cell&) noexcept = default;
		constexpr Cell& operator= (const Cell&) noexcept = default;
		constexpr Cell& operator= (Cell&&) noexcept = default;

		[[nodiscard]] constexpr value_type lat_dir(const unsigned char i) const noexcept {
			return lattice_[i];
		}
		[[nodiscard]] constexpr value_type& lat_dir(const unsigned char i) noexcept {
			return lattice_[i];
		}
		constexpr void create(const matrix_type& Mat, const bool is_FracToCart = true) noexcept {
			if (is_FracToCart) createFromFracToCart(Mat);
			else createFromCartToFrac(Mat);
		}
		constexpr void create(const value_type a = 10, const value_type b = 10, const value_type c = 10, const value_type alpha = 90, const value_type beta = 90, const value_type gamma = 90, const bool is_grad = true) {
			lattice_[0] = a;
			lattice_[1] = b;
			lattice_[2] = c;
			if (is_grad) {
				angleGrad_[0] = alpha;
				angleGrad_[1] = beta;
				angleGrad_[2] = gamma;
				for (int i = 0; i < 3; i++) {
					angleRad_[i] = GradtoRad(angleGrad_[i]);
				}
			}
			else {
				angleRad_[0] = alpha;
				angleRad_[1] = beta;
				angleRad_[2] = gamma;
				for (int i = 0; i < 3; i++) {
					angleGrad_[i] = RadtoGrad(angleRad_[i]);
				}
			}
			createMatrix();
		}

		[[nodiscard]] constexpr const value_type& getAngleRad(const unsigned char i) const noexcept {
			return angleRad_[i];
		}
		[[nodiscard]] constexpr value_type getAngleRad(const unsigned char i) noexcept {
			return angleRad_[i];
		}
		[[nodiscard]] constexpr const value_type& getAngleGrad(const unsigned char i) const noexcept {
			return angleGrad_[i];
		}
		[[nodiscard]] constexpr value_type getAngleGrad(const unsigned char i) noexcept {
			return angleGrad_[i];
		}
		[[nodiscard]] constexpr const matrix_type& fracToCart() const noexcept {
			return fracToCart_;
		}
		[[nodiscard]] constexpr const matrix_type& cartToFrac() const noexcept {
			return cartToFrac_;
		}

		value_type distance_in_01(const PointType& a, const PointType& b) const {
			PointType d = (a - b).MoveToCell();
			if (d[0] > 0.5) d[0] = 1 - d[0];
			if (d[1] > 0.5) d[1] = 1 - d[1];
			if (d[2] > 0.5) d[2] = 1 - d[2];
			return (fracToCart_*d).r();
		}

		template <class I>
		[[nodiscard]] constexpr geometry::Point<I> findOptimalSupercell(const value_type cutoff, I minimum) const noexcept {
			constexpr char shortContactsInTriclinic = 10;

			constexpr std::array<geometry::Point<I>, shortContactsInTriclinic> directions{
				geometry::Point<I>(1, 0, 1),
				geometry::Point<I>(0, 1,-1),
				geometry::Point<I>(1, 0,-1),
				geometry::Point<I>(1, 1, 0),
				geometry::Point<I>(1,-1, 0),
				geometry::Point<I>(0, 1, 1),
				geometry::Point<I>(1, 1, 1),
				geometry::Point<I>(1, 1,-1),
				geometry::Point<I>(1,-1, 1),
				geometry::Point<I>(1,-1,-1) };

			geometry::Point<I> currentSuperCell(
				std::max(static_cast<I>(std::ceil(cutoff / lattice_[0])), minimum),
				std::max(static_cast<I>(std::ceil(cutoff / lattice_[1])), minimum),
				std::max(static_cast<I>(std::ceil(cutoff / lattice_[2])), minimum));

			char i = 0;
			while (i < shortContactsInTriclinic) {
				if ((fracToCart_ * (directions[i] * currentSuperCell)).r() > cutoff) {
					i++;
					continue;
				}

				unsigned char minDimention = 0;
				switch (i) {
				case 0:
				case 1:
					if (lattice_[1] * currentSuperCell[1] > lattice_[2] * currentSuperCell[2]) minDimention = 2;
					else minDimention = 1;
					break;
				case 2:
				case 3:
					if (lattice_[0] * currentSuperCell[0] > lattice_[2] * currentSuperCell[2]) minDimention = 2;
					else minDimention = 0;
					break;
				case 4:
				case 5:
					if (lattice_[0] * currentSuperCell[0] > lattice_[1] * currentSuperCell[1]) minDimention = 1;
					else minDimention = 0;
					break;
				default:
					if (lattice_[minDimention] * currentSuperCell[minDimention] > lattice_[1] * currentSuperCell[1]) minDimention = 1;
					if (lattice_[minDimention] * currentSuperCell[minDimention] > lattice_[2] * currentSuperCell[2]) minDimention = 2;
					break;
				}
				currentSuperCell[minDimention] = currentSuperCell[minDimention] + 1;
			}
			return currentSuperCell;
		}
	};

	template<class T>
	struct Symm
	{
		using matrix_t = geometry::Matrix<T>;
		using point_t = geometry::Point<T>;
		matrix_t mat;
		point_t point;
		Symm() = default;
		Symm MirrorSymm() const {
			Symm out;
			out.mat = mat.Invert();
			out.point = -point;
			return out;
		}
		explicit Symm(const char* str) noexcept : mat(), point(static_cast<T>(0.0), static_cast<T>(0.0), static_cast<T>(0.0))
		{
			// separate into 3 domains
			std::array<size_t, 3> bracket;
			bracket[0] = findcomma(str);
			bracket[1] = findcomma(str + bracket[0] + 1);
			bracket[2] = findcomma(str + bracket[0] + 1 + bracket[1] + 1);

			for (unsigned char i = 0; i < 3; i++)
			{
				auto ret = parse(str, bracket[i]);
				for (unsigned char j = 0; j < 3; j++)
				{
					mat.El(i, j) = ret.first[j];
				}
				point[i] = ret.second;
				str += (bracket[i] + 1);
			}
		}
		inline point_t GenSymm(const point_t& in) const
		{
			return point + (mat * in);
		}
		inline point_t GenSymmNorm(const point_t& in) const
		{
			point_t res = GenSymm(in);
			res.MoveToCell();
			return res;
		}
		inline bool is_Eq() const noexcept
		{
			return mat.El(0, 0) == 1 && mat.El(1, 0) == 0 && mat.El(2, 0) == 0 &&
				mat.El(0, 1) == 0 && mat.El(1, 1) == 1 && mat.El(2, 1) == 0 &&
				mat.El(0, 2) == 0 && mat.El(1, 2) == 0 && mat.El(2, 2) == 1 &&
				point[0] == 0 && point[1] == 0 && point[2] == 0;
		}
	private:
		size_t findcomma(const char* str) const {
			size_t n = 0;
			for (; str[n] != ',' && str[n] != '\0'; n++);
			return n;
		}
		std::pair<point_t, T> parse(const char* str, const size_t len) const {
			point_t p{ 0,0,0 };
			T shift = 0;
			bool minus = false;
			for (unsigned int i = 0; i < len; i++) // iterator "i" modifies in parseshift function
			{
				switch (str[i]) {
				case 'x': [[fallthrough]];
				case 'X':
					if (minus == true) p[0] = -1;
					else p[0] = 1;
					minus = false;
					break;
				case 'y': [[fallthrough]];
				case 'Y':
					if (minus == true) p[1] = -1;
					else p[1] = 1;
					minus = false;
					break;
				case 'z': [[fallthrough]];
				case 'Z':
					if (minus == true) p[2] = -1;
					else p[2] = 1;
					minus = false;
					break;
				case ' ': [[fallthrough]];
				case '\'': [[fallthrough]];
				case '\"': [[fallthrough]];
				case '+':
					break;
				case '-':
					minus = !minus;
					break;
				default:

					T partshift = parseshift(str, i, len); // Modifies "i"
					if (minus) shift -= partshift;
					else shift += partshift;
					minus = false;
					break;
				}
			}
			return std::make_pair(p, shift);

		}
		T parseshift(const char* str, unsigned int& iter, const size_t len) const {
			bool dot = false;
			bool slash = false;

			int upper = 0;
			int lower = 1;
			for (; iter < len; iter++)
			{
				if (str[iter] < '0' || str[iter] > '9')
					switch (str[iter]) {
					case '.':
						dot = true;
						break;
					case '+':
					case '-':
						iter--;
						return upper / static_cast<T>(lower);
					case '/':
						slash = true;
						iter++;
						lower = (str[iter] - '0');
						break;
					default:
						// unexpected symbol
						break;
					}
				else {
					int num = (str[iter] - '0');
					if (dot) {
						upper = upper * 10 + num;
						lower *= 10;
					}
					else if (slash) {
						lower = lower * 10 + num;
					}
					else {
						upper = upper * 10 + num;
					}
				}
			}
			iter--;
			return upper / static_cast<T>(lower);
		}

	};

	template<class T, class AI>
	class HashedSpace {
	public:
		using AtomIndex = AI;
		using FloatingPointType = T;
		using PointType = Point<T>;
		using CellType = Cell<T>;
		using DimentionType = unsigned char;

		using SupListType = ::std::list<AtomIndex>;
		using SupType = ::std::vector<::std::vector<::std::vector<SupListType>>>;
		using SupPoint = geometry::Point<size_t>;

		/// <summary>
		/// Constructor
		/// </summary>
		/// <param name="cell"> - Unit cell</param>
		/// <param name="maxbond"> - calculated in Distances</param>
		constexpr HashedSpace(const CellType& cell, FloatingPointType maxbond) noexcept :
			cell_(cell) {
			calculateSep(maxbond);
		}

		/// <summary>
		/// Constant function to estimate theoretic effectivness of HashedSpace
		/// </summary>
		/// <returns>true if effective</returns>
		constexpr bool is_effective() const noexcept {
			return sep_[0] > 3 || sep_[1] > 3 || sep_[2] > 3;
		}

		/// <summary>
		/// Creates theoretical overestimated vector of Bonds
		/// </summary>
		/// <param name="points"> - vector of Points</param>
		/// <returns>vector with all bonds in boxes and between adjacent ones</returns>
		template<BondConcept BondType>
		::std::vector<BondType> create_hash_bonds(const ::std::vector<PointType>& points) const {
			::std::vector<BondType> ret;
			size_t estimated_size = points.size() * points.size() * sizemod();
			ret.reserve(estimated_size);
			SupType supply_table(sep_[0],
								 SupType::value_type(sep_[1],
													 SupType::value_type::value_type(sep_[2])));

			// Fill supply_table
			for (AtomIndex i = 0; i < points.size(); i++)
			{
				auto c = coordintate_of_point(points[i]);
				supply_table[c[0]][c[1]][c[2]].emplace_back(i);
			}

			// Create all bonds
			for (size_t i = 0; i < sep_[0]; i++) {
				for (size_t j = 0; j < sep_[1]; j++) {
					for (size_t k = 0; k < sep_[2]; k++) {
						box_working(ret, supply_table, i, j, k);
					}
				}
			}
			return ret;
		}

	private:
		static constexpr FloatingPointType modifier_ = 1.05;

		template<BondConcept BondType>
		void box_working(::std::vector<BondType>& ret, const SupType& supply_table, size_t i, size_t j, size_t k) const {
			create_bonds_in_box(ret, supply_table[i][j][k]);

			bool is_x = sep_[0] != 1;
			bool is_y = sep_[1] != 1;
			bool is_z = sep_[2] != 1;

			auto dx = (i + 1 == sep_[0]) ? 0 : i + 1;
			auto dy = (j + 1 == sep_[1]) ? 0 : j + 1;
			auto dz = (k + 1 == sep_[2]) ? 0 : k + 1;

			// dx
			if (is_x) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][j][k]);
			}
			// dy
			if (is_y) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][dy][k]);
			}
			// dz
			if (is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][j][dz]);
			}

			// dxdy
			if (is_x && is_y) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][dy][k]);
			}
			// dxdz
			if (is_x && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][j][dz]);
			}
			// dydz
			if (is_y && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[i][dy][dz]);
			}

			// dxdydz
			if (is_x && is_y && is_z) {
				create_bonds_between_boxes(ret, supply_table[i][j][k], supply_table[dx][dy][dz]);
			}
		}
		template<BondConcept BondType>
		void create_bonds_in_box(::std::vector<BondType>& ret, const SupListType& l) const {
			for (auto iter1 = l.begin(); iter1 != l.end(); iter1++)
			{
				auto iter2 = iter1;
				iter2++;
				for (; iter2 != l.end(); iter2++)
				{
					ret.emplace_back(*iter1, *iter2);
				}
			}
		}
		template<BondConcept BondType>
		void create_bonds_between_boxes(::std::vector<BondType>& ret, const SupListType& l1, const SupListType& l2) const {
			for (auto v1 : l1)
			{
				for (auto v2 : l2)
				{
					if (v2 > v1) {
						ret.emplace_back(v1, v2);
					}
					else {
						ret.emplace_back(v2, v1);
					}
				}
			}
		}

		constexpr SupPoint coordintate_of_point(const PointType& p) const noexcept {
			return {
				static_cast<size_t>(std::floor(p[0] * sep_[0])),
				static_cast<size_t>(std::floor(p[1] * sep_[1])),
				static_cast<size_t>(std::floor(p[2] * sep_[2]))
			};
		}
		constexpr FloatingPointType sizemod() const noexcept {
			FloatingPointType ret = 1;
			for (DimentionType i = 0; i < 3; i++)
			{
				if (sep_[i] > 3) {
					ret *= 3;
					ret /= sep_[i];
					ret *= modifier_;
				}
			}
			return ret;
		}
		constexpr void calculateSep(FloatingPointType maxbond) {
			maxbond *= modifier_;
			for (DimentionType i = 0; i < 3; i++) {
				sep_[i] = static_cast<size_t>(floor(cell_.lat_dir(i) / maxbond));
				if (sep_[i] == 0) sep_[i] = 1;
			}
		}

	private:
		// Data
		const CellType& cell_;
		::std::array<size_t, 3> sep_;
	};

}

template<class T>
struct std::hash<cpplib::geometry::Point<T>>
{
	std::size_t operator()(const cpplib::geometry::Point<T>& s) const noexcept
	{
		std::size_t h1 = std::hash<T>{}(s[0]);
		std::size_t h2 = std::hash<T>{}(s[1]);
		std::size_t h3 = std::hash<T>{}(s[2]);
		return h1 ^ (h2 << 2) ^ (h3 << 4);
	}
};