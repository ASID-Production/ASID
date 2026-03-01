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
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

/**
			 * Convert an angle from degrees to radians.
			 *
			 * @tparam T Floating or integral numeric type of the input and result.
			 * @param a Angle in degrees.
			 * @return Angle in radians.
			 */
			/**
			 * Convert an angle from radians to degrees.
			 *
			 * @tparam T Floating or integral numeric type of the input and result.
			 * @param a Angle in radians.
			 * @return Angle in degrees.
			 */
			/**
			 * Tolerance for comparing positions in real (Cartesian) space, in Angstroms.
			 *
			 * Value: 0.01
			 */
			 /**
			 * Tolerance for comparing positions in fractional (unit-cell) space.
			 *
			 * Computed as crystallography_eq_position_eps_realspace divided by 100.
			 */
			 /**
			  * Three-dimensional point or vector with components of type T.
			  *
			  * Provides storage for three components, common vector operations, and a
			  * lightweight hash functor specialized for floating and integral component types.
			  *
			  * The static member `eq_pos` is a small epsilon appropriate for fractional-space
			  * position comparisons.
			  *
			  * @tparam T Component value type.
			  */
			namespace cpplib::geometry {
	template <class T> inline T GradtoRad(T a) {
		return a * static_cast<T>(0.0174532925199432957692);
	}
	template <class T> inline T RadtoGrad(T a) {
		return a * static_cast<T>(57.295779513082320877);
	}

	constexpr double crystallography_eq_position_eps_realspace = 0.01; // Angstrom
	constexpr double crystallography_eq_position_eps_fractalspace = crystallography_eq_position_eps_realspace / 100;

	template<class T>
	struct Point {
	public:
		static constexpr T eq_pos = static_cast<T>(crystallography_eq_position_eps_fractalspace);

		using value_type = T;
		using array_type = ::std::array<value_type, 3>;

		struct Hash {
			std::size_t operator()(const Point<T>& point) const {
				if constexpr (std::is_floating_point_v<T>) {
					std::size_t hx = std::hash<T>{}(point[0]);
					std::size_t hy = std::hash<T>{}(point[1]);
					std::size_t hz = std::hash<T>{}(point[2]);

					return hx ^ (hy << 1) ^ (hz << 2) ^ (hx >> 31);
				} else {
					return
						(point[0] * std::size_t(73856093)) ^
						(point[1] * std::size_t(19349663)) ^
						(point[2] * std::size_t(83492791));
				}
			}
		};


	public:
		array_type a = {0,0,0};

	public:
		// Constructors
		constexpr Point() noexcept = default;
		constexpr Point(value_type x, value_type y, value_type z) noexcept : a{x, y, z} {};
		/**
 * Construct a Point from a 3-element array.
 * @param other Array containing the components (x, y, z) in that order; values are copied into the Point.
 */
explicit constexpr Point(const array_type& other) noexcept : a(other) {};
		/**
 * Construct a Point by taking ownership of a 3-element array of components.
 * @param other Array of three coordinates whose contents are moved into the new Point.
 */
explicit constexpr Point(array_type&& other) noexcept : a(::std::move(other)) {};

		template <typename T2>
			requires ((::std::integral<T2> || ::std::floating_point<T2>) && ::std::is_convertible<T2, T>::value)
		explicit constexpr Point(const Point<T2>& other) noexcept {
			a[0] = static_cast<T>(other[0]);
			a[1] = static_cast<T>(other[1]);
			a[2] = static_cast<T>(other[2]);
		}
		/// @brief Calculate distance to [0,0,0]. Don't use as (point[a]-point[b]).r().
		/**
		 * Compute the Euclidean norm of the point (distance from the origin).
		 * @returns Euclidean distance (sqrt of sum of squared components).
		 */
		 
		/**
		 * Compute the squared Euclidean distance of the point (distance^2 from the origin).
		 * Useful for comparisons when avoiding a square root.
		 * @returns Sum of squared components (distance squared).
		 */
		constexpr value_type r() const noexcept {
			return sqrt(fma(a[0], a[0], fma(a[1], a[1], a[2] * a[2])));
		}
		/// @brief Calculate square of distance to [0,0,0]. Use for comparisons
		/// @return Distance^2
		constexpr value_type rSq() const noexcept {
			return fma(a[0], a[0], fma(a[1], a[1], a[2] * a[2]));
		}

		/**
		 * Wrap the point's components into the unit cell [0,1) by subtracting each component's floor.
		 *
		 * @return Reference to this Point after its components have been wrapped into [0,1).
		 */
		constexpr Point& MoveToCell() noexcept {
			a[0] -= ::std::floor(a[0]);
			a[1] -= ::std::floor(a[1]);
			a[2] -= ::std::floor(a[2]);
			return *this;
		}

		/**
		 * Compute the dot (scalar) product of two 3D points/vectors.
		 *
		 * @param left Left operand vector.
		 * @param right Right operand vector.
		 * @returns The scalar dot product left.x*right.x + left.y*right.y + left.z*right.z.
		 */
		[[nodiscard]] static constexpr value_type Scalar(const Point& left, const Point& right) noexcept {
			return (left.a[0] * right.a[0] + left.a[1] * right.a[1] + left.a[2] * right.a[2]);
		}
		/**
		 * Compute the cross product of two 3D vectors.
		 * @param left Left-hand vector.
		 * @param right Right-hand vector.
		 * @returns A Point containing the cross product (left × right).
		 */
		[[nodiscard]] static constexpr Point Vector(const Point& left, const Point& right) noexcept {
			return Point(left.a[1] * right.a[2] - left.a[2] * right.a[1], left.a[2] * right.a[0] - left.a[0] * right.a[2], left.a[0] * right.a[1] - left.a[1] * right.a[0]);
		}
		/**
		 * Compute Euclidean distance between two 3D points.
		 *
		 * @param a First point.
		 * @param b Second point.
		 * @returns Euclidean distance between `a` and `b`.
		 */
		static constexpr value_type distance(const Point& a, const Point& b) noexcept {
			value_type d0 = a.a[0] - b.a[0];
			value_type d1 = a.a[1] - b.a[1];
			value_type d2 = a.a[2] - b.a[2];
			return sqrt(fma(d0, d0, fma(d1, d1, d2 * d2)));
		}
		/**
		 * Compute the squared Euclidean distance between two points.
		 * @param a First point.
		 * @param b Second point.
		 * @returns Squared distance (sum of squared component differences) between `a` and `b`.
		 */
		static constexpr value_type distanceSq(const Point& a, const Point& b) noexcept {
			value_type d0 = a.a[0] - b.a[0];
			value_type d1 = a.a[1] - b.a[1];
			value_type d2 = a.a[2] - b.a[2];
			return fma(d0, d0, fma(d1, d1, d2 * d2));
		}
		/**
		 * Compute the Euclidean distance between two points using the minimum-image convention
		 * in a unit cubic periodic cell.
		 *
		 * @param a First point (fractional coordinates in the cubic cell).
		 * @param b Second point (fractional coordinates in the cubic cell).
		 * @returns `value_type` distance between `a` and `b` accounting for periodic boundaries (minimum-image).
		 */
		static constexpr value_type distanceInCubicCell(const Point& a, const Point& b) noexcept {
			value_type d0 = fmod(a.a[0] - b.a[0] + T(0.5), T(1.0)) - T(0.5);
			value_type d1 = fmod(a.a[1] - b.a[1] + T(0.5), T(1.0)) - T(0.5);
			value_type d2 = fmod(a.a[2] - b.a[2] + T(0.5), T(1.0)) - T(0.5);
			return sqrt(fma(d0, d0, fma(d1, d1, d2 * d2)));
		}

		/**
		 * Determine whether two points are equal within a tolerance when wrapped into a cubic periodic cell.
		 *
		 * @param a First point in fractional coordinates.
		 * @param b Second point in fractional coordinates.
		 * @param epsilon Maximum allowed absolute difference per coordinate after applying periodic wrapping.
		 * @returns `true` if each coordinate difference, wrapped into the interval [-0.5, 0.5), is less than or equal to `epsilon`, `false` otherwise.
		 */
		static constexpr value_type isSameInCubicCell(const Point& a, const Point& b, T epsilon) noexcept {
			return abs(fmod(a.a[0] - b.a[0] + T(0.5), T(1.0)) - T(0.5)) <= epsilon &&
				abs(fmod(a.a[1] - b.a[1] + T(0.5), T(1.0)) - T(0.5)) <= epsilon &&
				abs(fmod(a.a[2] - b.a[2] + T(0.5), T(1.0)) - T(0.5)) <= epsilon;
		}
		/**
		 * Check whether two points are equal within a per-coordinate tolerance.
		 *
		 * @param epsilon Maximum allowed absolute difference for each coordinate (same units as point components).
		 * @returns `true` if for every coordinate the absolute difference is less than or equal to `epsilon`, `false` otherwise.
		 */
		static constexpr bool isSame(const Point& a, const Point& b, T epsilon) noexcept {
			return abs(a.a[0] - b.a[0]) <= epsilon &&
				abs(a.a[1] - b.a[1]) <= epsilon &&
				abs(a.a[2] - b.a[2]) <= epsilon;
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
		/**
		 * Round each component of the point to the nearest integer.
		 * @returns A Point whose components are rounded to the nearest integer.
		 */
		/**
		 * Floor each component of the point to the greatest integer less than or equal to it.
		 * @returns A Point whose components are floored.
		 */
		/**
		 * Quantize a point's components to the nearest multiple of epsilon.
		 * @param p Point to quantize.
		 * @param epsilon Grid spacing used for quantization; if `epsilon` is zero an all-zero Point is returned.
		 * @returns A Point whose components are rounded to the nearest multiple of `epsilon`.
		 */
		/**
		 * Access a component by index.
		 * @param i Component index (0, 1, or 2).
		 * @returns The value of the requested component.
		 */
		/**
		 * Access a component by index for modification.
		 * @param i Component index (0, 1, or 2).
		 * @returns A reference to the requested component.
		 */
		/**
		 * Negate all components of the point.
		 * @returns A Point with each component negated.
		 */
		/**
		 * Component-wise addition of two points.
		 * @param left Left-hand point operand.
		 * @param right Right-hand point operand.
		 * @returns A Point whose components are the component-wise sums; resulting value type is deduced from operand types.
		 */
		/**
		 * Add a scalar to each component of a point.
		 * @param left Point operand.
		 * @param b Scalar to add to each component.
		 * @returns A Point with each component increased by `b`; resulting value type is deduced from the operands.
		 */
		constexpr Point round() const {
			return Point(std::round(a[0]), std::round(a[1]), std::round(a[2]));
		}
		constexpr Point floor() const {
			return Point(std::floor(a[0]), std::floor(a[1]), std::floor(a[2]));
		}
		static constexpr Point quantize(const Point& p, T epsilon) noexcept {
			if (epsilon == 0) return {};
			return {
				std::round(p[0] / epsilon) * epsilon,
				std::round(p[1] / epsilon) * epsilon,
				std::round(p[2] / epsilon) * epsilon
			};
		}

		constexpr value_type operator[](const uint8_t i) const noexcept {
			return a[i];
		}
		constexpr value_type& operator[](const uint8_t i) noexcept {
			return a[i];
		}

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
			a[0] -= static_cast<T>(right.a[0]);
			a[1] -= static_cast<T>(right.a[1]);
			a[2] -= static_cast<T>(right.a[2]);
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
		array_type A{{ {0,0,0},{0,0,0},{0,0,0} }};
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
		/**
 * Constructs a 3x3 diagonal matrix with the given scalar on the main diagonal and zeros elsewhere.
 * @param v Scalar value placed on each diagonal element (A[0][0], A[1][1], A[2][2]).
 */
explicit constexpr Matrix(const T v) noexcept : A{{ {v,0,0},{0,v,0},{0,0,v} }} {}
		/**
		 * Construct a 3x3 matrix by copying values from a C-style 2D array.
		 * @param input_massive Pointer to an array of at least three pointers, each pointing to at least three elements; elements are read as input_massive[row][col] and copied into the matrix in row-major order.
		 */
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
		/**
		 * Compute the product of this matrix and another 3×3 matrix.
		 * @tparam T2 Type of elements in the right-hand matrix.
		 * @param right Right-hand operand matrix.
		 * @returns Matrix representing the standard matrix product; each element's type is the result of multiplying `T` by `T2`.
		 */
		[[nodiscard]] constexpr Matrix<decltype(T()* T2())> operator*(const Matrix<T2>& right) const noexcept {
			using resv = decltype(T()* T2());
			std::array<resv, 3> a1 = {A[0][0] * right.A[0][0] + A[0][1] * right.A[1][0] + A[0][2] * right.A[2][0],
										A[0][0] * right.A[0][1] + A[0][1] * right.A[1][1] + A[0][2] * right.A[2][1],
										A[0][0] * right.A[0][2] + A[0][1] * right.A[1][2] + A[0][2] * right.A[2][2]};
			std::array<resv, 3> a2 = {A[1][0] * right.A[0][0] + A[1][1] * right.A[1][0] + A[1][2] * right.A[2][0],
										A[1][0] * right.A[0][1] + A[1][1] * right.A[1][1] + A[1][2] * right.A[2][1],
										A[1][0] * right.A[0][2] + A[1][1] * right.A[1][2] + A[1][2] * right.A[2][2]};
			std::array<resv, 3> a3 = {A[2][0] * right.A[0][0] + A[2][1] * right.A[1][0] + A[2][2] * right.A[2][0],
										A[2][0] * right.A[0][1] + A[2][1] * right.A[1][1] + A[2][2] * right.A[2][1],
										A[2][0] * right.A[0][2] + A[2][1] * right.A[1][2] + A[2][2] * right.A[2][2]};
			std::array<std::array<resv, 3>, 3> b{a1, a2, a3};
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
			std::array<T, 3> a1 = {A[0][0] / right, A[0][1] / right, A[0][2] / right};
			std::array<T, 3> a2 = {A[1][0] / right, A[1][1] / right, A[1][2] / right};
			std::array<T, 3> a3 = {A[2][0] / right, A[2][1] / right, A[2][2] / right};
			array_type b{a1, a2, a3};
			return Matrix(std::move(b));
		}
		[[nodiscard]] constexpr Matrix<T> Transponate() const noexcept {
			std::array<T, 3> a1 = {A[0][0],A[1][0],A[2][0]};
			std::array<T, 3> a2 = {A[0][1],A[1][1],A[2][1]};
			std::array<T, 3> a3 = {A[0][2],A[1][2],A[2][2]};
			array_type b{a1, a2, a3};
			return Matrix(b);
		}
		[[nodiscard]] constexpr Matrix<T> Invert() const {
			const T det = Det();
			std::array<T, 3> a1 = {(A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det, (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det, (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det};
			std::array<T, 3> a2 = {(A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det, (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det, (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / det};
			std::array<T, 3> a3 = {(A[1][0] * A[2][1] - A[2][0] * A[1][1]) / det, (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / det, (A[1][1] * A[0][0] - A[1][0] * A[0][1]) / det};
			array_type b = {a1,a2,a3};
			return Matrix(b);
		}
		[[nodiscard]] constexpr Matrix<T> Modul() const noexcept {
			constexpr T zero = 0;
			std::array<T, 3> a1 = {A[0][0] < zero?A[0][0]:-A[0][0], A[0][1] < zero?A[0][1]:-A[0][1], A[0][2] < zero?A[0][2]:-A[0][2]};
			std::array<T, 3> a2 = {A[1][0] < zero?A[1][0]:-A[1][0], A[1][1] < zero?A[1][1]:-A[1][1], A[1][2] < zero?A[1][2]:-A[1][2]};
			std::array<T, 3> a3 = {A[2][0] < zero?A[2][0]:-A[2][0], A[2][1] < zero?A[2][1]:-A[2][1], A[2][2] < zero?A[2][2]:-A[2][2]};
			array_type b{a1, a2, a3};
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

		template<class T2>
		constexpr Point<decltype(T()* T2())> TransposeMultiply(const Point<T2>& right) const noexcept {
			Point<decltype(T()* T2())> res;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					res[i] = std::fma(A[j][i], right[j], res[i]);
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
		std::array<T, 4> a = {T(0), T(0), T(0), T(0)}; // [ A, B, C, D]

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
			assert(check_corectness());
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
			} else {
				angleRad_[0] = alpha;
				angleRad_[1] = beta;
				angleRad_[2] = gamma;
				for (int i = 0; i < 3; i++) {
					angleGrad_[i] = RadtoGrad(angleRad_[i]);
				}
			}
			createMatrix();
			assert(check_corectness());
		}

		/**
		 * Access the cell angle in radians by index.
		 *
		 * @param i Index of the angle: 0 = alpha, 1 = beta, 2 = gamma.
		 * @returns const reference to the requested angle value (radians).
		 */
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
		/**
		 * Access the 3x3 transformation matrix that maps fractional (direct lattice) coordinates to Cartesian coordinates.
		 * @returns Reference to the matrix used to convert fractional coordinates into Cartesian coordinates.
		 */
		 
		/**
		 * Access the 3x3 transformation matrix that maps Cartesian coordinates to fractional (direct lattice) coordinates.
		 * @returns Reference to the matrix used to convert Cartesian coordinates into fractional coordinates.
		 */
		[[nodiscard]] constexpr const matrix_type& fracToCart() const noexcept {
			return fracToCart_;
		}
		[[nodiscard]] constexpr const matrix_type& cartToFrac() const noexcept {
			return cartToFrac_;
		}
		/**
			 * Check that Cartesian lengths of selected neighbor fractional shifts meet minimum lattice constraints.
			 *
			 * Tests ten predefined neighboring fractional shifts by transforming them to Cartesian coordinates and
			 * comparing their squared lengths to the square of the smallest lattice edge. A tolerance of 1e-5 is
			 * applied to the comparison.
			 *
			 * @returns `true` if all tested neighbor displacements have squared Cartesian length greater than or
			 *          equal to the smallest lattice edge squared minus 1e-5, `false` otherwise.
			 */
			bool check_corectness() const noexcept {
			constexpr std::array<Point<int8_t>, 10> shiftTable = {{
				{-1, -1, -1}, { 0, -1, -1}, { 1, -1, -1},
				{-1,  0, -1}, { 1,  0, -1}, {-1,  1, -1}, 
				{ 0,  1, -1}, { 1,  1, -1}, {-1, -1,  0}, 
				{ 1, -1,  0}
			}};

			T minlatSq = std::min(lattice_[0], std::min(lattice_[1], lattice_[2]));
			minlatSq *= minlatSq;

			for (char i = 0; i < 10; i++) {
				auto rsq = (fracToCart_ * shiftTable[i]).rSq();
                if (minlatSq - rsq > 1e-5)
					return false;
			}
			return true;
		}

		/**
		 * Compute the shortest Cartesian distance between two points using periodic boundary conditions on fractional coordinates.
		 *
		 * @param a First point given in fractional (direct) coordinates; expected to be within [0,1) in each component.
		 * @param b Second point given in fractional (direct) coordinates; expected to be within [0,1) in each component.
		 * @returns The Euclidean distance in Cartesian space after applying the minimum-image convention in the unit cell.
		 */
		value_type distance_in_01(const PointType& a, const PointType& b) const {
			PointType d = (a - b).MoveToCell();
			if (d[0] > 0.5) d[0] = 1 - d[0];
			if (d[1] > 0.5) d[1] = 1 - d[1];
			if (d[2] > 0.5) d[2] = 1 - d[2];
			return (fracToCart_ * d).r();
		}

		template <class I>
		/**
		 * Compute minimal integer supercell dimensions that satisfy a distance cutoff for neighbor directions.
		 *
		 * Ensures each of a predefined set of neighbor direction vectors, when scaled by the returned
		 * supercell and converted to Cartesian coordinates, has length greater than `cutoff`. Each
		 * returned component is at least `minimum`.
		 *
		 * @param cutoff Distance cutoff (same units as the cell lattice lengths).
		 * @param minimum Minimum allowed repeat count for each supercell dimension.
		 * @returns A Point<I> whose three components are the integer repeat counts along the cell axes
		 *          that meet the cutoff requirement and are each >= `minimum`.
		 */
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
				geometry::Point<I>(1,-1,-1)};

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
	struct Symm {
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
		bool isValid() const noexcept {
			return abs(point[0]) < 1 &&
				abs(point[1]) < 1 &&
				abs(point[2]) < 1 &&
				abs(abs(mat.Det()) - 1.0) < 0.0001;
		}
	private:
		size_t findcomma(const char* str) const {
			size_t n = 0;
			for (; str[n] != ',' && str[n] != '\0'; n++);
			return n;
		}
		std::pair<point_t, T> parse(const char* str, const size_t len) const {
			/**
		 * Parse a symmetry domain substring into a 3D signed axis mask and a numeric shift.
		 *
		 * Parses characters in `str[0..len-1]` interpreting letters 'x','y','z' (case-insensitive)
		 * as axis indicators that set the corresponding component of the returned point to -1 or 1
		 * depending on preceding '-' signs, and numeric tokens (integers or fractions) as additive
		 * components of the returned shift value. Spaces and quote characters are ignored; '+' is ignored;
		 * each '-' toggles the sign for the next axis or numeric token. The iteration index `i`
		 * is advanced by the numeric parser when a numeric token is consumed.
		 *
		 * @param str Pointer to the character sequence containing the domain substring to parse.
		 * @param len Number of characters available in `str`.
		 * @returns std::pair<point_t, T> where the first element is a point_t with components in {-1,0,1}
		 *          indicating axis directions found, and the second element is the cumulative shift value.
		 */
		point_t p{0,0,0};
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
		/**
		 * Parse a numeric shift value from a character buffer supporting integer, decimal, and simple fractional formats.
		 *
		 * Reads characters from `str` starting at `iter` up to `len` and parses a value of type `T` in one of these forms:
		 * - integer: "12"
		 * - decimal: "1.25"
		 * - fraction: "3/4" (numerator/denominator; denominator may have multiple digits)
		 *
		 * The function advances `iter` while consuming digits and format separators; before returning it decrements `iter`
		 * so that `iter` points to the last consumed character. If a '+' or '-' sign is encountered, the function steps
		 * the iterator back and returns the value parsed so far. Any unexpected non-numeric character terminates parsing.
		 *
		 * @tparam T Numeric type to return (e.g., float, double).
		 * @param str Null-terminated character buffer containing the textual number.
		 * @param iter Index in `str` where parsing starts; updated to the position of the last consumed character (decremented once before return).
		 * @param len Length of the buffer `str` (upper bound for `iter`).
		 * @returns Parsed numeric value converted to type `T`.
		 */
		T parseshift(const char* str, unsigned int& iter, const size_t len) const {
			bool dot = false;
			bool slash = false;

			int upper = 0;
			int lower = 1;
			for (; iter < len; iter++)
			{
				if (str[iter] < '0' || str[iter] > '9') {
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
				} else {
					int num = (str[iter] - '0');
					if (dot) {
						upper = upper * 10 + num;
						lower *= 10;
					} else if (slash) {
						lower = lower * 10 + num;
					} else {
						upper = upper * 10 + num;
					}
				}

			}
			iter--;
			return upper / static_cast<T>(lower);
		}
	};

	/**
		 * Compactly encodes a 3D lattice shift in the range [-1,1]^3 as a single byte code (0..26).
		 *
		 * Provides construction from an existing code, from a 3-component shift vector, or default-initialization.
		 *
		 * @note The encoded space includes all 27 shifts from (-1,-1,-1) to (1,1,1).
		 *
		 * @param c Code value in the range 0..26; assertion fires if the value is outside this range.
		 * @param sp 3-component shift point with each component in the range -1..1; assertion fires if any component is outside that range.
		 */
		struct ShiftCode {
		using ShiftPoint = Point<int8_t>;

		constexpr ShiftCode() noexcept = default;
		constexpr explicit ShiftCode(uint8_t c) noexcept : code(c) {
			assert(c <= 26);
		}
		constexpr explicit ShiftCode(ShiftPoint sp) noexcept {
			assert(sp[0] >= -1 && sp[0] <= 1);
			assert(sp[1] >= -1 && sp[1] <= 1);
			assert(sp[2] >= -1 && sp[2] <= 1);
			code = compress_shift(sp);
		}

		ShiftCode& operator=(uint8_t c) noexcept {
			assert(c <= 26);
			code = c;
			return *this;
		}

		ShiftCode& operator=(ShiftPoint sp) noexcept {
			assert(sp[0] >= -1 && sp[0] <= 1);
			assert(sp[1] >= -1 && sp[1] <= 1);
			assert(sp[2] >= -1 && sp[2] <= 1);
			code = compress_shift(sp);
			return *this;
		}

		constexpr ShiftPoint get_shift() const noexcept {
			return shiftTable[code];
		}
		constexpr uint8_t get_code() const noexcept {
			return code;
		}
		/**
		 * Map a compact shift code to its 3D shift vector.
		 *
		 * @param c Shift code in the range 0..26 that encodes a 3D shift with components in {-1, 0, 1}.
		 *             The mapping follows the internal shift table ordering from (-1,-1,-1) to (1,1,1).
		 *             The function asserts if `c` is greater than 26.
		 * @returns ShiftPoint corresponding to code `c`; each component is in {-1, 0, 1}.
		 */
		static constexpr ShiftPoint get_shift(uint8_t c) noexcept {
			assert( c <= 26);
			return shiftTable[c];
		}
		/**
		 * Replace the stored 3D shift code with its inverse (the opposite shift).
		 *
		 * The operation is performed in-place on the object's internal code.
		 */
		constexpr void inverse() noexcept {
			code = 26 - code;
		}
		/**
		 * Compute the inverse shift code in the 27-entry shift table.
		 *
		 * @param c Shift code in the range 0..26.
		 * @returns `uint8_t` equal to 26 - `c`, the inverse code.
		 */
		static constexpr uint8_t inverse(uint8_t c) noexcept {
			assert(c <= 26);
			return 26 - c;
		}
		/**
		 * Compute the inverse (negated) shift code for a given encoded shift.
		 *
		 * @param sc Encoded shift to invert; valid codes are 0..26 representing shifts from (-1,-1,-1) to (1,1,1).
		 * @returns A ShiftCode whose internal code is 26 - sc.code, representing the opposite shift.
		 */
		static constexpr ShiftCode inverse(ShiftCode sc) noexcept {
			assert(sc.code <= 26);
			return ShiftCode(26 - sc.code);
		}
				
		static constexpr std::array<ShiftPoint, 27> shiftTable = {{
			{-1, -1, -1}, { 0, -1, -1}, { 1, -1, -1}, // code 0, 1, 2
			{-1,  0, -1}, { 0,  0, -1}, { 1,  0, -1}, // code 3, 4, 5
			{-1,  1, -1}, { 0,  1, -1}, { 1,  1, -1}, // code 6, 7, 8

			{-1, -1,  0}, { 0, -1,  0}, { 1, -1,  0}, // code 9, 10, 11
			{-1,  0,  0}, { 0,  0,  0}, { 1,  0,  0}, // code 12, 13 (Center), 14
			{-1,  1,  0}, { 0,  1,  0}, { 1,  1,  0}, // code 15, 16, 17

			{-1, -1,  1}, { 0, -1,  1}, { 1, -1,  1}, // code 18, 19, 20
			{-1,  0,  1}, { 0,  0,  1}, { 1,  0,  1}, // code 21, 22, 23
			{-1,  1,  1}, { 0,  1,  1}, { 1,  1,  1}  // code 24, 25, 26
		}};
	private:
		uint8_t code = 13;

		/// @brief Compress ShiftType[-1,+1] (usually Point<int8_t>) to char
		/// @param s The shift
		/**
		 * Encode a 3D shift with components in {-1, 0, 1} into a compact 0–26 code.
		 * @param s ShiftPoint with each component equal to -1, 0, or 1.
		 * @returns Compressed code in range 0..26 computed as (s[0]+1) + (s[1]+1)*3 + (s[2]+1)*9.
		 */
		static constexpr uint8_t compress_shift(ShiftPoint s) {
			return (s[0] + uint8_t(1)) +
				(s[1] + uint8_t(1)) * uint8_t(3) +
				(s[2] + uint8_t(1)) * uint8_t(9);
		}
	};


	/// @brief Class for using spartial hashing algorithm to find all bonds in 3D periodic space.
	/// @tparam T Floating point type
	template <class T>
	struct SpatialGrid {
		using ShiftType = Point<int8_t>;
		using PointType = Point<T>;
		using CellType = Cell<T>;
		struct BondWithShift {
			int first = 0;
			int second = 0;
			ShiftCode shiftcode = 13;
			BondWithShift() = default;
			BondWithShift(int a, int b)
				: first(a), second(b) {}
			BondWithShift(int a, int b, ShiftCode scode)
				: first(a), second(b), shiftcode(scode) {}

		};

		struct VirtualNeighbour {
			int realBoxIndex;
			ShiftCode shiftcode;
		};

		std::vector<int> pointIndices;          // [N] All point indices in box order
		std::vector<int> boxOffsets;            // [C+1] Offsets in pointIndices array where each virtual box starts
		std::vector<int> realBoxOffsets;        // [C_real+1] Offsets for real boxes
		std::vector<VirtualNeighbour> virtMap;  // Mapping from virtual box to real box

		std::array<uint8_t, 3> gridDim = {1,1,1};      // Real grid dimensions
		std::array<uint8_t, 3> gridDimVirt = {1,1,1};  // Virtual grid dimensions
		int numBoxes = 1;       // number of real boxes
		int numBoxesVirt = 1;   // number of virtual boxes

		std::array<int, 13> left_boxes_shifts;

		/// @brief Build the spatial grid from a set of points. 
		/// 
		/// @note All points must have coordinates in[0, 1). Unnormalized coordinates
		/// will produce undefined behavior(out - of - bounds access).
		/// 
		/// @param points Vector of points with coordinates normalized to[0, 1) in fractional space
		/// @param cell   The unit cell definition
		/// @param cutoff Distance cutoff for bonding
		void build(const std::vector<Point<T>>& points, const CellType& cell, const T cutoff) {
			// Calculate grid dimensions
			calculateGridDim(cell, cutoff);

			// Prepare vectors
			std::vector<int> realBoxCount(numBoxes, 0);
			std::vector<int> virtBoxCount(numBoxesVirt, 0);
			realBoxOffsets.resize(numBoxes + 1, 0);
			boxOffsets.resize(numBoxesVirt + 1, 0);

			// Initialize virtual box mapping
			virtMap.resize(numBoxesVirt);

			// Build virtual-to-real mapping with periodic boundary conditions
			build_virtual_mapping();

			// Temporary array for storing point-to-virtual-box assignment
			std::vector<int> temp_virt_box_IDx(points.size(), 0);
			pointIndices.resize(points.size(), 0);

			// Count atoms in virtual boxes
			for (size_t i = 0; i < points.size(); i++) {
				auto temp = get_virtual_box_index(points[i]);
				temp_virt_box_IDx[i] = temp;
				virtBoxCount[temp]++;

				// Also count for real boxes (for fast access)
				realBoxCount[get_real_box_index(points[i])]++;
			}

			// Build offsets for virtual boxes
			int currentOffset = 0;
			for (int i = 0; i < numBoxesVirt; i++) {
				boxOffsets[i] = currentOffset;
				currentOffset += virtBoxCount[i];
				virtBoxCount[i] = 0;  // Reset for filling
			}
			boxOffsets[numBoxesVirt] = currentOffset;

			// Build offsets for real boxes
			currentOffset = 0;
			for (int i = 0; i < numBoxes; i++) {
				realBoxOffsets[i] = currentOffset;
				currentOffset += realBoxCount[i];
			}
			realBoxOffsets[numBoxes] = currentOffset;

			// Fill pointIndices in virtual box order
			for (size_t i = 0; i < points.size(); i++) {
				int vIdx = temp_virt_box_IDx[i];
				int destPos = boxOffsets[vIdx] + virtBoxCount[vIdx];
				pointIndices[destPos] = i;
				virtBoxCount[vIdx]++;
			}

			// Calculate shifts for 13 left boxes
			auto baseshift = get_box_by_index(1, 1, 1, gridDimVirt);
			for (int i = 0; i < 13; i++) {
				auto temp = get_box_by_index(ShiftCode::shiftTable[i][0] + 1,
											 ShiftCode::shiftTable[i][1] + 1,
											 ShiftCode::shiftTable[i][2] + 1,
											 gridDimVirt);
				left_boxes_shifts[i] = temp - baseshift;
			}
		}
		/// @brief Create bonds, based on the spatial grid.
		/// @param double_sided Default: false. Boolian, whether to create bonds in both directions.
		/// @return vector of all possible bonds
		std::vector<BondWithShift> get_bonds(bool double_sided = false) {
			std::vector<BondWithShift> bonds;
			// Preliminary memory reservation to reduce reallocations
			bonds.reserve(pointIndices.size() * (double_sided?26:13));

			// Iterate through real grid dimensions
			for (int rz = 0; rz < gridDim[2]; ++rz) {
				for (int ry = 0; ry < gridDim[1]; ++ry) {
					for (int rx = 0; rx < gridDim[0]; ++rx) {
						process_box_bonds(rx, ry, rz, bonds, double_sided);
					}
				}
			}
			return bonds;
		}



	private:
		// Build virtual box mapping
		void build_virtual_mapping() {
			const int virtDimX = gridDimVirt[0];
			const int virtDimY = gridDimVirt[1];
			const int virtDimXY = virtDimX * virtDimY;

			for (int virtIdx = 0; virtIdx < numBoxesVirt; virtIdx++) {
				// Decompose linear index
				std::array<int, 3> v;
				v[2] = virtIdx / virtDimXY;
				v[1] = (virtIdx % virtDimXY) / virtDimX;
				v[0] = virtIdx % virtDimX;

				// Map virtual coordinates to real coordinates
				std::array<int, 3> r;
				ShiftType s(0, 0, 0);

				for (char j = 0; j < 3; j++)
				{
					const int dim = gridDim[j];
					r[j] = v[j] - 1;
					if (r[j] < 0) {
						s[j] = -1;
						r[j] += dim;
					} else if (r[j] >= dim) {
						r[j] -= dim;
						s[j] = 1;
					}

				}

				// Calculate real box index				
				int realIdx = get_box_by_index(r[0], r[1], r[2], gridDim);

				virtMap[virtIdx] = {
					.realBoxIndex = realIdx,
					.shiftcode = geometry::ShiftCode(s)
				};
			}
		}

		/**
		 * Process and append all bonds originating from the real grid box at (rx,ry,rz).
		 *
		 * Scans points contained in the specified real box and its 13 left-neighbour boxes (via the virtual mapping)
		 * and appends corresponding BondWithShift entries to `bonds`. Bonds between points inside the same real box
		 * are produced with no shift; bonds to neighbouring boxes carry the neighbour's stored shift code.
		 *
		 * @param rx Real-grid x index of the source box.
		 * @param ry Real-grid y index of the source box.
		 * @param rz Real-grid z index of the source box.
		 * @param bonds Output vector that will be appended with generated BondWithShift entries.
		 * @param double_sided If true, add both directions for each bond (also adds the inverse bond with inverted shift).
		 */
		void process_box_bonds(int rx, int ry, int rz, std::vector<BondWithShift>& bonds, bool double_sided) {
			// Current box in virtual grid (center of the 3x3x3 neighborhood)
			int vIdx = get_box_by_index(rx + 1, ry + 1, rz + 1, gridDimVirt);
			int rIdx = get_box_by_index(rx, ry, rz, gridDim);

			int start_a = realBoxOffsets[rIdx];
			int end_a = realBoxOffsets[rIdx + 1];

			// 1. Internal bonds: Shift code is always 13 (0,0,0)
			for (int i = start_a; i < end_a; ++i) {
				for (int j = i + 1; j < end_a; ++j) {
					add_bond_pair(pointIndices[i], pointIndices[j], ShiftCode(), bonds, double_sided);
				}
			}

			// 2. External bonds: Get shift code from the neighbor's virtual mapping
			for (int s = 0; s < 13; ++s) {
				int neighborVIdx = vIdx + left_boxes_shifts[s];
				int neighborRIdx = virtMap[neighborVIdx].realBoxIndex;

				// The shiftcode is stored in virtMap for each virtual cell
				ShiftCode sCode = virtMap[neighborVIdx].shiftcode;


				int start_b = realBoxOffsets[neighborRIdx];
				int end_b = realBoxOffsets[neighborRIdx + 1];

				for (int i = start_a; i < end_a; ++i) {
					for (int j = start_b; j < end_b; ++j) {
						add_bond_pair(pointIndices[i], pointIndices[j], sCode, bonds, double_sided);
					}
				}
			}
		}

		/**
		 * Append a directed bond from idxA to idxB to bonds and, if requested, also append the opposite bond.
		 * @param idxA Index of the source point.
		 * @param idxB Index of the destination point.
		 * @param shiftCode ShiftCode describing the periodic-image displacement from source to destination.
		 * @param bonds Vector to which the created BondWithShift(s) are appended.
		 * @param double_sided If true, also append the inverse bond (idxB -> idxA) with the shift code inverted.
		 */
		inline void add_bond_pair(int idxA, int idxB, ShiftCode shiftCode, std::vector<BondWithShift>& bonds, bool double_sided) const {
			// Basic bond a -> b
			bonds.push_back(BondWithShift{idxA, idxB, shiftCode});

			if (double_sided) {
				// Inverse bond b -> a
				// The shift for the opposite direction must be inverted
				ShiftCode sc = shiftCode;
				sc.inverse();
				bonds.push_back(BondWithShift{idxB, idxA, sc});
			}
		}

		/**
		 * Compute the linear index of the virtual grid cell that contains a point given in fractional coordinates.
		 *
		 * @param p Point with normalized fractional coordinates in [0,1); components correspond to positions along the three grid axes.
		 * @returns Linear index into the virtual-box indexing scheme that identifies the virtual cell containing `p`.
		 *
		 * @note Assumes `p` is already normalized to the half-open interval [0,1). Behavior is undefined if components fall outside this range.
		 */
		inline int get_virtual_box_index(const PointType& p) const {
			auto ix = static_cast<int>(p[0] * gridDim[0]) + 1;
			auto iy = static_cast<int>(p[1] * gridDim[1]) + 1;
			auto iz = static_cast<int>(p[2] * gridDim[2]) + 1;

			return get_box_by_index(ix, iy, iz, gridDimVirt);
		}
		/**
		 * Map a point with fractional coordinates in [0,1) to the linear index of its containing real grid box.
		 *
		 * The point `p` is interpreted in fractional (direct) coordinates; each component must be in the range [0,1).
		 * No validation is performed. The returned value is the linear real-box index computed from the point's
		 * grid coordinates and the current `gridDim`.
		 *
		 * @param p Point with normalized fractional coordinates in [0,1).
		 * @returns The linear index of the real grid box that contains `p`.
		 */
		inline int get_real_box_index(const PointType& p) const {
			auto ix = static_cast<int>(p[0] * gridDim[0]);
			auto iy = static_cast<int>(p[1] * gridDim[1]);
			auto iz = static_cast<int>(p[2] * gridDim[2]);
			return get_box_by_index(ix, iy, iz, gridDim);
		}

		// Grid dimensions are constrained by physical unit cell sizes (typically < 1000 Å).
		/**
		 * Compute grid dimensions and related counts for spatial partitioning based on the cell and cutoff.
		 *
		 * Sets `gridDim[i]` to floor(cell.lat_dir(i) / cutoff) with a minimum of 1, sets `gridDimVirt[i] = gridDim[i] + 2`
		 * (two extra layers for virtual boxes), and updates `numBoxes` and `numBoxesVirt` as the product of the
		 * respective dimensions.
		 *
		 * @param cell Cell describing lattice directions used to determine box sizes.
		 * @param cutoff Maximum distance determinating desired box size; must be greater than 0.
		 */
		constexpr void calculateGridDim(const CellType& cell, T cutoff) {
			assert(cutoff > static_cast<T>(0));
			for (uint8_t i = 0; i < 3; i++) {
				gridDim[i] = static_cast<uint8_t>(std::floor(cell.lat_dir(i) / cutoff));
				if (gridDim[i] == 0) gridDim[i] = 1;
				gridDimVirt[i] = gridDim[i] + 2;  // +2 for virtual grid
			}
			numBoxes = gridDim[0] * gridDim[1] * gridDim[2];
			numBoxesVirt = gridDimVirt[0] * gridDimVirt[1] * gridDimVirt[2];
		}

		/**
		 * Compute the linear index of a 3D box given its (ix, iy, iz) coordinates.
		 * 
		 * @param ix X-coordinate of the box (0-based).
		 * @param iy Y-coordinate of the box (0-based).
		 * @param iz Z-coordinate of the box (0-based).
		 * @param grid Array of grid dimensions {nx, ny, nz}.
		 * @returns Linear index in row-major order (x fastest): ix + iy * nx + iz * nx * ny.
		 */
		inline int get_box_by_index(int ix, int iy, int iz,
									const std::array<uint8_t, 3>& grid) const {
			return ix + iy * grid[0] + iz * grid[0] * grid[1];
		}
	};

} // namespace cpplib::geometry

template<class T>
struct std::hash<cpplib::geometry::Point<T>> {
	std::size_t operator()(const cpplib::geometry::Point<T>& s) const noexcept
	{
		std::size_t h1 = std::hash<T>{}(s[0]);
		std::size_t h2 = std::hash<T>{}(s[1]);
		std::size_t h3 = std::hash<T>{}(s[2]);
		return h1 ^ (h2 << 2) ^ (h3 << 4);
	}
};