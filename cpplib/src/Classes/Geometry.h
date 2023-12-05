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
#include <cstdint>
#include <type_traits>
#include <array>
#include <vector>
#include <cmath> // for std::sqrt and std::floor
#include <utility> //for std::move

namespace geometry {

	template<class T>
	class Point {
	public:
		using value_type = T;
		using array_type = std::array<T, 3>;
	private:
		array_type a_ = { 0,0,0 };

	public:
		// Constructors
		constexpr Point() noexcept = default;
		constexpr Point(T x, T y, T z) noexcept : a_{ x, y, z } {};
		explicit constexpr Point(const array_type& other) noexcept : a_(other) {};
		explicit constexpr Point(array_type&& other) noexcept : a_(std::move(other)) {};

		constexpr T r() const noexcept {
			return std::sqrt(a_[0] * a_[0] + a_[1] * a_[1] + a_[2] * a_[2]);
		}
		constexpr Point& MoveToCell() noexcept {
			a_[0] -= std::floor(a_[0]);
			a_[1] -= std::floor(a_[1]);
			a_[2] -= std::floor(a_[2]);
			return *this;
		}

		// Static constexpr functions
		static constexpr T Scalar(const Point& left, const Point& right) noexcept {
			return (left.a_[0] * right.a_[0] + left.a_[1] * right.a_[1] + left.a_[2] * right.a_[2]);
		}
		static constexpr Point Vector(const Point& left, const Point& right) noexcept {
			return Point(left.a_[1] * right.a_[2] - left.a_[2] * right.a_[1], left.a_[2] * right.a_[0] - left.a_[0] * right.a_[2], left.a_[0] * right.a_[1] - left.a_[1] * right.a_[0]);
		}

		//Operators
		constexpr Point operator-() const noexcept {
			return Point(-a_[0], -a_[1], -a_[2]);
		}
		constexpr Point operator+(const Point& right) const noexcept {
			return Point(a_[0] + right.a_[0], a_[1] + right.a_[1], a_[2] + right.a_[2]);
		}
		constexpr Point operator+(const T b) const noexcept {
			return Point(a_[0] + b, a_[1] + b, a_[2] + b);
		}
		constexpr Point operator-(const Point& right) const noexcept {
			return Point(a_[0] - right.a_[0], a_[1] - right.a_[1], a_[2] - right.a_[2]);
		}
		constexpr Point operator-(const T b) const noexcept {
			return Point(a_[0] - b, a_[1] - b, a_[2] - b);
		}
		constexpr Point operator*(const Point& right) const noexcept {
			return Point(a_[0] * right.a_[0], a_[1] * right.a_[1], a_[2] * right.a_[2]);
		}
		constexpr Point operator*(const T b) const noexcept {
			return Point(a_[0] * b, a_[1] * b, a_[2] * b);
		}
		constexpr Point operator/(const Point& right) const noexcept {
			return Point(a_[0] / right.a_[0], a_[1] / right.a_[1], a_[2] / right.a_[2]);
		}
		constexpr Point operator/(const T b) const noexcept {
			return Point(a_[0] / b, a_[1] / b, a_[2] / b);
		}
		Point& operator+=(const Point& right) noexcept {
			a_[0] += right.a_[0];
			a_[1] += right.a_[1];
			a_[2] += right.a_[2];
			return *this;
		}
		Point& operator+=(const T right) noexcept {
			a_[0] += right;
			a_[1] += right;
			a_[2] += right;
			return *this;
		}
		Point& operator-=(const Point& right) noexcept {
			a_[0] -= right.a_[0];
			a_[1] -= right.a_[1];
			a_[2] -= right.a_[2];
			return *this;
		}
		Point& operator-=(const T right) noexcept {
			a_[0] -= right;
			a_[1] -= right;
			a_[2] -= right;
			return *this;
		}
		Point& operator*=(const T right) noexcept {
			a_[0] *= right;
			a_[1] *= right;
			a_[2] *= right;
			return *this;
		}
		Point& operator/=(const T right) noexcept {
			a_[0] /= right;
			a_[1] /= right;
			a_[2] /= right;
			return *this;
		}
		constexpr bool operator==(const Point& other) const noexcept {
			return a_[0] == other.a_[0] && a_[1] == other.a_[1] && a_[2] == other.a_[2];
		}

		constexpr T get(const unsigned char i) const noexcept {
			return a_[i];
		}

		constexpr void set(const unsigned char i, const T newval) noexcept {
			a_[i] = newval;
		}
	};

	template<class T>
	class Matrix {
	public:
		using size_t = unsigned char;
		using value_type = T;
		using array_type = std::array< std::array<T, 3>, 3>;
		using const_array_type = const array_type;
	private:
		array_type A{ { {0,0,0},{0,0,0},{0,0,0} } };
		template<class T2> friend class Matrix; // for constructors
	public:
		constexpr Matrix() noexcept = default;

		constexpr Matrix(const Matrix<T>& r) noexcept {
			if (this != &r) A = r.A;
		}
		constexpr Matrix<T>& operator=(const Matrix<T>& r) noexcept {
			if (this != &r) A = r.A;
			return *this;
		}
		constexpr Matrix<T>& operator=(Matrix<T>&& r) noexcept {
			if (this != &r) A = std::move(r.A);
			return *this;
		}
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
		constexpr Matrix(Matrix<T>&& r) noexcept {
			if (this == &r) return;
			A = std::move(r.A);
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
	};
	template<class T1, class T2>
	Point<decltype(T1()* T2())> operator*(const Matrix<T1>& left, const Point<T2>& right) noexcept {
		Point<decltype(T1()* T2())> res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				res.set(i, res.get(i) + left.El(i, j) * right.get(j));
			}
		}
		return res;
	}

	template<class T>
	constexpr Matrix<T> EqualMatrix(static_cast<T>(1));

	template<class T>
	struct Cell {
	public:
		using value_type = T;
		using array_type = std::array<T, 3>;
		using matrix_type = geometry::Matrix<T>;
	private:
		static constexpr value_type GradtoRad = static_cast<value_type>(0.0174532925199432957692);
		static constexpr value_type RadtoGrad = static_cast<value_type>(57.295779513082320877);
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
				angleGrad_[i] = angleRad_[i] * RadtoGrad;
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
					angleRad_[i] = angleGrad_[i] * GradtoRad;
				}
			}
			else {
				angleRad_[0] = alpha;
				angleRad_[1] = beta;
				angleRad_[2] = gamma;
				for (int i = 0; i < 3; i++) {
					angleGrad_[i] = angleRad_[i] * RadtoGrad;
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

		template <class I>
		[[nodiscard]] constexpr geometry::Point<I> findOptimalSupercell(const value_type cutoff) const noexcept {
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
				static_cast<I>(std::ceil(cutoff / lattice_[0])),
				static_cast<I>(std::ceil(cutoff / lattice_[1])),
				static_cast<I>(std::ceil(cutoff / lattice_[2])));

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
					if (lattice_[1] * currentSuperCell.get(1) > lattice_[2] * currentSuperCell.get(2)) minDimention = 2;
					else minDimention = 1;
					break;
				case 2:
				case 3:
					if (lattice_[0] * currentSuperCell.get(0) > lattice_[2] * currentSuperCell.get(2)) minDimention = 2;
					else minDimention = 0;
					break;
				case 4:
				case 5:
					if (lattice_[0] * currentSuperCell.get(0) > lattice_[1] * currentSuperCell.get(1)) minDimention = 1;
					else minDimention = 0;
					break;
				default:
					if (lattice_[minDimention] * currentSuperCell.get(minDimention) > lattice_[1] * currentSuperCell.get(1)) minDimention = 1;
					if (lattice_[minDimention] * currentSuperCell.get(minDimention) > lattice_[2] * currentSuperCell.get(2)) minDimention = 2;
					break;
				}
				currentSuperCell.set(minDimention, currentSuperCell.get(minDimention) + 1);
			}
			return currentSuperCell;
		}
	};

	template<class T>
	struct Symm
	{
		using matrix_t = geometry::Matrix<int8_t>;
		using point_t = geometry::Point<T>;
		matrix_t mat;
		point_t point;
		uint8_t mult = 1;
		Symm() {}
		Symm(const Symm<T>& in) noexcept : mult(in.mult), point(in.point), mat(in.mat) {}
		Symm(Symm<T>&& in) noexcept : mult(std::move(in.mult)), point(std::move(in.point)), mat(std::move(in.mat)) {}
		Symm<T>& operator= (const Symm<T>& other) noexcept {
			if (&other != this) {
				mat = other.mat;
				point = other.point;
				mult = other.mult;
			}
			return *this;
		}
		Symm<T>& operator= (Symm<T>&& other) noexcept {
			if (&other != this) {
				mat = other.mat;
				point = other.point;
				mult = other.mult;
			}
			return *this;
		}
		Symm<T> MirrorSymm() const {
			Symm<T> out;
			out.mat = this->mat.Invert();
			out.point = -(this->point);
			out.mult = this->mult;
			return out;
		}
		void Multiplicity() {
			point_t c(static_cast<T>(0.1), static_cast<T>(0.15), static_cast<T>(0.2));
			point_t t(GenSymmNorm(c));
			for (mult = 1; (t - c).r() > (static_cast<T>(0.001)); mult++) {
				t = GenSymmNorm(t);
			}
		}
		explicit Symm(const char* str, bool isLATT = false) noexcept : mat(), point(static_cast<T>(0.0), static_cast<T>(0.0), static_cast<T>(0.0))
		{
			bool minus = false;
			bool correct = true;
			unsigned int j(0);
			for (unsigned int i = 0; str[i] != '\0'; i++) {
				switch (str[i]) {
				case 'x': [[fallthrough]];
				case 'X':
					if (minus == true) mat.El(j, 0) = -1;
					else mat.El(j, 0) = 1;
					minus = false;
					break;
				case 'y': [[fallthrough]];
				case 'Y':
					if (minus == true) mat.El(j, 1) = -1;
					else mat.El(j, 1) = 1;
					minus = false;
					break;
				case 'z': [[fallthrough]];
				case 'Z':
					if (minus == true) mat.El(j, 2) = -1;
					else mat.El(j, 2) = 1;
					minus = false;
					break;
				case ',':
					j++;
					[[fallthrough]];
				case ' ':
					[[fallthrough]];
				case '+':
					break;
				case '-':
					minus = true;
					break;
				default:
				{
					unsigned int k(i);
					for (; (str[i] >= '0' && str[i] <= '9') || str[i] == '.' || str[i] == '/'; i++);

					switch (i - k) {
					case 3:
						if (str[k] == '1') {
							switch (str[k + 2]) {
							case '2':
								point.set(j, static_cast<T>(0.5));
								break;
							case '3':
								if (minus == true) point.set(j, static_cast<T>(2.0 / 3.0));
								else point.set(j, static_cast<T>(1.0 / 3.0));
								break;
							case '4':
								if (minus == true) point.set(j, static_cast<T>(0.75));
								else point.set(j, static_cast<T>(0.25));
								break;
							case '6':
								if (minus == true) point.set(j, static_cast<T>(5.0 / 6.0));
								else point.set(j, static_cast<T>(1.0 / 6.0));
								break;
							}
						}
						else {
							switch (str[k + 2]) {
							case '3':
								if (minus == true) point.set(j, static_cast<T>(1.0 / 3.0));
								else point.set(j, static_cast<T>(2.0 / 3.0));
								break;
							case '4':
								if (minus == true) point.set(j, static_cast<T>(0.25));
								else point.set(j, static_cast<T>(0.75));
								break;
							case '5':
								point.set(j, static_cast<T>(0.5));
								break;
							case '6':
								if (minus == true) point.set(j, static_cast<T>(1.0 / 6.0));
								else point.set(j, static_cast<T>(5.0 / 6.0));
								break;
							}

						}
						break;
					case 4:
						if ((minus == true && str[k + 2] == '7') || (minus == false && str[k + 2] != '7')) point.set(j, static_cast<T>(0.25));
						else point.set(j, static_cast<T>(0.75));
						break;
					default:
						// SHOULD NOT BE HERE
						correct = false;
						break;
					}
					minus = false;
					i--;
				}
				}
			}
			if (correct == true)
				Multiplicity();
			else
				mult = 0;
		}
		point_t GenSymm(const point_t& in) const
		{
			return (point + (mat * in));
		}
		point_t GenSymmNorm(const point_t& in) const
		{
			point_t res = GenSymm(in);
			res.MoveToCell();
			return res;
		}
	};
}