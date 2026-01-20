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
#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>
#include <numbers>
#include <vector>
#include <utility>

#include <Geometry.h>
#include <BaseTypes.h>

using namespace std;
using namespace cpplib;
using namespace cpplib::geometry;
using namespace cpplib::basic_types;

// Common constants for tests
constexpr Point<FloatingPointType> TEST_POINT_A(0.01423f, 0.27322f, 0.01346f);
constexpr Point<FloatingPointType> TEST_POINT_B(1.0f, 2.0f, 3.0f);
constexpr Point<FloatingPointType> TEST_POINT_C(4.0f, 5.0f, 6.0f);
constexpr Matrix<FloatingPointType> TEST_MATRIX_M({10.4804f, -5.2402f, 0.f, 0.f, 9.07629264f, 0.f, 0.f, 0.f, 31.8116f});
constexpr Matrix<FloatingPointType> TEST_MATRIX_N({1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f});
constexpr FloatingPointType EPSILON = 0.00001;

// ==================== TESTS FOR Point ====================
TEST(PointTest, CreationAndOperations) {
	// Testing constructors
	ASSERT_NO_THROW({
		Point<FloatingPointType> p0;
		Point<FloatingPointType> pf(p0);
		Point<FloatingPointType> p1(0.01423f, 0.27322f, 0.01346f);
		Point<FloatingPointType> p2(0.0, 1.0, 2.0);
		EXPECT_NEAR(p0[0], 0.0, EPSILON);
		EXPECT_NEAR(p0[1], 0.0, EPSILON);
		EXPECT_NEAR(p0[2], 0.0, EPSILON);
					});

	// Testing point-matrix multiplication
	Point<FloatingPointType> res = TEST_MATRIX_M * TEST_POINT_A;
	EXPECT_NEAR(res[0], -1.28259146, EPSILON);
	EXPECT_NEAR(res[1], 2.47982478, EPSILON);
	EXPECT_NEAR(res[2], 0.428184122, EPSILON);

	// Testing r() method
	EXPECT_NEAR(TEST_POINT_A.r(), 0.273921, EPSILON);
}

TEST(PointTest, VectorOperations) {
	// Test scalar product
	FloatingPointType scalar = Point<FloatingPointType>::Scalar(TEST_POINT_B, TEST_POINT_C);
	EXPECT_NEAR(scalar, 32.0, EPSILON); // 1*4 + 2*5 + 3*6 = 32

	// Test vector product
	Point<FloatingPointType> cross = Point<FloatingPointType>::Vector(TEST_POINT_B, TEST_POINT_C);
	EXPECT_NEAR(cross[0], -3.0, EPSILON); // 2*6 - 3*5 = -3
	EXPECT_NEAR(cross[1], 6.0, EPSILON);  // 3*4 - 1*6 = 6
	EXPECT_NEAR(cross[2], -3.0, EPSILON); // 1*5 - 2*4 = -3

	// Test distance calculation
	FloatingPointType dist = Point<FloatingPointType>::distance(TEST_POINT_B, TEST_POINT_C);
	EXPECT_NEAR(dist, sqrt(27.0), EPSILON);
}

TEST(PointTest, AngleCalculations) {
	Point<FloatingPointType> a(1.0f, 0.0f, 0.0f);
	Point<FloatingPointType> b(0.0f, 0.0f, 0.0f);
	Point<FloatingPointType> c(0.0f, 1.0f, 0.0f);

	// Test angle calculation (90 degrees)
	FloatingPointType angle = Point<FloatingPointType>::angleRad(a, b, c);
	EXPECT_NEAR(angle, numbers::pi_v<FloatingPointType> / 2, EPSILON);

	// Test angle in degrees
	FloatingPointType angleDeg = Point<FloatingPointType>::angleGrad(a, b, c);
	EXPECT_NEAR(angleDeg, 90.0, EPSILON);
}

TEST(PointTest, TorsionSymmetry) {
	// Testing torsion symmetry properties
	Point<FloatingPointType> points[] = {
		{0.865301423f, 0.21727322f, 0.032461346f},
		{0.65630135423f, 0.23467322f, 0.835801346f},
		{0.35601423f, 0.346346227322f, 0.601346f},
		{0.65301423f, 0.2373334622f, 0.568501346f}
	};

	const auto& a1 = points[0];
	const auto& a2 = points[1];
	const auto& a3 = points[2];
	const auto& a4 = points[3];

	// Testing symmetry of torsion calculation
	EXPECT_NEAR(Point<FloatingPointType>::torsionRad(a1, a2, a3, a4),
				Point<FloatingPointType>::torsionRad(a4, a3, a2, a1), EPSILON);

	EXPECT_NEAR(Point<FloatingPointType>::torsionRad(a2, a1, a3, a4),
				Point<FloatingPointType>::torsionRad(a4, a3, a1, a2), EPSILON);
}

TEST(PointTest, ComparisonOperations) {
	Point<FloatingPointType> p1(1.0f, 2.0f, 3.0f);
	Point<FloatingPointType> p2(1.000001f, 2.000001f, 3.000001f);
	Point<FloatingPointType> p3(4.0f, 5.0f, 6.0f);

	EXPECT_TRUE(p1 == p2); // Within epsilon tolerance
	EXPECT_TRUE(p1 != p3);
	EXPECT_TRUE(p1 < p3);
	EXPECT_TRUE(p3 > p1);
}

TEST(PointTest, CellOperations) {
	Point<FloatingPointType> p(1.5f, 2.7f, -0.3f);

	// Test MoveToCell
	Point<FloatingPointType> p_cell = p;
	p_cell.MoveToCell();
	EXPECT_GE(p_cell[0], 0.0f);
	EXPECT_LT(p_cell[0], 1.0f);
	EXPECT_GE(p_cell[1], 0.0f);
	EXPECT_LT(p_cell[1], 1.0f);
	EXPECT_GE(p_cell[2], 0.0f);
	EXPECT_LT(p_cell[2], 1.0f);

	// Test round and floor
	Point<FloatingPointType> p_round = p.round();
	Point<FloatingPointType> p_floor = p.floor();
	EXPECT_NEAR(p_round[0], 2.0f, EPSILON);
	EXPECT_NEAR(p_round[1], 3.0f, EPSILON);
	EXPECT_NEAR(p_round[2], 0.0f, EPSILON);
	EXPECT_NEAR(p_floor[0], 1.0f, EPSILON);
	EXPECT_NEAR(p_floor[1], 2.0f, EPSILON);
	EXPECT_NEAR(p_floor[2], -1.0f, EPSILON);
}

// ==================== TESTS FOR Matrix ====================
TEST(MatrixTest, BasicOperations) {
	// Test constructor from array
	Matrix<FloatingPointType> m1(TEST_MATRIX_M);
	Matrix<FloatingPointType> m2(TEST_MATRIX_N);

	// Test element access
	EXPECT_NEAR(m1.El(0, 0), 10.4804f, EPSILON);
	EXPECT_NEAR(m1.El(1, 1), 9.07629264f, EPSILON);
	EXPECT_NEAR(m1.El(2, 2), 31.8116f, EPSILON);

	// Test determinant
	FloatingPointType det1 = m1.Det();
	FloatingPointType det2 = m2.Det();
	EXPECT_NEAR(det1, double(10.4804f) * double(9.07629264f) * double(31.8116f), EPSILON);
	EXPECT_NEAR(det2, 0.0f, EPSILON); // Matrix N is singular

	// Test trace
	FloatingPointType trace2 = m2.Trace();
	EXPECT_NEAR(trace2, 5.0f, EPSILON); // (1+5+9)/3 = 5
}

TEST(MatrixTest, MatrixOperations) {
	Matrix<FloatingPointType> m1(TEST_MATRIX_M);
	Matrix<FloatingPointType> identity(1.0f);

	// Test multiplication with identity
	auto m1_times_id = m1 * identity;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			EXPECT_NEAR(m1_times_id.El(i, j), m1.El(i, j), EPSILON);
		}
	}

	// Test transposition
	auto m1_transposed = m1.Transponate();
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			EXPECT_NEAR(m1_transposed.El(i, j), m1.El(j, i), EPSILON);
		}
	}

	// Test inversion (for non-singular matrix)
	auto m1_inverted = m1.Invert();
	auto product = m1 * m1_inverted;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (i == j) {
				EXPECT_NEAR(product.El(i, j), 1.0f, EPSILON);
			}
			else {
				EXPECT_NEAR(product.El(i, j), 0.0f, EPSILON);
			}
		}
	}
}

// ==================== TESTS FOR Plane ====================
TEST(PlaneTest, CreationAndOperations) {
	// Create plane from three points
	Point<FloatingPointType> p1(0.0f, 0.0f, 0.0f);
	Point<FloatingPointType> p2(1.0f, 0.0f, 0.0f);
	Point<FloatingPointType> p3(0.0f, 1.0f, 0.0f);

	Plane<FloatingPointType> plane1(p1, p2, p3);

	// All points should be in the plane
	EXPECT_NEAR(plane1.distance(p1), 0.0f, EPSILON);
	EXPECT_NEAR(plane1.distance(p2), 0.0f, EPSILON);
	EXPECT_NEAR(plane1.distance(p3), 0.0f, EPSILON);

	// Normal should point in positive z direction
	auto normal = plane1.normal();
	EXPECT_NEAR(normal[0], 0.0f, EPSILON);
	EXPECT_NEAR(normal[1], 0.0f, EPSILON);
	EXPECT_NEAR(abs(normal[2]), 1.0f, EPSILON);

	// Create plane from point and normal
	Point<FloatingPointType> point(1.0f, 2.0f, 3.0f);
	Point<FloatingPointType> norm(0.0f, 0.0f, 1.0f);
	Plane<FloatingPointType> plane2(point, norm);

	// Test projection
	Point<FloatingPointType> proj = plane2.make_projection(Point<FloatingPointType>(1.0f, 2.0f, 10.0f));
	EXPECT_NEAR(proj[0], 1.0f, EPSILON);
	EXPECT_NEAR(proj[1], 2.0f, EPSILON);
	EXPECT_NEAR(proj[2], 3.0f, EPSILON);

	// Test side calculation
	FloatingPointType side1 = plane2.side(Point<FloatingPointType>(1.0f, 2.0f, 5.0f));
	FloatingPointType side2 = plane2.side(Point<FloatingPointType>(1.0f, 2.0f, 1.0f));
	EXPECT_GT(side1, 0.0f);
	EXPECT_LT(side2, 0.0f);
}

// ==================== TESTS FOR Cell ====================
TEST(CellTest, CreationAndOperations) {
	// Create cubic cell
	Cell<FloatingPointType> cubicCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f);

	EXPECT_NEAR(cubicCell.lat_dir(0), 10.0f, EPSILON);
	EXPECT_NEAR(cubicCell.lat_dir(1), 10.0f, EPSILON);
	EXPECT_NEAR(cubicCell.lat_dir(2), 10.0f, EPSILON);
	EXPECT_NEAR(cubicCell.getAngleGrad(0), 90.0f, EPSILON);
	EXPECT_NEAR(cubicCell.getAngleGrad(1), 90.0f, EPSILON);
	EXPECT_NEAR(cubicCell.getAngleGrad(2), 90.0f, EPSILON);

	// Test fractional to cartesian conversion
	Point<FloatingPointType> fracPoint(0.5f, 0.5f, 0.5f);
	Point<FloatingPointType> cartPoint = cubicCell.fracToCart() * fracPoint;
	EXPECT_NEAR(cartPoint[0], 5.0f, EPSILON);
	EXPECT_NEAR(cartPoint[1], 5.0f, EPSILON);
	EXPECT_NEAR(cartPoint[2], 5.0f, EPSILON);

	// Test cartesian to fractional conversion
	Point<FloatingPointType> fracPoint2 = cubicCell.cartToFrac() * cartPoint;
	EXPECT_NEAR(fracPoint2[0], 0.5f, EPSILON);
	EXPECT_NEAR(fracPoint2[1], 0.5f, EPSILON);
	EXPECT_NEAR(fracPoint2[2], 0.5f, EPSILON);

	// Test distance calculation within unit cell
	Point<FloatingPointType> p1(0.1f, 0.1f, 0.1f);
	Point<FloatingPointType> p2(0.9f, 0.9f, 0.9f);
	FloatingPointType dist = cubicCell.distance_in_01(p1, p2);
	EXPECT_GT(dist, 0.0f);
	EXPECT_LT(dist, sqrt(3.0f) * 10.0f); // Maximum distance in cube
}

TEST(CellTest, TriclinicCell) {
	// Create triclinic cell
	Cell<FloatingPointType> triCell(10.0f, 12.0f, 15.0f, 80.0f, 85.0f, 75.0f);

	// Verify parameters
	EXPECT_NEAR(triCell.lat_dir(0), 10.0f, EPSILON);
	EXPECT_NEAR(triCell.lat_dir(1), 12.0f, EPSILON);
	EXPECT_NEAR(triCell.lat_dir(2), 15.0f, EPSILON);
	EXPECT_NEAR(triCell.getAngleGrad(0), 80.0f, EPSILON);
	EXPECT_NEAR(triCell.getAngleGrad(1), 85.0f, EPSILON);
	EXPECT_NEAR(triCell.getAngleGrad(2), 75.0f, EPSILON);

	// Verify that matrices are inverses of each other
	auto fracToCart = triCell.fracToCart();
	auto cartToFrac = triCell.cartToFrac();
	auto product = fracToCart * cartToFrac;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (i == j) {
				EXPECT_NEAR(product.El(i, j), 1.0f, EPSILON);
			}
			else {
				EXPECT_NEAR(product.El(i, j), 0.0f, EPSILON);
			}
		}
	}
}

// ==================== TESTS FOR Symm ====================
std::vector<const char*> SYMMETRY_STRINGS = {
"x, y, z"
,"-y, x-y, z"
,"-x+y, -x, z"
,"y, x, -z"
,"x-y, -y, -z"
,"-x, -x+y, -z"
,"-x, -y, -z"
,"y, -x+y, -z"
,"x-y, x, -z"
,"-y, -x, z"
,"-x+y, y, z"
,"x, x-y, z"
,"z, y, -x"
,"y, x, -z"
,"x, z, -y"
,"z, x, -y"
,"y, z, -x"
,"x, y, -z"
,"z, -y, x"
,"y, -x, z"
,"x, -z, y"
,"z, -x, y"
,"y, -z, x"
,"x, -y, z"
,"-z, y, x"
,"-y, x, z"
,"-x, z, y"
,"-z, x, y"
,"-y, z, x"
,"-x, y, z"
,"-z, -y, -x"
,"-y, -x, -z"
,"-x, -z, -y"
,"-z, -x, -y"
,"-y, -z, -x"
,"-x, -y, -z"
,"-z, -y, x"
,"-y, -x, z"
,"-x, -z, y"
,"-z, -x, y"
,"-y, -z, x"
,"-x, -y, z"
,"-z, y, -x"
,"-y, x, -z"
,"-x, z, -y"
,"-z, x, -y"
,"-y, z, -x"
,"-x, y, -z"
,"z, -y, -x"
,"y, -x, -z"
,"x, -z, -y"
,"z, -x, -y"
,"y, -z, -x"
,"x, -y, -z"
,"z, y, x"
,"y, x, z"
,"x, z, y"
,"z, x, y"
,"y, z, x"
,"x, y, z"
,"z, y+1/2, -x+1/2"
,"z+1/2, y, -x+1/2"
,"z+1/2, y+1/2, -x"
,"y, x+1/2, -z+1/2"
,"y+1/2, x, -z+1/2"
,"y+1/2, x+1/2, -z"
,"x, z+1/2, -y+1/2"
,"x+1/2, z, -y+1/2"
,"x+1/2, z+1/2, -y"
,"z, x+1/2, -y+1/2"
,"z+1/2, x, -y+1/2"
,"z+1/2, x+1/2, -y"
,"y, z+1/2, -x+1/2"
,"y+1/2, z, -x+1/2"
,"y+1/2, z+1/2, -x"
,"x, y+1/2, -z+1/2"
,"x+1/2, y, -z+1/2"
,"x+1/2, y+1/2, -z"
,"z, -y+1/2, x+1/2"
,"z+1/2, -y, x+1/2"
,"z+1/2, -y+1/2, x"
,"y, -x+1/2, z+1/2"
,"y+1/2, -x, z+1/2"
,"y+1/2, -x+1/2, z"
,"x, -z+1/2, y+1/2"
,"x+1/2, -z, y+1/2"
,"x+1/2, -z+1/2, y"
,"z, -x+1/2, y+1/2"
,"z+1/2, -x, y+1/2"
,"z+1/2, -x+1/2, y"
,"y, -z+1/2, x+1/2"
,"y+1/2, -z, x+1/2"
,"y+1/2, -z+1/2, x"
,"x, -y+1/2, z+1/2"
,"x+1/2, -y, z+1/2"
,"x+1/2, -y+1/2, z"
,"-z, y+1/2, x+1/2"
,"-z+1/2, y, x+1/2"
,"-z+1/2, y+1/2, x"
,"-y, x+1/2, z+1/2"
,"-y+1/2, x, z+1/2"
,"-y+1/2, x+1/2, z"
,"-x, z+1/2, y+1/2"
,"-x+1/2, z, y+1/2"
,"-x+1/2, z+1/2, y"
,"-z, x+1/2, y+1/2"
,"-z+1/2, x, y+1/2"
,"-z+1/2, x+1/2, y"
,"-y, z+1/2, x+1/2"
,"-y+1/2, z, x+1/2"
,"-y+1/2, z+1/2, x"
,"-x, y+1/2, z+1/2"
,"-x+1/2, y, z+1/2"
,"-x+1/2, y+1/2, z"
,"-z, -y+1/2, -x+1/2"
,"-z+1/2, -y, -x+1/2"
,"-z+1/2, -y+1/2, -x"
,"-y, -x+1/2, -z+1/2"
,"-y+1/2, -x, -z+1/2"
,"-y+1/2, -x+1/2, -z"
,"-x, -z+1/2, -y+1/2"
,"-x+1/2, -z, -y+1/2"
,"-x+1/2, -z+1/2, -y"
,"-z, -x+1/2, -y+1/2"
,"-z+1/2, -x, -y+1/2"
,"-z+1/2, -x+1/2, -y"
,"-y, -z+1/2, -x+1/2"
,"-y+1/2, -z, -x+1/2"
,"-y+1/2, -z+1/2, -x"
,"-x, -y+1/2, -z+1/2"
,"-x+1/2, -y, -z+1/2"
,"-x+1/2, -y+1/2, -z"
,"-z, -y+1/2, x+1/2"
,"-z+1/2, -y, x+1/2"
,"-z+1/2, -y+1/2, x"
,"-y, -x+1/2, z+1/2"
,"-y+1/2, -x, z+1/2"
,"-y+1/2, -x+1/2, z"
,"-x, -z+1/2, y+1/2"
,"-x+1/2, -z, y+1/2"
,"-x+1/2, -z+1/2, y"
,"-z, -x+1/2, y+1/2"
,"-z+1/2, -x, y+1/2"
,"-z+1/2, -x+1/2, y"
,"-y, -z+1/2, x+1/2"
,"-y+1/2, -z, x+1/2"
,"-y+1/2, -z+1/2, x"
,"-x, -y+1/2, z+1/2"
,"-x+1/2, -y, z+1/2"
,"-x+1/2, -y+1/2, z"
,"-z, y+1/2, -x+1/2"
,"-z+1/2, y, -x+1/2"
,"-z+1/2, y+1/2, -x"
,"-y, x+1/2, -z+1/2"
,"-y+1/2, x, -z+1/2"
,"-y+1/2, x+1/2, -z"
,"-x, z+1/2, -y+1/2"
,"-x+1/2, z, -y+1/2"
,"-x+1/2, z+1/2, -y"
,"-z, x+1/2, -y+1/2"
,"-z+1/2, x, -y+1/2"
,"-z+1/2, x+1/2, -y"
,"-y, z+1/2, -x+1/2"
,"-y+1/2, z, -x+1/2"
,"-y+1/2, z+1/2, -x"
,"-x, y+1/2, -z+1/2"
,"-x+1/2, y, -z+1/2"
,"-x+1/2, y+1/2, -z"
,"z, -y+1/2, -x+1/2"
,"z+1/2, -y, -x+1/2"
,"z+1/2, -y+1/2, -x"
,"y, -x+1/2, -z+1/2"
,"y+1/2, -x, -z+1/2"
,"y+1/2, -x+1/2, -z"
,"x, -z+1/2, -y+1/2"
,"x+1/2, -z, -y+1/2"
,"x+1/2, -z+1/2, -y"
,"z, -x+1/2, -y+1/2"
,"z+1/2, -x, -y+1/2"
,"z+1/2, -x+1/2, -y"
,"y, -z+1/2, -x+1/2"
,"y+1/2, -z, -x+1/2"
,"y+1/2, -z+1/2, -x"
,"x, -y+1/2, -z+1/2"
,"x+1/2, -y, -z+1/2"
,"x+1/2, -y+1/2, -z"
,"z, y+1/2, x+1/2"
,"z+1/2, y, x+1/2"
,"z+1/2, y+1/2, x"
,"y, x+1/2, z+1/2"
,"y+1/2, x, z+1/2"
,"y+1/2, x+1/2, z"
,"x, z+1/2, y+1/2"
,"x+1/2, z, y+1/2"
,"x+1/2, z+1/2, y"
,"z, x+1/2, y+1/2"
,"z+1/2, x, y+1/2"
,"z+1/2, x+1/2, y"
,"y, z+1/2, x+1/2"
,"y+1/2, z, x+1/2"
,"y+1/2, z+1/2, x"
,"x, y+1/2, z+1/2"
,"x+1/2, y, z+1/2"
,"x+1/2, y+1/2, z"
,"x-y, -y, -z+1/2"
,"-x, -x+y, -z+1/2"
,"y, x, -z+1/2"
,"x-y, x, -z"
,"y, -x+y, -z"
,"-x, -y, -z"
,"-x+y, y, z+1/2"
,"x, x-y, z+1/2"
,"-y, -x, z+1/2"
,"-x+y, -x, z"
,"-y, x-y, z"
,"x, y, z"
,"x-y+2/3, -y+1/3, -z+5/6"
,"x-y+1/3, -y+2/3, -z+1/6"
,"-x+2/3, -x+y+1/3, -z+5/6"
,"-x+1/3, -x+y+2/3, -z+1/6"
,"y+2/3, x+1/3, -z+5/6"
,"y+1/3, x+2/3, -z+1/6"
,"x-y+2/3, x+1/3, -z+1/3"
,"x-y+1/3, x+2/3, -z+2/3"
,"y+2/3, -x+y+1/3, -z+1/3"
,"y+1/3, -x+y+2/3, -z+2/3"
,"-x+2/3, -y+1/3, -z+1/3"
,"-x+1/3, -y+2/3, -z+2/3"
,"-x+y+2/3, y+1/3, z+5/6"
,"-x+y+1/3, y+2/3, z+1/6"
,"x+2/3, x-y+1/3, z+5/6"
,"x+1/3, x-y+2/3, z+1/6"
,"-y+2/3, -x+1/3, z+5/6"
,"-y+1/3, -x+2/3, z+1/6"
,"-x+y+2/3, -x+1/3, z+1/3"
,"-x+y+1/3, -x+2/3, z+2/3"
,"-y+2/3, x-y+1/3, z+1/3"
,"-y+1/3, x-y+2/3, z+2/3"
,"x+2/3, y+1/3, z+1/3"
,"x+1/3, y+2/3, z+2/3"
};

TEST(SymmTest, CreationAndParsing) {
	// Testing creation of symmetry objects from strings
	for (const auto& str : SYMMETRY_STRINGS) {
		Symm<FloatingPointType> symm(str);
		EXPECT_TRUE(symm.isValid());
	}
}

TEST(SymmTest, SymmetryOperations) {
	// Test identity symmetry
	Symm<FloatingPointType> identity("x, y, z");
	Point<FloatingPointType> testPoint(0.1f, 0.2f, 0.3f);

	// Identity should return the same point
	Point<FloatingPointType> result = identity.GenSymm(testPoint);
	EXPECT_NEAR(result[0], 0.1f, EPSILON);
	EXPECT_NEAR(result[1], 0.2f, EPSILON);
	EXPECT_NEAR(result[2], 0.3f, EPSILON);

	// Test with normalization
	Point<FloatingPointType> normResult = identity.GenSymmNorm(testPoint);
	EXPECT_NEAR(normResult[0], 0.1f, EPSILON);
	EXPECT_NEAR(normResult[1], 0.2f, EPSILON);
	EXPECT_NEAR(normResult[2], 0.3f, EPSILON);

	// Test mirror symmetry
	Symm<FloatingPointType> mirrorX("-x, y, z");
	result = mirrorX.GenSymm(testPoint);
	EXPECT_NEAR(result[0], -0.1f, EPSILON);
	EXPECT_NEAR(result[1], 0.2f, EPSILON);
	EXPECT_NEAR(result[2], 0.3f, EPSILON);

	// Test with shift
	Symm<FloatingPointType> shifted("x+1/2, y+1/2, z+1/2");
	result = shifted.GenSymm(testPoint);
	EXPECT_NEAR(result[0], 0.6f, EPSILON);
	EXPECT_NEAR(result[1], 0.7f, EPSILON);
	EXPECT_NEAR(result[2], 0.8f, EPSILON);

	// Test normalized result with shift
	normResult = shifted.GenSymmNorm(testPoint);
	EXPECT_GE(normResult[0], 0.0f);
	EXPECT_LT(normResult[0], 1.0f);
	EXPECT_GE(normResult[1], 0.0f);
	EXPECT_LT(normResult[1], 1.0f);
	EXPECT_GE(normResult[2], 0.0f);
	EXPECT_LT(normResult[2], 1.0f);
}

TEST(SymmTest, MirrorSymmetry) {
	Symm<FloatingPointType> symm("x+1/2, -y, z+1/3");
	Symm<FloatingPointType> mirror = symm.MirrorSymm();

	// Test that applying symm then its mirror gives identity (within cell)
	Point<FloatingPointType> testPoint(0.25f, 0.5f, 0.75f);
	Point<FloatingPointType> transformed = symm.GenSymm(testPoint);
	Point<FloatingPointType> backTransformed = mirror.GenSymmNorm(transformed);

	// Should get back to original point (within unit cell)
	EXPECT_NEAR(backTransformed[0], testPoint[0], EPSILON);
	EXPECT_NEAR(backTransformed[1], testPoint[1], EPSILON);
	EXPECT_NEAR(backTransformed[2], testPoint[2], EPSILON);
}

TEST(SymmTest, ValidityChecks) {
	// Valid symmetries
	Symm<FloatingPointType> valid1("x, y, z");
	Symm<FloatingPointType> valid2("-y, x-y, z");
	Symm<FloatingPointType> valid3("x+1/2, y+1/2, z+1/2");

	EXPECT_TRUE(valid1.isValid());
	EXPECT_TRUE(valid2.isValid());
	EXPECT_TRUE(valid3.isValid());

	// Check identity detection
	EXPECT_TRUE(valid1.is_Eq());
	EXPECT_FALSE(valid2.is_Eq());
	EXPECT_FALSE(valid3.is_Eq());
}

// ==================== TESTS FOR Polygon ====================
class PolygonTest : public ::testing::Test {
protected:
	using PointType = Point<FloatingPointType>;

	void SetUp() override {
		// Standard polygons for testing
		square = {PointType(0,0,0), PointType(1,0,0),
				   PointType(1,1,0), PointType(0,1,0)};

		triangle = {PointType(0,0,0), PointType(1,0,0), PointType(0,1,0)};

		concaveShape = {PointType(0,0,0), PointType(2,0,0),
						 PointType(1,1,0), PointType(2,2,0), PointType(0,2,0)};
	}

	// Helper method to validate polygon convexity
	void validateConvexity(const Polygon<FloatingPointType>& poly, bool expectedConvex) {
		EXPECT_EQ(poly.isConvex(), expectedConvex);

		if (poly.size() >= 3 && expectedConvex) {
			// Additional validation for convex polygons
			for (size_t i = 0; i < poly.size(); ++i) {
				size_t j = (i + 1) % poly.size();
				size_t k = (i + 2) % poly.size();

				PointType edge1 = poly[j] - poly[i];
				PointType edge2 = poly[k] - poly[j];
				PointType cross = PointType::Vector(edge1, edge2);

				// All cross products should have the same sign for convex polygons
				EXPECT_GT(cross[2], -EPSILON); // For 2D polygons in XY plane
			}
		}
	}

	std::vector<PointType> square;
	std::vector<PointType> triangle;
	std::vector<PointType> concaveShape;
};

TEST_F(PolygonTest, CreationFromPoints) {
	// Test square creation
	Polygon<FloatingPointType> polySquare(square);
	ASSERT_EQ(polySquare.size(), 4);
	validateConvexity(polySquare, true);

	// Test triangle creation
	Polygon<FloatingPointType> polyTriangle(triangle);
	ASSERT_EQ(polyTriangle.size(), 3);
	validateConvexity(polyTriangle, true);

	// Test concave shape creation
	Polygon<FloatingPointType> polyConcave(concaveShape);
	validateConvexity(polyConcave, false);
}

TEST_F(PolygonTest, CreationFromPlaneAndCenter) {
	// Test polygon creation from plane and center
	PointType center(0.5, 0.5, 0);
	Plane<FloatingPointType> plane(center, PointType(0, 0, 1));

	Polygon<FloatingPointType> poly(plane, center, 1.0);
	ASSERT_EQ(poly.size(), 4);

	// All points should lie in the plane
	for (size_t i = 0; i < poly.size(); ++i) {
		EXPECT_NEAR(plane.distance(poly[i]), 0.0, EPSILON);
	}

	validateConvexity(poly, true);
}

TEST_F(PolygonTest, AccessOperations) {
	Polygon<FloatingPointType> poly(square);

	// Test element access
	EXPECT_EQ(poly[0], square[0]);
	EXPECT_EQ(poly[1], square[1]);
	EXPECT_EQ(poly[2], square[2]);
	EXPECT_EQ(poly[3], square[3]);
}

TEST_F(PolygonTest, ClippingSimple) {
	Polygon<FloatingPointType> polySquare(square);
	Plane<FloatingPointType> plane1(PointType(0.5, 0, 0), PointType(1, 0, 0));
	polySquare.clipByPlane(plane1);

	for (size_t i = 0; i < polySquare.size(); ++i) {
		EXPECT_GE(polySquare[i][0], 0.5 - EPSILON);
	}
	validateConvexity(polySquare, true);
}
TEST_F(PolygonTest, ClippingOrderProblem) {
	Polygon<FloatingPointType> polySquare1(square);
	Plane<FloatingPointType> plane1(PointType(0.5, 0.5, 0), PointType(-1, -1, 0));
	polySquare1.clipByPlane(plane1);

	EXPECT_EQ(polySquare1.size(), 3);
	validateConvexity(polySquare1, true);

	Polygon<FloatingPointType> polySquare2(square);
	Plane<FloatingPointType> plane2(PointType(0.0, 0.0, 0), PointType(-1, -1, 0));
	polySquare2.clipByPlane(plane2);

	EXPECT_EQ(polySquare2.size(), 0);
	validateConvexity(polySquare2, false);

	Polygon<FloatingPointType> polySquare3(square);
	Plane<FloatingPointType> plane3(PointType(0.0, 0.0, 0), PointType(1, 1, 0));
	polySquare3.clipByPlane(plane3);

	EXPECT_EQ(polySquare3.size(), 4);
	validateConvexity(polySquare3, true);

}
TEST_F(PolygonTest, ClippingCorrectSides) {
	Polygon<FloatingPointType> polySquare5(square);
	Plane<FloatingPointType> plane5(PointType(0.75, 0.75, 0), PointType(-1, -1, 0));
	polySquare5.clipByPlane(plane5);

	EXPECT_EQ(polySquare5.size(), 5);
	validateConvexity(polySquare5, true);

	Polygon<FloatingPointType> polySquare3(square);
	Plane<FloatingPointType> plane3(PointType(0.75, 0.75, 0), PointType(1, 1, 0));
	polySquare3.clipByPlane(plane3);

	EXPECT_EQ(polySquare3.size(), 3);
	validateConvexity(polySquare3, true);
}

TEST_F(PolygonTest, EdgeCases) {
	// Insufficient number of points
	std::vector<PointType> insufficient = {PointType(0,0,0), PointType(1,0,0)};
	Polygon<FloatingPointType> smallPoly(insufficient);
	EXPECT_FALSE(smallPoly.isConvex());

	// Empty polygon
	Polygon<FloatingPointType> emptyPoly;
	EXPECT_EQ(emptyPoly.size(), 0);
}

//// ==================== TESTS FOR Voronoi ====================
//class VoronoiTest : public ::testing::Test {
//protected:
//	using PointType = Point<FloatingPointType>;
//
//    void SetUp() override {
//        // Common points for Voronoi tests
//        commonPoints = {
//            PointType(0.1, 0.1, 0.1),
//            PointType(0.4, 0.4, 0.4),
//            PointType(0.254, 0.4, 0.364),
//            PointType(0.954, 0.866, 0.23),
//            PointType(0.7, 0.7, 0.7)
//        };
//        singlePoint = {
//            PointType(0.0, 0.0, 0.0)
//        };
//    }
//
//	// Helper method to validate Voronoi cell properties
//	void validateVoronoiCell(const VoronoiCell<FloatingPointType>& cell) {
//		const auto& seed = cell.getSeed();
//		const auto& faces = cell.getFaces();
//
//        for (const auto& face : faces) {
//            EXPECT_TRUE(face.poly.isConvex());
//
//			// Verify that seed is on the correct side of the face
//			VoronoiCell<FloatingPointType>::Face::PlaneType plane(face.poly[0], face.poly[1], face.poly[2]);
//			EXPECT_GE(plane.side(seed), 0.0 - EPSILON);
//		}
//	}
//
//    std::vector<PointType> commonPoints;
//    std::vector<PointType> singlePoint;
//};
//
//
//TEST_F(VoronoiTest, AlexTest) {
//    cpplib::geometry::Cell cell (10.9815, 6.8214, 8.9974, 90.0, 101.511, 90.0);
//    std::vector<const char*> symms{"x, y, z", "-x, y+1/2, -z+1/2", "-x, -y, -z", "x, -y-1/2, z-1/2"};
//
//
//
//
//
//    cpplib::geometry::VoronoiDiagram<double>::PointVector data = {{0.85991, 0.6099,  0.54993},
//                                                                  {0.9802,  0.2973,  0.68819},
//                                                                  {0.45147, 0.2722, -0.01408},
//                                                                  {0.4257,  0.369,  -0.0692 },
//                                                                  {0.54605, 0.33304, 0.10639},
//                                                                  {0.59711, 0.18811, 0.18507},
//                                                                  {0.5699,  0.0585,  0.1577 },
//                                                                  {0.6972,  0.22032, 0.31696},
//                                                                  {0.73155, 0.40837, 0.37126},
//                                                                  {0.6897,  0.5202,  0.323  },
//                                                                  {0.8262,  0.42976, 0.49478},
//                                                                  {0.88775, 0.26835, 0.56724},
//                                                                  {0.85479, 0.08196, 0.51712},
//                                                                  {0.8965, -0.0288,  0.5674 },
//                                                                  {0.75907, 0.05834, 0.39111},
//                                                                  {0.7356, -0.0699,  0.3552 }};
//
//
//
//
//    cpplib::geometry::VoronoiDiagram<double>::BoolVector bools{false, false, false, false, 
//                                                               false, true,  false, true, 
//                                                               false, false, false, false, 
//                                                               false, false, false, false};
//    VoronoiDiagram<FloatingPointType> vd(data, bools);
//    HashedSpace<FloatingPointType, AtomIndex> hs(cell, 6.0);
//    auto bonds = hs.create_hash_bonds<std::pair<AtomIndex, AtomIndex>>(data);
//    vd.calculateFaces(bonds);
//    auto cells = vd.extractCells();
//    for (const auto& c : cells) {
//        validateVoronoiCell(c);
//    }
//    VoronoiFused<double> vf;
//    vf.AddCells(cells);
//}
//
//TEST_F(VoronoiTest, CellConstruction) {
//	// Default constructor
//	VoronoiCell<FloatingPointType> defaultCell;
//	EXPECT_EQ(defaultCell.getSeed(), PointType(0, 0, 0));
//	validateVoronoiCell(defaultCell);
//
//	// Constructor with specified seed
//	PointType seed(1.0, 2.0, 3.0);
//	VoronoiCell<FloatingPointType> customCell(seed);
//	EXPECT_EQ(customCell.getSeed(), seed);
//	validateVoronoiCell(customCell);
//}
//
//TEST_F(VoronoiTest, CellInteraction) {
//	// Normal cell interaction
//	VoronoiCell<FloatingPointType> cell1(PointType(0.1, 0.1, 0.1));
//	VoronoiCell<FloatingPointType> cell2(PointType(0.4, 0.4, 0.4));
//
//	int result = VoronoiCell<FloatingPointType>::interact(cell1, cell2);
//	EXPECT_EQ(result, 0); // Normal interaction
//
//	// Cells with close seeds
//	VoronoiCell<FloatingPointType> cell3(PointType(0.1, 0.1, 0.1));
//	VoronoiCell<FloatingPointType> cell4(PointType(0.1000001, 0.1000001, 0.1000001));
//
//	result = VoronoiCell<FloatingPointType>::interact(cell3, cell4);
//	EXPECT_EQ(result, 1); // Too close seeds
//}
//
//TEST_F(VoronoiTest, DiagramConstructionCommon) {
//    // Basic constructor
//    VoronoiDiagram<FloatingPointType> vd(commonPoints);
//    auto cells = vd.extractCells();
//
//	EXPECT_EQ(cells.size(), commonPoints.size());
//	for (const auto& cell : cells) {
//		validateVoronoiCell(cell);
//	}
//
//    // Constructor with flags
//    std::vector<bool> flags(commonPoints.size(), true);
//    VoronoiDiagram<FloatingPointType> vdWithFlags(commonPoints, flags);
//
//	cells = vdWithFlags.extractCells();
//	EXPECT_EQ(cells.size(), commonPoints.size());
//}
//TEST_F(VoronoiTest, DiagramConstructionSingle) {
//    // Basic constructor
//    VoronoiDiagram<FloatingPointType> vd(singlePoint);
//    auto cells = vd.extractCells();
//
//    EXPECT_EQ(cells.size(), singlePoint.size());
//    for (const auto& cell : cells) {
//        validateVoronoiCell(cell);
//    }
//
//    // Constructor with flags
//    std::vector<bool> flags(singlePoint.size(), true);
//    VoronoiDiagram<FloatingPointType> vdWithFlags(singlePoint, flags);
//
//    cells = vdWithFlags.extractCells();
//    EXPECT_EQ(cells.size(), singlePoint.size());
//}
//
//TEST_F(VoronoiTest, DiagramOperations) {
//	VoronoiDiagram<FloatingPointType> vd;
//
//	// Adding points to diagram
//	std::vector<PointType> points = {PointType(0.2, 0.2, 0.2), PointType(0.5, 0.5, 0.5)};
//	std::vector<bool> flags(points.size(), true);
//	vd.addPoints(points, flags);
//
//	auto cells = vd.extractCells();
//	EXPECT_EQ(cells.size(), points.size());
//
//	// Calculating faces using bond list
//	std::vector<std::pair<int, int>> bondlist = {
//		{0,1}, {0,2}, {0,3}, {0,4},
//		{1,2}, {1,3}, {1,4},
//		{2,3}, {2,4},
//		{3,4}
//	};
//
//	vd.calculateFaces<int>(bondlist);
//	cells = vd.extractCells();
//
//	for (const auto& cell : cells) {
//		validateVoronoiCell(cell);
//	}
//}
//
//TEST_F(VoronoiTest, SpecializedOperations) {
//	// Calculating longest diagonal
//	std::vector<PointType> points = {PointType(0.1, 0.1, 0.1), PointType(0.9, 0.9, 0.9)};
//	VoronoiDiagram<FloatingPointType> vd(points);
//
//    Cell<FloatingPointType> cell(10, 10, 10, 90, 90, 120, true);
//    FloatingPointType diagonal = vd.calculateLongestDiagonal(cell.fracToCart());
//    EXPECT_GT(diagonal, 0.0);
//    HashedSpace<FloatingPointType,AtomIndex> hs(cell, diagonal);
//    auto bonds = hs.create_hash_bonds<std::pair<AtomIndex,AtomIndex>>(points);
//    vd.calculateFaces(bonds);
//    
//
//	// Working with VoronoiFused
//	auto cells = vd.extractCells();
//	VoronoiFused<FloatingPointType> vf;
//	vf.AddCells(cells);
//
//	// Verifying vertex uniqueness
//	for (size_t i = 0; i < vf.vertexes.size(); ++i) {
//		for (size_t j = i + 1; j < vf.vertexes.size(); ++j) {
//			EXPECT_GT((vf.vertexes[i] - vf.vertexes[j]).r(), 0.0001);
//		}
//	}
//
//	EXPECT_EQ(vf.centers.size(), cells.size());
//}
//
//TEST_F(VoronoiTest, BoundaryConditions) {
//	// Testing cell interactions at boundaries
//	VoronoiCell<FloatingPointType> cell1(PointType(0.0, 0.0, 0.0));
//	VoronoiCell<FloatingPointType> cell2(PointType(0.999, 0.999, 0.999));
//
//	int result = VoronoiCell<FloatingPointType>::interact(cell1, cell2);
//	EXPECT_EQ(result, 0); // Should interact normally
//}

// ==================== ENTRY POINT ====================
int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}