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


#include "../src/BaseHeaders/BaseTypes.h"
#include "../src/Classes/Cluster.h"
#include "../src/Classes/Geometry.h"
#include "../src/Classes/Voronoi.h"

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
	EXPECT_NEAR(trace2, 15.0f, EPSILON); // (1+5+9) = 15
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

// ==================== TESTS FOR Voronoi ====================
TEST(VoronoiTest, VertexEquality) {
	using namespace cpplib::voronoi;

	// Test vertex equality with epsilon
	Vertex v1(0, Point<FloatingPointType>(1.0, 2.0, 3.0));
	Vertex v2(1, Point<FloatingPointType>(1.0 + voronoi::EPSILON / 2,
										  2.0 + voronoi::EPSILON / 2,
										  3.0 + voronoi::EPSILON / 2));
	Vertex v3(2, Point<FloatingPointType>(1.0 + voronoi::EPSILON * 2,
										  2.0, 3.0));

	EXPECT_TRUE(v1 == v2);  // Within epsilon
	EXPECT_FALSE(v1 == v3); // Outside epsilon
}

TEST(VoronoiTest, VertexDistanceCalculation) {
	using namespace cpplib::voronoi;

	voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);
	geometry::Matrix<FloatingPointType> fractocart(1.0); // Identity

	const auto* maxVert = cell.update_vertices_distances(fractocart);

	EXPECT_NE(maxVert, nullptr);
	EXPECT_GT(maxVert->distance, 0.0);

	// Verify all vertices have computed distances
	for (const auto& vert : cell.vertices) {
		if (vert->get_state() != voronoi::State::DELETE) {
			EXPECT_GT(vert->distance, 0.0);
		}
	}
}

TEST(VoronoiTest, EdgeStateCalculation) {
	using namespace cpplib::voronoi;

	// Create vertices
	auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
	auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));

	Edge edge(0, v1.get(), v2.get());

	// Initially, edge should be in INVALID state (no faces)
	EXPECT_EQ(edge.calculateState(), State::INVALID);

	// Create faces
	auto f1 = std::make_unique<Face>(0);
	auto f2 = std::make_unique<Face>(1);

	// Add edge to faces
	f1->edges.insert(&edge);
	f2->edges.insert(&edge);
	edge.faces.insert(f1.get());
	edge.faces.insert(f2.get());

	// Now edge should be VALID (2 vertices, 2 faces, both vertices valid)
	EXPECT_EQ(edge.calculateState(), State::VALID);

	// Delete one vertex
	v1->set_state(State::DELETE);
	EXPECT_EQ(edge.calculateState(), State::MODIFICATION);

	// Delete both vertices
	v2->set_state(State::DELETE);
	EXPECT_EQ(edge.calculateState(), State::DELETE);
}

TEST(VoronoiTest, EdgeIntersection) {
	using namespace cpplib::voronoi;

	auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
	auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));
	v1->set_state(State::VALID);
	v2->set_state(State::DELETE);

	Edge edge(0, v1.get(), v2.get());
	auto f1 = std::make_unique<Face>(0);
	auto f2 = std::make_unique<Face>(1);
	f1->edges.insert(&edge);
	f2->edges.insert(&edge);
	edge.faces.insert(f1.get());
	edge.faces.insert(f2.get());

	// Edge in MODIFICATION state
	EXPECT_EQ(edge.calculateState(), State::MODIFICATION);

	// Create a plane cutting the edge at x=0.5
	geometry::Plane<FloatingPointType> plane(
		Point<FloatingPointType>(0.5, 0.0, 0.0),
		Point<FloatingPointType>(1.0, 0.0, 0.0)
	);

	auto intersection = edge.intersectSegmentPlane(plane);

	EXPECT_NEAR(intersection[0], 0.5, voronoi::EPSILON);
	EXPECT_NEAR(intersection[1], 0.0, voronoi::EPSILON);
	EXPECT_NEAR(intersection[2], 0.0, voronoi::EPSILON);
}

TEST(VoronoiTest, FaceStateCalculation) {
	using namespace cpplib::voronoi;

	Face face(0);

	// Create vertices
	auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
	auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));
	auto v3 = std::make_unique<Vertex>(2, Point<FloatingPointType>(1.0, 1.0, 0.0));

	// Create edges
	Edge e1(0, v1.get(), v2.get());
	Edge e2(1, v2.get(), v3.get());
	Edge e3(2, v3.get(), v1.get());

	face.edges.insert({&e1, &e2, &e3});
	face.vertices.insert({v1.get(), v2.get(), v3.get()});

	// All edges VALID -> face VALID
	e1.set_state(State::VALID);
	e2.set_state(State::VALID);
	e3.set_state(State::VALID);
	EXPECT_EQ(face.calculateState(), State::VALID);

	// One edge MODIFICATION -> face MODIFICATION
	e1.set_state(State::MODIFICATION);
	EXPECT_EQ(face.calculateState(), State::MODIFICATION);

	// All edges DELETE -> face DELETE
	e1.set_state(State::DELETE);
	e2.set_state(State::DELETE);
	e3.set_state(State::DELETE);
	EXPECT_EQ(face.calculateState(), State::DELETE);
}

TEST(VoronoiTest, CellClippingParallelPlane) {
	using namespace cpplib::voronoi;

	voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);

	size_t initialVertexCount = cell.vertices.size();
	size_t initialEdgeCount = cell.edges.size();
	size_t initialFaceCount = cell.faces.size();

	// Plane parallel to a face of the cube (shouldn't clip anything)
	geometry::Plane<FloatingPointType> parallelPlane(
		Point<FloatingPointType>(2.0, 0.5, 0.5),
		Point<FloatingPointType>(-1.0, 0.0, 0.0)
	);

	cell.clipByPlaneAndAddNewFace(parallelPlane, 1, ShiftCode(14));

	// Should not change the cell
	EXPECT_EQ(cell.vertices.size(), initialVertexCount);
	EXPECT_EQ(cell.edges.size(), initialEdgeCount);
	EXPECT_EQ(cell.faces.size(), initialFaceCount);
}

TEST(VoronoiTest, CellClippingTangent) {
	using namespace cpplib::voronoi;

	voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);

	// Plane tangent to the cube (just touches an edge)
	geometry::Plane<FloatingPointType> tangentPlane(
		Point<FloatingPointType>(1.0, 0.5, 0.5),
		Point<FloatingPointType>(-1.0, 0.0, 0.0)
	);

	size_t initialVertexCount = cell.vertices.size();

	cell.clipByPlaneAndAddNewFace(tangentPlane, 1, ShiftCode(14));

	// Should add some vertices at the tangent points
	EXPECT_GE(cell.vertices.size(), initialVertexCount);
}

TEST(VoronoiTest, VoronoiFusedVertexMerging) {
	using namespace cpplib::voronoi;

	// Create two cells with overlapping vertices
	std::vector<voronoi::Cell > cells;
	cells.emplace_back(Point<FloatingPointType>(0.3, 0.5, 0.5), 0);
	cells.emplace_back(Point<FloatingPointType>(0.7, 0.5, 0.5), 1);

	// Clip both cells so they share a face
	geometry::Plane<FloatingPointType> plane1(
		Point<FloatingPointType>(0.5, 0.5, 0.5),
		Point<FloatingPointType>(-1.0, 0.0, 0.0)
	);
	cells[0].clipByPlaneAndAddNewFace(plane1, 1, ShiftCode());

	geometry::Plane<FloatingPointType> plane2(
		Point<FloatingPointType>(0.5, 0.5, 0.5),
		Point<FloatingPointType>(1.0, 0.0, 0.0)
	);
	cells[1].clipByPlaneAndAddNewFace(plane2, 0, ShiftCode());

	geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
	VoronoiFused fused(cells, cell.fracToCart());

	// Verify that duplicate vertices were merged
	for (size_t i = 0; i < fused.vertices.size(); ++i) {
		for (size_t j = i + 1; j < fused.vertices.size(); ++j) {
			FloatingPointType dist = (fused.vertices[i] - fused.vertices[j]).r();
			EXPECT_GT(dist, voronoi::EPSILON);
		}
	}

	// Verify polyhedra have correct vertex references
	EXPECT_GT(fused.polyhedra[0].vert_ids.size(), 0);
	EXPECT_GT(fused.polyhedra[1].vert_ids.size(), 0);
}

TEST(VoronoiTest, VoronoiFusedEdgePolygonPopulation) {
	using namespace cpplib::voronoi;

	// Use the simple two-cell setup
	std::vector<voronoi::Cell> cells;
	cells.emplace_back(Point<FloatingPointType>(0.3, 0.5, 0.5), 0);
	cells.emplace_back(Point<FloatingPointType>(0.7, 0.5, 0.5), 1);

	geometry::Plane<FloatingPointType> plane1(
		Point<FloatingPointType>(0.5, 0.5, 0.5),
		Point<FloatingPointType>(-1.0, 0.0, 0.0)
	);
	cells[0].clipByPlaneAndAddNewFace(plane1, 1, ShiftCode());

	geometry::Plane<FloatingPointType> plane2(
		Point<FloatingPointType>(0.5, 0.5, 0.5),
		Point<FloatingPointType>(1.0, 0.0, 0.0)
	);
	cells[1].clipByPlaneAndAddNewFace(plane2, 0, ShiftCode());

	geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
	VoronoiFused fused(cells, cell.fracToCart());

	// Verify that polyhedra have edge_ids and poly_ids populated
	for (size_t i = 0; i < fused.polyhedra.size(); ++i) {
		EXPECT_GT(fused.polyhedra[i].edge_ids.size(), 0)
			<< "Polyhedron " << i << " should have edges";
		EXPECT_GT(fused.polyhedra[i].poly_ids.size(), 0)
			<< "Polyhedron " << i << " should have polygons";

		// Verify all edge IDs are valid
		for (auto edge_id : fused.polyhedra[i].edge_ids) {
			EXPECT_LT(edge_id, fused.edges.size());
		}

		// Verify all polygon IDs are valid
		for (auto poly_id : fused.polyhedra[i].poly_ids) {
			EXPECT_LT(poly_id, fused.polygons.size());
		}
	}
}

TEST(VoronoiTest, DegenerateCases) {
	using namespace cpplib::voronoi;

	// Test with points very close together
	geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

	voronoi::VoronoiDiagram::PointVector closePoints = {
		Point<FloatingPointType>(0.5, 0.5, 0.5),
		Point<FloatingPointType>(0.5001, 0.5, 0.5)
	};

	std::vector<geometry::Symm<FloatingPointType>> symms = {
		geometry::Symm<FloatingPointType>("x, y, z")
	};

	cpplib::cluster_detail::UnitCellBuilder ucb(symms);
	auto buildresult = ucb.build(closePoints,
								   std::vector<AtomTypeBase>(closePoints.size(), 1));

	geometry::SpatialGrid<FloatingPointType> space;
	space.build(buildresult.atoms.points, cell, 6.0);
	auto bonds = space.get_bonds();

	// Should not crash
	EXPECT_NO_THROW({
		voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds,
									 cell);
		auto cells = vd.extractCells();
	});
}

TEST(VoronoiTest, EmptyDiagram) {
	using namespace cpplib::voronoi;

	voronoi::VoronoiDiagram::PointVector emptyPoints;
	geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

	std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift> emptyBonds;

	EXPECT_NO_THROW({
		voronoi::VoronoiDiagram vd(emptyPoints, emptyBonds, cell);
		auto cells = vd.extractCells();
		EXPECT_EQ(cells.size(), 0);
	});
}

TEST(VoronoiTest, AlexTest) {
    cpplib::geometry::Cell cell (10.9815, 6.8214, 8.9974, 90.0, 101.511, 90.0);
    std::vector<const char*> symms{"x, y, z", "-x, y+1/2, -z+1/2", "-x, -y, -z", "x, -y-1/2, z-1/2"};

    voronoi::VoronoiDiagram::PointVector data = {{0.85991, 0.6099,  0.54993},
                                                 {0.9802,  0.2973,  0.68819},
                                                 {0.45147, 0.2722, -0.01408},
                                                 {0.4257,  0.369,  -0.0692 },
                                                 {0.54605, 0.33304, 0.10639},
                                                 {0.59711, 0.18811, 0.18507},
                                                 {0.5699,  0.0585,  0.1577 },
                                                 {0.6972,  0.22032, 0.31696},
                                                 {0.73155, 0.40837, 0.37126},
                                                 {0.6897,  0.5202,  0.323  },
                                                 {0.8262,  0.42976, 0.49478},
                                                 {0.88775, 0.26835, 0.56724},
                                                 {0.85479, 0.08196, 0.51712},
                                                 {0.8965, -0.0288,  0.5674 },
                                                 {0.75907, 0.05834, 0.39111},
                                                 {0.7356, -0.0699,  0.3552 }};

	voronoi::VoronoiDiagram::BoolVector bools{true, true, true, true,
                                              true, true,  true, true, 
                                              true, true, true, true, 
                                              true, true, true, true};

	std::vector<geometry::Symm<FloatingPointType>> symmvec;
	symmvec.reserve(symms.size());
	for (int i = 0; i < symms.size(); i++)
	{
		symmvec.emplace_back(symms[i]);
	}

	cpplib::cluster_detail::UnitCellBuilder ucb(symmvec);
	auto buildresult = ucb.build(data, std::vector<AtomTypeBase>(data.size(), AtomTypeBase(1)));

	cpplib::geometry::SpatialGrid<FloatingPointType> space;
	space.build(buildresult.atoms.points, cell, 6);
	auto bonds = space.get_bonds();

	voronoi::VoronoiDiagram vd(buildresult.atoms.points,bonds,cell, bools);
    auto cells = vd.extractCells();

	std::array< std::array<uint32_t, 3>, 16> dead = {{{93, 97, 25},
                                                      {125, 129, 30},
                                                      {65, 69, 20},
                                                      {96, 100, 24},
                                                      {83, 87, 21},
                                                      {79, 83, 24},
                                                      {90, 94, 25},
                                                      {87, 91, 24},
                                                      {76, 80, 23},
                                                      {94, 98, 24},
                                                      {90, 94, 24},
                                                      {71, 75, 20},
                                                      {71, 75, 22},
                                                      {80, 84, 22},
                                                      {82, 86, 23},
                                                      {91, 95, 24}}};

	std::array< std::array<uint32_t, 3>, 16> alive = {{{30, 45, 17},
													   {42, 63, 23},
													   {22, 33, 13},
													   {32, 48, 18},
													   {26, 39, 15},
													   {26, 39, 15},
													   {32, 48, 18},
													   {30, 45, 17},
													   {24, 36, 14},
													   {32, 48, 18},
													   {30, 45, 17},
													   {22, 33, 13},
													   {22, 33, 13},
													   {26, 39, 15},
													   {28, 42, 16},
													   {32, 48, 18}}};

	for (size_t i = 0; i < 16; i++)
	{
		EXPECT_EQ(cells[i].vertices.size(), dead[i][0]);
		EXPECT_EQ(cells[i].edges.size(), dead[i][1]);
		EXPECT_EQ(cells[i].faces.size(), dead[i][2]);


		int a=0, b=0, c=0;
		for (auto& v : cells[i].vertices)
		{
			if (v->get_state() == voronoi::State::VALID) a++;
		}
		for (auto& v : cells[i].edges)
		{
			if (v->get_state() == voronoi::State::VALID) b++;

		}
		for (auto& v : cells[i].faces)
		{
			if (v->get_state() == voronoi::State::VALID) c++;
		}

		EXPECT_EQ(a, alive[i][0]);
		EXPECT_EQ(b, alive[i][1]);
		EXPECT_EQ(c, alive[i][2]);
	}
	voronoi::VoronoiFused vf(cells, cell.fracToCart());
	// Verify polygon topology
	for (const auto& polygon : vf.polygons) {
		// Each polygon should have equal number of vertices and edges
		EXPECT_EQ(polygon.vert_ids.size(), polygon.edge_ids.size())
			<< "Polygon should have equal vertices and edges";

		// Each polygon should have at least 3 vertices
		EXPECT_GE(polygon.vert_ids.size(), 3)
			<< "Polygon should have at least 3 vertices";

		// Verify all vertex and edge IDs are valid
		for (auto vid : polygon.vert_ids) {
			EXPECT_LT(vid, vf.vertices.size());
		}
		for (auto eid : polygon.edge_ids) {
			EXPECT_LT(eid, vf.edges.size());
		}
	}

	// Verify edge topology
	for (const auto& edge : vf.edges) {
		// Each edge should connect two different vertices
		EXPECT_NE(edge.vert_ids[0], edge.vert_ids[1]);
		EXPECT_LT(edge.vert_ids[0], vf.vertices.size());
		EXPECT_LT(edge.vert_ids[1], vf.vertices.size());
	}
	// Verify sizes
	EXPECT_EQ(vf.vertices.size(), 305);
	EXPECT_EQ(vf.edges.size(), 515);
	EXPECT_EQ(vf.polygons.size(), 228);
	for (size_t i = 0; i < 16; i++)
	{
		EXPECT_EQ(vf.polyhedra[i].vert_ids.size(), alive[i][0]);
		EXPECT_EQ(vf.polyhedra[i].edge_ids.size(), alive[i][1]);
		EXPECT_EQ(vf.polyhedra[i].poly_ids.size(), alive[i][2]);
	}

}

TEST(VoronoiTest, Benzene) {
	using namespace cpplib::voronoi;

	geometry::Cell<FloatingPointType> cell(7.243, 9.310, 6.756, 90.0, 90.0, 90.0);
	std::vector<const char*> symms{
		"x, y, z",
		"-x+1/2, -y, z+1/2",
		"-x, y+1/2, -z+1/2",
		"x+1/2, -y+1/2, -z",
		"-x, -y, -z",
		"x+1/2, y, -z+1/2",
		"x, -y+1/2, z+1/2",
		"-x+1/2, y+1/2, z"
	};

	VoronoiDiagram::PointVector data = {{
		{-0.06070, 0.13930, -0.00690},
		{-0.13770, 0.04470,  0.12600},
		{ 0.07700, 0.09580, -0.13250},
		{-0.10460, 0.25020, -0.01230},
		{-0.24580, 0.07810,  0.22410},
		{ 0.13710, 0.16810, -0.23600}
	}};

	voronoi::VoronoiDiagram::BoolVector bools{false, true, false, false,
											  false, true};

	std::vector<geometry::Symm<FloatingPointType>> symmvec;
	symmvec.reserve(symms.size());
	for (int i = 0; i < symms.size(); i++)
	{
		symmvec.emplace_back(symms[i]);
	}

	cpplib::cluster_detail::UnitCellBuilder ucb(symmvec);
	auto buildresult = ucb.build(data, std::vector<AtomTypeBase>(data.size(), AtomTypeBase(1)));

	cpplib::geometry::SpatialGrid<FloatingPointType> space;
	space.build(buildresult.atoms.points, cell, 6);
	auto bonds = space.get_bonds();

	voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds, cell, bools);
	auto cells = vd.extractCells();
}


// ==================== ENTRY POINT ====================
int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}