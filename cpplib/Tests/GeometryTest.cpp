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
            } else {
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
            } else {
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
//TEST(VoronoiTest, VertexEquality) {
//    using namespace cpplib::voronoi;
//
//    // Test vertex equality with epsilon
//    Vertex v1(0, Point<FloatingPointType>(1.0, 2.0, 3.0));
//    Vertex v2(1, Point<FloatingPointType>(1.0 + voronoi::EPSILON / 2,
//                                          2.0 + voronoi::EPSILON / 2,
//                                          3.0 + voronoi::EPSILON / 2));
//    Vertex v3(2, Point<FloatingPointType>(1.0 + voronoi::EPSILON * 2,
//                                          2.0, 3.0));
//
//    EXPECT_TRUE(v1 == v2);  // Within epsilon
//    EXPECT_FALSE(v1 == v3); // Outside epsilon
//}
//
//TEST(VoronoiTest, VertexDistanceCalculation) {
//    using namespace cpplib::voronoi;
//
//    voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);
//    geometry::Matrix<FloatingPointType> fractocart(1.0); // Identity
//
//    const auto* maxVert = cell.update_vertices_distances(fractocart);
//
//    EXPECT_NE(maxVert, nullptr);
//    EXPECT_GT(maxVert->distance, 0.0);
//
//    // Verify all vertices have computed distances
//    for (const auto& vert : cell.vertices) {
//        if (vert->get_state() != voronoi::State::DELETE) {
//            EXPECT_GT(vert->distance, 0.0);
//        }
//    }
//}
//
//TEST(VoronoiTest, EdgeStateCalculation) {
//    using namespace cpplib::voronoi;
//
//    // Create vertices
//    auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
//    auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));
//
//    Edge edge(0, v1.get(), v2.get());
//
//    // Initially, edge should be in INVALID state (no faces)
//    EXPECT_EQ(edge.calculateState(), State::INVALID);
//
//    // Create faces
//    auto f1 = std::make_unique<Face>(0);
//    auto f2 = std::make_unique<Face>(1);
//
//    // Add edge to faces
//    f1->edges.insert(&edge);
//    f2->edges.insert(&edge);
//    edge.faces.insert(f1.get());
//    edge.faces.insert(f2.get());
//
//    // Now edge should be VALID (2 vertices, 2 faces, both vertices valid)
//    EXPECT_EQ(edge.calculateState(), State::VALID);
//
//    // Delete one vertex
//    v1->set_state(State::DELETE);
//    EXPECT_EQ(edge.calculateState(), State::MODIFICATION);
//
//    // Delete both vertices
//    v2->set_state(State::DELETE);
//    EXPECT_EQ(edge.calculateState(), State::DELETE);
//}
//
//TEST(VoronoiTest, EdgeIntersection) {
//    using namespace cpplib::voronoi;
//
//    auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
//    auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));
//    v1->set_state(State::VALID);
//    v2->set_state(State::DELETE);
//
//    Edge edge(0, v1.get(), v2.get());
//    auto f1 = std::make_unique<Face>(0);
//    auto f2 = std::make_unique<Face>(1);
//    f1->edges.insert(&edge);
//    f2->edges.insert(&edge);
//    edge.faces.insert(f1.get());
//    edge.faces.insert(f2.get());
//
//    // Edge in MODIFICATION state
//    EXPECT_EQ(edge.calculateState(), State::MODIFICATION);
//
//    // Create a plane cutting the edge at x=0.5
//    geometry::Plane<FloatingPointType> plane(
//        Point<FloatingPointType>(0.5, 0.0, 0.0),
//        Point<FloatingPointType>(1.0, 0.0, 0.0)
//    );
//
//    auto intersection = edge.intersectSegmentPlane(plane);
//
//    EXPECT_NEAR(intersection[0], 0.5, voronoi::EPSILON);
//    EXPECT_NEAR(intersection[1], 0.0, voronoi::EPSILON);
//    EXPECT_NEAR(intersection[2], 0.0, voronoi::EPSILON);
//}
//
//TEST(VoronoiTest, FaceStateCalculation) {
//    using namespace cpplib::voronoi;
//
//    Face face(0);
//
//    // Create vertices
//    auto v1 = std::make_unique<Vertex>(0, Point<FloatingPointType>(0.0, 0.0, 0.0));
//    auto v2 = std::make_unique<Vertex>(1, Point<FloatingPointType>(1.0, 0.0, 0.0));
//    auto v3 = std::make_unique<Vertex>(2, Point<FloatingPointType>(1.0, 1.0, 0.0));
//
//    // Create edges
//    Edge e1(0, v1.get(), v2.get());
//    Edge e2(1, v2.get(), v3.get());
//    Edge e3(2, v3.get(), v1.get());
//
//    face.edges.insert({&e1, &e2, &e3});
//    face.vertices.insert({v1.get(), v2.get(), v3.get()});
//
//    // All edges VALID -> face VALID
//    e1.set_state(State::VALID);
//    e2.set_state(State::VALID);
//    e3.set_state(State::VALID);
//    EXPECT_EQ(face.calculateState(), State::VALID);
//
//    // One vertex DELETE -> face MODIFICATION
//    v1->set_state(State::DELETE);
//    EXPECT_EQ(face.calculateState(), State::MODIFICATION);
//
//    // All vertexes DELETE -> face DELETE
//    v1->set_state(State::DELETE);
//    v2->set_state(State::DELETE);
//    v3->set_state(State::DELETE);
//    EXPECT_EQ(face.calculateState(), State::DELETE);
//}
//
//TEST(VoronoiTest, CellClippingParallelPlane) {
//    using namespace cpplib::voronoi;
//
//    voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);
//
//    size_t initialVertexCount = cell.vertices.size();
//    size_t initialEdgeCount = cell.edges.size();
//    size_t initialFaceCount = cell.faces.size();
//
//    // Plane parallel to a face of the cube (shouldn't clip anything)
//    geometry::Plane<FloatingPointType> parallelPlane(
//        Point<FloatingPointType>(2.0, 0.5, 0.5),
//        Point<FloatingPointType>(-1.0, 0.0, 0.0)
//    );
//
//    cell.clipByPlaneAndAddNewFace(parallelPlane, 1, ShiftCode(14));
//
//    // Should not change the cell
//    EXPECT_EQ(cell.vertices.size(), initialVertexCount);
//    EXPECT_EQ(cell.edges.size(), initialEdgeCount);
//    EXPECT_EQ(cell.faces.size(), initialFaceCount);
//}
//
//TEST(VoronoiTest, CellClippingTangent) {
//    using namespace cpplib::voronoi;
//
//    voronoi::Cell cell(Point<FloatingPointType>(0.5, 0.5, 0.5), 0);
//
//    // Plane tangent to the cube (just touches an edge)
//    geometry::Plane<FloatingPointType> tangentPlane(
//        Point<FloatingPointType>(1.0, 0.5, 0.5),
//        Point<FloatingPointType>(-1.0, 0.0, 0.0)
//    );
//
//    size_t initialVertexCount = cell.vertices.size();
//
//    cell.clipByPlaneAndAddNewFace(tangentPlane, 1, ShiftCode(14));
//
//    // Should add some vertices at the tangent points
//    EXPECT_GE(cell.vertices.size(), initialVertexCount);
//}
//
//TEST(VoronoiTest, VoronoiFusedVertexMerging) {
//    using namespace cpplib::voronoi;
//
//    // Create two cells with overlapping vertices
//    std::vector<voronoi::Cell > cells;
//    cells.emplace_back(Point<FloatingPointType>(0.3, 0.5, 0.5), 0);
//    cells.emplace_back(Point<FloatingPointType>(0.7, 0.5, 0.5), 1);
//
//    // Clip both cells so they share a face
//    geometry::Plane<FloatingPointType> plane1(
//        Point<FloatingPointType>(0.5, 0.5, 0.5),
//        Point<FloatingPointType>(-1.0, 0.0, 0.0)
//    );
//    cells[0].clipByPlaneAndAddNewFace(plane1, 1, ShiftCode());
//
//    geometry::Plane<FloatingPointType> plane2(
//        Point<FloatingPointType>(0.5, 0.5, 0.5),
//        Point<FloatingPointType>(1.0, 0.0, 0.0)
//    );
//    cells[1].clipByPlaneAndAddNewFace(plane2, 0, ShiftCode());
//
//    geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
//    VoronoiFused fused(cells, cell.fracToCart());
//
//    // Verify that duplicate vertices were merged
//    for (size_t i = 0; i < fused.vertices.size(); ++i) {
//        for (size_t j = i + 1; j < fused.vertices.size(); ++j) {
//            FloatingPointType dist = (fused.vertices[i] - fused.vertices[j]).r();
//            EXPECT_GT(dist, voronoi::EPSILON);
//        }
//    }
//
//    // Verify polyhedra have correct vertex references
//    EXPECT_GT(fused.polyhedra[0].vert_ids.size(), 0);
//    EXPECT_GT(fused.polyhedra[1].vert_ids.size(), 0);
//}
//
//TEST(VoronoiTest, VoronoiFusedEdgePolygonPopulation) {
//    using namespace cpplib::voronoi;
//
//    // Use the simple two-cell setup
//    std::vector<voronoi::Cell> cells;
//    cells.emplace_back(Point<FloatingPointType>(0.3, 0.5, 0.5), 0);
//    cells.emplace_back(Point<FloatingPointType>(0.7, 0.5, 0.5), 1);
//
//    geometry::Plane<FloatingPointType> plane1(
//        Point<FloatingPointType>(0.5, 0.5, 0.5),
//        Point<FloatingPointType>(-1.0, 0.0, 0.0)
//    );
//    cells[0].clipByPlaneAndAddNewFace(plane1, 1, ShiftCode());
//
//    geometry::Plane<FloatingPointType> plane2(
//        Point<FloatingPointType>(0.5, 0.5, 0.5),
//        Point<FloatingPointType>(1.0, 0.0, 0.0)
//    );
//    cells[1].clipByPlaneAndAddNewFace(plane2, 0, ShiftCode());
//
//    geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
//    VoronoiFused fused(cells, cell.fracToCart());
//
//    // Verify that polyhedra have edge_ids and poly_ids populated
//    for (size_t i = 0; i < fused.polyhedra.size(); ++i) {
//        EXPECT_GT(fused.polyhedra[i].edge_ids.size(), 0)
//            << "Polyhedron " << i << " should have edges";
//        EXPECT_GT(fused.polyhedra[i].poly_ids.size(), 0)
//            << "Polyhedron " << i << " should have polygons";
//
//        // Verify all edge IDs are valid
//        for (auto edge_id : fused.polyhedra[i].edge_ids) {
//            EXPECT_LT(edge_id, fused.edges.size());
//        }
//
//        // Verify all polygon IDs are valid
//        for (auto poly_id : fused.polyhedra[i].poly_ids) {
//            EXPECT_LT(poly_id, fused.polygons.size());
//        }
//    }
//}
//
//TEST(VoronoiTest, DegenerateCases) {
//    using namespace cpplib::voronoi;
//
//    // Test with points very close together
//    geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
//
//    voronoi::VoronoiDiagram::PointVector closePoints = {
//        Point<FloatingPointType>(0.5, 0.5, 0.5),
//        Point<FloatingPointType>(0.5001, 0.5, 0.5)
//    };
//
//    std::vector<geometry::Symm<FloatingPointType>> symms = {
//        geometry::Symm<FloatingPointType>("x, y, z")
//    };
//
//    cpplib::cluster_detail::UnitCellBuilder ucb(symms);
//    auto buildresult = ucb.build(closePoints,
//                                 std::vector<AtomTypeBase>(closePoints.size(), 1));
//
//    geometry::SpatialGrid<FloatingPointType> space;
//    space.build(buildresult.atoms.points, cell, 6.0);
//    auto bonds = space.get_bonds();
//
//    // Should not crash
//    EXPECT_NO_THROW({
//        voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds,
//                                     cell);
//        auto cells = vd.extractCells();
//                    });
//}
//
//TEST(VoronoiTest, EmptyDiagram) {
//    using namespace cpplib::voronoi;
//
//    voronoi::VoronoiDiagram::PointVector emptyPoints;
//    geometry::Cell<FloatingPointType> cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
//
//    std::vector<geometry::SpatialGrid<FloatingPointType>::BondWithShift> emptyBonds;
//
//    EXPECT_NO_THROW({
//        voronoi::VoronoiDiagram vd(emptyPoints, emptyBonds, cell);
//        auto cells = vd.extractCells();
//        EXPECT_EQ(cells.size(), 0);
//                    });
//}
//
//TEST(VoronoiTest, AlexTest) {
//    cpplib::geometry::Cell cell(10.9815, 6.8214, 8.9974, 90.0, 101.511, 90.0);
//    std::vector<const char*> symms{"x, y, z", "-x, y+1/2, -z+1/2", "-x, -y, -z", "x, -y-1/2, z-1/2"};
//
//    voronoi::VoronoiDiagram::PointVector data = {{0.85991, 0.6099,  0.54993},
//                                                 {0.9802,  0.2973,  0.68819},
//                                                 {0.45147, 0.2722, -0.01408},
//                                                 {0.4257,  0.369,  -0.0692 },
//                                                 {0.54605, 0.33304, 0.10639},
//                                                 {0.59711, 0.18811, 0.18507},
//                                                 {0.5699,  0.0585,  0.1577 },
//                                                 {0.6972,  0.22032, 0.31696},
//                                                 {0.73155, 0.40837, 0.37126},
//                                                 {0.6897,  0.5202,  0.323  },
//                                                 {0.8262,  0.42976, 0.49478},
//                                                 {0.88775, 0.26835, 0.56724},
//                                                 {0.85479, 0.08196, 0.51712},
//                                                 {0.8965, -0.0288,  0.5674 },
//                                                 {0.75907, 0.05834, 0.39111},
//                                                 {0.7356, -0.0699,  0.3552 }};
//
//    voronoi::VoronoiDiagram::BoolVector bools{true, true, true, true,
//                                              true, true,  true, true,
//                                              true, true, true, true,
//                                              true, true, true, true};
//
//    std::vector<geometry::Symm<FloatingPointType>> symmvec;
//    symmvec.reserve(symms.size());
//    for (int i = 0; i < symms.size(); i++)
//    {
//        symmvec.emplace_back(symms[i]);
//    }
//
//    cpplib::cluster_detail::UnitCellBuilder ucb(symmvec);
//    auto buildresult = ucb.build(data, std::vector<AtomTypeBase>(data.size(), AtomTypeBase(1)));
//
//    cpplib::geometry::SpatialGrid<FloatingPointType> space;
//    space.build(buildresult.atoms.points, cell, 6);
//    auto bonds = space.get_bonds();
//
//    voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds, cell, bools);
//    auto cells = vd.extractCells();
//
//    std::array< std::array<uint32_t, 3>, 16> dead = {{{93, 97, 25},
//                                                      {125, 129, 30},
//                                                      {65, 69, 20},
//                                                      {96, 100, 24},
//                                                      {83, 87, 21},
//                                                      {79, 83, 24},
//                                                      {90, 94, 25},
//                                                      {87, 91, 24},
//                                                      {76, 80, 23},
//                                                      {94, 98, 24},
//                                                      {90, 94, 24},
//                                                      {71, 75, 20},
//                                                      {71, 75, 22},
//                                                      {80, 84, 22},
//                                                      {82, 86, 23},
//                                                      {91, 95, 24}}};
//
//    std::array< std::array<uint32_t, 3>, 16> alive = {{{30, 45, 17},
//                                                       {42, 63, 23},
//                                                       {22, 33, 13},
//                                                       {32, 48, 18},
//                                                       {26, 39, 15},
//                                                       {26, 39, 15},
//                                                       {32, 48, 18},
//                                                       {30, 45, 17},
//                                                       {24, 36, 14},
//                                                       {32, 48, 18},
//                                                       {30, 45, 17},
//                                                       {22, 33, 13},
//                                                       {22, 33, 13},
//                                                       {26, 39, 15},
//                                                       {28, 42, 16},
//                                                       {32, 48, 18}}};
//
//    for (size_t i = 0; i < 16; i++)
//    {
//        EXPECT_EQ(cells[i].vertices.size(), dead[i][0]);
//        EXPECT_EQ(cells[i].edges.size(), dead[i][1]);
//        EXPECT_EQ(cells[i].faces.size(), dead[i][2]);
//
//
//        int a = 0, b = 0, c = 0;
//        for (auto& v : cells[i].vertices)
//        {
//            if (v->get_state() == voronoi::State::VALID) a++;
//        }
//        for (auto& v : cells[i].edges)
//        {
//            if (v->get_state() == voronoi::State::VALID) b++;
//
//        }
//        for (auto& v : cells[i].faces)
//        {
//            if (v->get_state() == voronoi::State::VALID) c++;
//        }
//
//        EXPECT_EQ(a, alive[i][0]);
//        EXPECT_EQ(b, alive[i][1]);
//        EXPECT_EQ(c, alive[i][2]);
//    }
//    voronoi::VoronoiFused vf(cells, cell.fracToCart());
//    // Verify polygon topology
//    for (const auto& polygon : vf.polygons) {
//        // Each polygon should have equal number of vertices and edges
//        EXPECT_EQ(polygon.vert_ids.size(), polygon.edge_ids.size())
//            << "Polygon should have equal vertices and edges";
//
//        // Each polygon should have at least 3 vertices
//        EXPECT_GE(polygon.vert_ids.size(), 3)
//            << "Polygon should have at least 3 vertices";
//
//        // Verify all vertex and edge IDs are valid
//        for (auto vid : polygon.vert_ids) {
//            EXPECT_LT(vid, vf.vertices.size());
//        }
//        for (auto eid : polygon.edge_ids) {
//            EXPECT_LT(eid, vf.edges.size());
//        }
//    }
//
//    // Verify edge topology
//    for (const auto& edge : vf.edges) {
//        // Each edge should connect two different vertices
//        EXPECT_NE(edge.vert_ids[0], edge.vert_ids[1]);
//        EXPECT_LT(edge.vert_ids[0], vf.vertices.size());
//        EXPECT_LT(edge.vert_ids[1], vf.vertices.size());
//    }
//    // Verify sizes
//    EXPECT_EQ(vf.vertices.size(), 305);
//    EXPECT_EQ(vf.edges.size(), 515);
//    EXPECT_EQ(vf.polygons.size(), 228);
//    for (size_t i = 0; i < 16; i++)
//    {
//        EXPECT_EQ(vf.polyhedra[i].vert_ids.size(), alive[i][0]);
//        EXPECT_EQ(vf.polyhedra[i].edge_ids.size(), alive[i][1]);
//        EXPECT_EQ(vf.polyhedra[i].poly_ids.size(), alive[i][2]);
//    }
//
//}
//
//TEST(VoronoiTest, Benzene) {
//    using namespace cpplib::voronoi;
//
//    geometry::Cell<FloatingPointType> cell(7.243, 9.310, 6.756, 90.0, 90.0, 90.0);
//    std::vector<const char*> symms{
//        "x, y, z",
//        "-x+1/2, -y, z+1/2",
//        "-x, y+1/2, -z+1/2",
//        "x+1/2, -y+1/2, -z",
//        "-x, -y, -z",
//        "x+1/2, y, -z+1/2",
//        "x, -y+1/2, z+1/2",
//        "-x+1/2, y+1/2, z"
//    };
//
//    VoronoiDiagram::PointVector data = {{
//        {-0.06070, 0.13930, -0.00690},
//        {-0.13770, 0.04470,  0.12600},
//        { 0.07700, 0.09580, -0.13250},
//        {-0.10460, 0.25020, -0.01230},
//        {-0.24580, 0.07810,  0.22410},
//        { 0.13710, 0.16810, -0.23600}
//    }};
//
//    voronoi::VoronoiDiagram::BoolVector bools{false, true, false, false,
//                                              false, true};
//
//    std::vector<geometry::Symm<FloatingPointType>> symmvec;
//    symmvec.reserve(symms.size());
//    for (int i = 0; i < symms.size(); i++)
//    {
//        symmvec.emplace_back(symms[i]);
//    }
//
//    cpplib::cluster_detail::UnitCellBuilder ucb(symmvec);
//    auto buildresult = ucb.build(data, std::vector<AtomTypeBase>(data.size(), AtomTypeBase(1)));
//
//    cpplib::geometry::SpatialGrid<FloatingPointType> space;
//    space.build(buildresult.atoms.points, cell, 6);
//    auto bonds = space.get_bonds();
//
//    voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds, cell, bools);
//    auto cells = vd.extractCells();
//}
//
//TEST(VoronoiTest, AlexTestMass) {
//    using namespace cpplib::voronoi;
//
//    geometry::Cell<FloatingPointType> cell(48.616676330566406, 30.1771297454834, 47.705909729003906, 90.0, 90.0, 90.0);
//    std::vector<const char*> symms{
//            "x, y, z"
//    };
//    VoronoiDiagram::PointVector data = {{
//            {0.0050663175, 0.61703503, 0.5423728},
//            {0.19670719, 0.6375811, 0.2896085},
//            {0.054837663, 0.4195399, 0.032163966},
//            {0.1944302, 0.66899794, 0.29418913},
//            {0.021821218, 0.74489874, 0.2662075},
//            {0.047880176, 0.3729707, 0.34968776},
//            {0.11467207, 0.6363995, -0.18090159},
//            {0.11941582, 0.5623753, 0.31390575},
//            {0.090023585, 0.42582327, 0.4648009},
//            {0.48942775, 0.47910678, 0.28635558},
//            {0.090023585, 0.61677086, 0.4648009},
//            {0.21232992, 0.63403624, 0.27676284},
//            {0.2661556, 0.53026336, 0.50974417},
//            {-0.11467207, 0.7537261, 0.5128311},
//            {-0.15932651, 0.63591295, 0.31118393},
//            {0.13820103, 0.63521785, 0.28864592},
//            {0.050220422, 0.9005233, 0.09074953},
//            {-0.21024267, 0.75379556, 0.45866022},
//            {0.20024918, 0.6208301, 0.30683565},
//            {0.13099055, 0.6678858, 0.34948862},
//            {-0.30315676, 0.6166597, 0.2260108},
//            {0.23598538, 0.8597927, 0.20151441},
//            {-0.034471195, 0.8688285, 0.5638819},
//            {0.00594549, 0.8767522, 0.5589029},
//            {-0.24737035, 0.8688285, 0.15192415},
//            {0.24338561, 0.40459603, 0.113752246},
//            {0.22535941, 0.48960224, 0.29803953},
//            {0.23480892, 0.54221845, 0.59823656},
//            {0.16510755, 0.6165763, 0.2727465},
//            {-0.022833215, 0.5212971, 0.1624463},
//            {-0.022833215, 0.37693256, 0.28077918},
//            {0.13402654, 0.63591295, 0.34772936},
//            {0.05818991, 0.56237525, 0.6467315},
//            {0.40796188, 0.5943482, -0.06505818},
//            {-0.25280985, 0.4061252, 0.35008606},
//            {0.05836701, 0.42587888, 0.3150343},
//            {-0.14059821, 0.6164998, -0.04530838},
//            {0.4731725, 0.40668124, 0.31118393},
//            {0.119770005, 0.6164582, -0.01729353},
//            {0.09089011, 0.8688285, 0.46559754},
//            {-0.2661556, 0.87779474, 0.15411487},
//            {0.15382376, 0.41286728, 0.3482273},
//            {0.03870894,  0.636469, 0.2910026},
//            {-0.17317823,  0.215956, 0.12039083},
//            {0.20739642, 0.5113577, -0.026853103},
//            {0.15382376, 0.6297269, 0.3482273},
//            {0.33908272, 0.72446394, 0.11481442},
//            {-0.35467383, 0.54395616, 0.35954607},
//            {0.124981806, 0.62048256, 0.36329687},
//            {0.09089011, 0.9005233, 0.46400428},
//            {0.20739642, 0.5113577, 0.30507642},
//            {0.2229559, 0.6211081, 0.5498081},
//            {0.22453715, 0.79577744, 0.3969877},
//            {0.038424317, 0.8461695, 0.52551085},
//            {-0.06343965, 0.75365657, 0.5160509},
//            {0.20739642, 0.5312365, 0.30507642},
//            {0.15382376, 0.62972695, 0.68015677},
//            {0.03870894,  0.636469, -0.040926915},
//            {-0.004427493, 1.0122894, 0.2025766},
//            {-0.015432976, 0.21498291, 0.06034479},
//            {0.11977002, 0.6164582,  0.314636},
//            {-0.01922797, 0.41043457, 0.6384664},
//            {0.16407025, 0.56209725, -0.059714116},
//            {0.16463949,  0.636052, 0.24134597},
//            {-0.021947715, 0.83713365, 0.27892038},
//            {-0.050093923, 0.5123308, -0.011849892},
//            {-0.17804848, 0.4073763, 0.20924836},
//            {0.037064448, 0.87779474, 0.106217444},
//            {0.0066412394, 0.41015655, 0.46128246},
//            {-0.18285546, 0.83713365, 0.23626743},
//            {0.1646395, 0.63605195, -0.090583555},
//            {0.16407025, 0.56209725, 0.2722154},
//            {0.069258645, 0.73808724, 0.17734995},
//            {0.1686875, 0.40709826, 0.46184677},
//            {0.1547219, 0.54451215, 0.4122233},
//            {0.02536321, 0.41168568, 0.34998646},
//            {0.27918506, 0.53026336, 0.27218223},
//            {-0.037254192, 0.6270856, -0.14438935},
//            {0.015432976, 0.48007986, -0.06034479},
//            {0.091712356, 0.5943482, 0.23102295},
//            {-0.047880176, 1.0171549, -0.017758248},
//            {0.076595634, 0.83713365, -0.037574425},
//            {0.2712156,  0.985043, 0.5145903},
//            {0.11941582, 0.4802189, 0.31390575},
//            {0.40574813, 0.40702876, 0.18767294},
//            {-0.031688202, 0.6674688, 0.20848495},
//            {0.045033928, 0.6375116, 0.34862557},
//            {0.034344696, 0.7740219, 0.05340746},
//            {0.034471195, 0.5212971, 0.09997717},
//            {0.02536321, 0.6309085, 0.3499865},
//            {-0.012523481, 0.76755786, 0.22119783},
//            {0.3327577, 0.40668124, 0.25780967},
//            {0.047880176, 0.6696235, 0.34968776},
//            {-0.01309273, 0.7734659, 0.060046043},
//            {-0.105247825, 0.66830283, 0.2892102},
//            {-0.23079889, 0.66809434, 0.40209946},
//            {-0.030676201, 0.63521785, 0.20622781},
//            {0.19252637, 0.5212971, 0.1549779},
//            {0.054837663, 0.6230543, 0.36409348},
//            {0.35350367, 0.76303995, 0.35350496},
//            {-0.054837663, 0.97058564, 0.29976556},
//            {0.16315944, 0.5212971, -0.15033753},
//            {0.05836701, 0.6167153, 0.3150343},
//            {-0.00594549, 0.5133734, 0.43688565},
//            {0.111762576, 0.4053606, 0.42247987},
//            {0.2253594, 0.5529919, 0.29803953},
//            {-0.026628207, 0.5529919, 0.4398398},
//            {0.013092729, 0.6166597, 0.2718835},
//            {-0.27754056, 0.4061252, -0.12503785},
//            {0.08949861, 0.63556546, 0.31022134},
//            {0.01309273, 0.4259345, 0.2718835},
//            {-0.07684863, 0.8688285, 0.27739352},
//            {0.18285546, 0.48960224, 0.42759162},
//            {0.08545062, 0.37449983, 0.4277244},
//            {0.3642752, 0.5408284, 0.4243387},
//            {0.2229559, 0.42148608, -0.11405098},
//            {-0.0067677395, 0.5212971, 0.65532845},
//            {0.119770005, 0.42613602,  0.314636},
//            {-0.25280985, 0.4061252, 0.01815654},
//            {0.05836701, 0.42587888, 0.64696383},
//            {0.090320855, 0.56279236, -0.19806236},
//            {0.05818991, 0.5623753, 0.31480196},
//            {-0.022833215, 0.66566163, 0.28077918},
//            {-0.03870894, 0.9840004, 0.040926915},
//            {0.13339405, 0.9005233, 0.26162684},
//            {0.16510755, 0.42601788, 0.2727465},
//            {0.13099055, 0.37470832, 0.017559057},
//            {0.17279872, 0.5212971, 0.57419497},
//            {-0.010562734, 0.5461803, 0.0011617483},
//            {0.016508223, 0.63591295, 0.24008463},
//            {0.5272512, 0.3742913, -0.042719338},
//            {-0.034344696, 0.42649055, 0.27852207},
//            {-0.23193738, 0.87876785, 0.22156295},
//            {0.19974321, 0.53019387, 0.41643876},
//            {0.21549243, 0.41231126, 0.18571456},
//            {-0.35350367, 0.41550857, -0.021575417},
//            {0.015432976, 0.5625143, 0.27158472},
//            {-0.045033928, 0.75261396, -0.016696054},
//            {-0.18728295, 0.87876785, 0.5490446},
//            {0.27918506, 0.5123308, 0.27218223},
//            {-0.1944302, 1.0165293, 0.3696699},
//            {0.034344696, 0.9636351, 0.38533697},
//            {0.08880919, 0.5212971, 0.6455332},
//            {0.27539006, 0.5312365, 0.23955354},
//            {-0.170838, 0.32021543, 0.20051864},
//            {0.0316882, 1.0150002, 0.45537412},
//            {-0.2253594, 0.8688285, 0.36422628},
//            {0.2753901, 0.5113577, 0.23955354},
//            {0.13086404, 0.5212971, 0.09556251},
//            {-0.30967152, 0.56320935, 0.28642198},
//            {0.033022773, 0.5212971, 0.6290828},
//            {0.004427493, 0.37783617, 0.46128246},
//            {-0.022833215, 0.37693256, 0.6127087},
//            {0.16510755, 0.6165763,  0.604676},
//            {-0.022833215, 0.5212971, 0.49437582},
//            {0.28905204, 0.45109576, 0.14146836},
//            {0.23480894, 0.5422185, 0.26630706},
//            {0.12959905, 0.76887846, 0.29823866},
//            {0.21644115, 0.5212971, 0.24360308},
//            {0.35467383, 0.8461695, -0.027616538},
//            {0.1944302, 0.37359628, -0.03774038},
//            {0.4412313, 0.62048256, 0.13459744},
//            {0.26754707, 0.5212971, 0.25691345},
//            {0.05818991, 0.4802189,  0.314802},
//            {0.07615288,  0.575651, 0.2039375},
//            {0.08026413, 0.5122613, 0.36747918},
//            {0.13402654, 0.63591295, 0.67965895},
//            {-0.35141644, 0.5321401, 0.07478373},
//            {0.5164987, 0.62083006, 0.19105865},
//            {0.22535941, 0.5212971, 0.2996328},
//            {0.23851538, 0.53026336, 0.31350744},
//            {-0.026628207, 0.48960224, 0.4398398},
//            {-0.05818991, 0.9099067, 0.017127557},
//            {-0.14059821, 0.42609435, -0.04530838},
//            {0.29430178, 0.9005233, 0.2189739},
//            {0.23480892, 0.50037575, 0.26630706},
//            {0.04802565, 0.5408284, -0.25837395},
//            {0.2660291, 0.5529919, 0.25671428},
//            {-0.0027829956, 0.8589586, 0.6085264},
//            {0.09651935, 0.5212971,  0.261826},
//            {0.2896213, 0.8597927, 0.19915773},
//            {0.116506316, 0.85993165, -0.08145551},
//            {0.2839288, 0.41863632, 0.08902351},
//            {0.2660291, 0.48960224, 0.25671428},
//            {0.107714586, 0.47952378, -0.12258158},
//            {0.23851536, 0.5123308, 0.31350744},
//            {-0.083173625, 0.5123308, 0.48020247},
//            {0.11919444, 0.5212971, 0.47249505},
//            {0.054837663, 0.4195399, 0.36409348},
//            {0.1944302, 0.6689979, 0.62611866},
//            {0.0065779896, 0.56320935, 0.54340184},
//            {0.14606933, 0.42582327, 0.48820195},
//            {0.20024918, 0.42176414, 0.30683565},
//            {0.24616863, 0.48960224, 0.1790096},
//            {-0.20739642, 0.87876785, 0.026853103},
//            {0.20448694,  0.752892, 0.25651512},
//            {-0.010562734, 0.5461803, 0.3330913},
//            {0.016508223, 0.63591295, 0.57201415},
//            {0.13099055, 0.37470832, 0.3494886},
//            {0.17279872, 0.5212971, 0.2422654},
//            {0.30960828,  0.979969, 0.2953177},
//            {0.14217313, 0.5212971, 0.2933593},
//            {0.23193738, 0.5113577, -0.22156295},
//            {0.16407025, 0.48049688, 0.2722154},
//            {-0.15692301, 0.9834444, 0.47714868},
//            {-0.024983712, 0.5212971, 0.4389104},
//            {0.124981806, 0.42211163, 0.36329687},
//            {0.21549243, 0.63028294, 0.5176441},
//            {0.26141185, 0.7670713, 0.19812873},
//            {-0.0066412394,  0.979969, 0.5345061},
//            {-0.07937863, 0.8688285, 0.0941684},
//            {0.08159237, 0.6187449, 0.08414414},
//            {0.15698627, 0.5212971, -0.28562537},
//            {0.24737038, 0.5212971, 0.18000537},
//            {0.13402654, 0.40668124, 0.34772938},
//            {-0.1944302, 1.0165293, 0.03774038},
//            {0.034344696, 0.9636351, 0.05340746},
//            {-0.02536321, 0.9784399, -0.018056955},
//            {0.088809185, 0.5212971, 0.3136037},
//            {0.27918506, 0.5123308, 0.60411173},
//            {-0.037064448, 0.53026336, 0.22571209},
//            {0.033022776, 0.5212971, 0.29715326},
//            {0.004427493, 0.37783617, 0.12935293},
//            {0.0316882, 1.0150002, 0.123444594},
//            {-0.2253594, 0.8688285, 0.03229674},
//            {0.2753901, 0.5113577, 0.5714831},
//            {0.13086404, 0.5212971, 0.42749202},
//            {0.24072912, 0.5212971, -0.20466773},
//            {0.045033928, 0.40508258, 0.34862557},
//            {0.18443671, 0.5212971, 0.09589443},
//            {0.14123702, 0.5312365, 0.3678443},
//            {-0.17279872, 0.8688285, 0.42159367},
//            {-0.23851536, 0.8598622, 0.01842208},
//            {0.010461533, 0.5212971, 0.2413692},
//            {-0.077734135, 0.53026336, 0.18438683},
//            {-0.10714534, 0.32570642, 0.10860734},
//            {0.015432976, 0.48007986, 0.27158472},
//            {0.091712356, 0.5943482, -0.10090658},
//            {-0.047880176, 1.0171549, 0.31417128},
//            {-0.20024918, 0.7692955, 0.02509388},
//            {-0.28363785, 0.4259345, 0.3424517},
//            {0.1646395, 0.40654224, 0.24134597},
//            {0.03261165, 0.4259345, 0.48737213},
//            {0.07008089, 0.9005233, 0.34497437},
//            {0.4623188, 0.61677086, 0.34162188},
//            {0.14541152, 1.0152782, -0.03455386},
//            {0.016508223, 0.40668124, 0.24008462},
//            {-0.08994135, 0.12608439, 0.11361948},
//            {0.38683638, 0.90059286, 0.29837143},
//            {-0.11941582, 0.82775027, 0.34995332},
//            {-0.07158624, 0.5212971, 0.4649005},
//            {0.13339405, 0.83713365, 0.26162684},
//            {0.18665045, 0.6212471, 0.13227391},
//            {0.21024267, 0.63632995, 0.20519884},
//            {0.26141185, 0.7670713, 0.53005826},
//            {-0.07937863, 0.8688285, 0.42609793},
//            {-0.0066412394,  0.979969, 0.2025766},
//            {0.08159237, 0.6187449, 0.41607365},
//            {-0.1912677, 0.42211163, -0.1973321},
//            {0.124981806, 0.42211163, 0.031367324},
//            {0.21549243, 0.63028294, 0.18571456},
//            {0.2229559, 0.6211081, 0.21787855},
//            {0.20739642, 0.5113577, 0.6370059},
//            {-0.010562734, 0.49641386, 0.0011617483},
//            {-0.090320855, 0.91032374, 0.19806236},
//            {0.08545062, 0.66809434, 0.4277244},
//            {-0.28363785, 0.6166597, 0.010522162},
//            {0.21100166, 0.6683029, 0.20868409},
//            {0.34287772, 0.9005233, 0.27387506},
//            {0.18525895, 0.72223973, 0.51545334},
//            {0.35467383, 0.8461695, 0.30431297},
//            {-0.09247135, 0.75247496, 0.5672676},
//            {0.21644118, 0.5212971, -0.08832644},
//            {0.12959905, 0.76887846, -0.033690847},
//            {0.1944302, 0.37359628, 0.29418916},
//            {0.30372605, 0.76755786, 0.27669644},
//            {0.14606933, 0.61677086, 0.15627243},
//            {-0.091712356, 0.79577744, 0.10090658},
//            {-0.08159237, 0.9662763, 0.24778539},
//            {0.1756513, 0.61649984, 0.21127315},
//            {0.11467207, 0.6363995, 0.15102793},
//            {0.11941582, 0.56237525, -0.018023778},
//            {-0.08026413, 0.18280151, 0.29637986},
//            {0.024287963, 0.40772384, 0.5171462},
//            {0.17343123, 0.53019387, 0.36864093},
//            {-0.1756513, 0.77362573, 0.4525859},
//            {0.15692301, 0.63591295, 0.18671036},
//            {0.14477903, 0.56300086, 0.15783249},
//            {0.21024267, 0.40626422, 0.5371283},
//            {0.1686875, 0.6354959, 0.12991722},
//            {0.022833215, 0.8688285, 0.50141275},
//            {-0.083173625, 0.53026336, 0.48020247},
//            {0.08159237, 0.42384928, 0.08414414},
//            {0.25964084, 0.53026336, 0.4591249},
//            {-0.0027829956, 0.87869835, 0.6085264},
//            {-0.050220422, 0.5529919, -0.09074953},
//            {0.17317823, 0.5634874, 0.21153869},
//            {0.16116075, 0.62805873, 0.11146193},
//            {-0.060024157, 0.8688285,  0.191364},
//            {0.1912677,  0.769643, 0.1973321},
//            { 0.170838, 0.66774684, 0.1314109},
//            {0.091712356, 0.44824603, -0.10090658},
//            {0.19670719, 0.6375811, 0.62153804},
//            {0.0050663175, 0.61703503, 0.21044332},
//            {-0.13402654, 0.9834444, 0.31612965},
//            {-0.090023585, 0.77335465, 0.19905813},
//            {0.03261165, 0.6166597, 0.1554426},
//            {0.026628207, 0.5122613, 0.033192955},
//            {0.024287961, 0.63487035, 0.18521668},
//            {0.09247135, 0.40494356, 0.42852104},
//            {0.20024918, 0.42176414, 0.63876516},
//            {0.24616863, 0.48960224, 0.5109391},
//            {0.0065779896, 0.56320935, 0.21147229},
//            {0.14606933, 0.42582327, 0.15627243},
//            {0.21232992, 0.40855792, -0.055166684},
//            {-0.124981806,  0.769643, 0.3005622},
//            {0.21232992, 0.63403624, 0.60869235},
//            {-0.054837663, 0.07200853, -0.032163966},
//            {0.090023585, 0.61677086, 0.13287139},
//            {0.06343965,  0.636469, 0.14780821},
//            {0.1547219,  0.498082, 0.08029375},
//            {-0.24616863, 0.9005233, 0.48484945},
//            {0.15654351, 0.48960224, 0.37979376},
//            { 0.170838, 0.37484738, 0.46334043},
//            {0.034724195, 0.5625143, 0.15670393},
//            {-0.15698627, 0.8688285, 0.28562537},
//            {0.23965387, 0.48960224, 0.46031988},
//            {0.21100168, 0.3742913, 0.5406136},
//            {-0.26754707, 0.8688285, 0.4069456},
//            {-0.107714586, 0.21553898, 0.12258158},
//            {0.0066412394, 0.6324376, 0.12935293},
//            {0.07684863, 0.5212971, 0.38646555},
//            {-0.15692301, 0.7542126, 0.14521916},
//            {0.13402654, 0.40668124, 0.67965895},
//            {0.090320855, 0.56279236, 0.13386717},
//            {0.05818991, 0.56237525, -0.017127557},
//            {-0.08026413, 0.8778643, 0.29637986},
//            {0.075646885, 0.87869835, 0.031898417},
//            {-0.016508223, 0.9834444, 0.42377445},
//            {0.0916491,  0.531167, 0.39818266},
//            {-0.026628207, 0.8778643, 0.29873657},
//            {0.012523481, 0.6225678, 0.11073169},
//            {0.29088628, 0.7592171, 0.18402171},
//            {0.091712356, 0.44824603, 0.23102295},
//            {-0.2817783, 0.5212971, 0.39791712},
//            { 0.170838, 0.66774684, -0.20051864},
//            {0.004427493,  0.664758, 0.12935293},
//            {-0.010461533, 0.8688285, 0.09056034},
//            {0.19683369, 0.82775027, 0.47987053},
//            {-0.07008089, 0.48960224, -0.013044839},
//            {-0.16510755, 0.9641077, 0.059183042},
//            {0.09247135, 0.63765055, 0.096591495},
//            {0.111762576, 0.63723356, 0.09055037},
//            {0.4623188, 0.42582327, 0.34162188},
//            {-0.28363785, 0.6166597, 0.3424517},
//            {0.21100168, 0.66830283, 0.5406136},
//            {-0.010562734, 0.49641386, 0.33309126},
//            {-0.090320855, 0.91032374, 0.52999187},
//            {-0.121819325, 0.37359628, -0.12822439},
//            {0.08545062, 0.66809434, 0.095794864},
//            {-0.035166945, 0.5104541, 0.42311054},
//            {-0.09247135,  0.985182, 0.23533802},
//            {-0.011005484, 0.4241273, 0.13406633},
//            {0.0065779896, 0.4793848, 0.54340184},
//            {0.52364594, 0.5312365, 0.19281788},
//            {-0.042187683, 0.4216946, 0.22046758},
//            {-0.14123702, 0.87876785, 0.29601476},
//            {-0.069258645, 0.6520384, 0.15457958},
//            {0.29664204, 0.63473135, 0.21183743},
//            {-0.038424317, 0.54395616, 0.13834822},
//            {0.37968916,  0.636469, 0.35008606},
//            {0.106006846, 0.9838614, 0.039234072},
//            {0.22179842, 0.54061985, 0.15653796},
//            {0.18443671, 0.5212971, -0.2360351},
//            {0.24072912, 0.5212971, 0.12726177},
//            {0.19670719, 0.40501308, 0.62153804},
//            {-0.16116074, 0.76206684, 0.22046758},
//            {-0.07912563, 0.5105236, 0.44760698},
//            {0.1756513, 0.42609435, -0.12065638},
//            {0.22179843, 0.50197434, 0.15653796},
//            {0.050220422, 0.83713365, 0.42267907},
//            {0.24616863, 0.5529919, 0.1790096},
//            {0.19974321, 0.5124003, -0.24742028},
//            {-0.0066412394, 0.7576879, 0.2025766},
//            {0.23965387, 0.5529919, 0.12839034},
//            {0.04802565, 0.5408284, 0.40548512},
//            {-0.0916491, 0.8589586, 0.2656764},
//            {0.24060263,  0.531167, 0.19786319},
//            {0.105247825, 0.7218227, 0.042719338},
//            {-0.01309273, 0.9641911, 0.060046043},
//            { 0.314099, 0.40334493, 0.098151565},
//            {0.13820103, 0.63521785, 0.6205754},
//            {0.2661556, 0.53026336, 0.17781463},
//            {-0.11467207, 0.7537261, 0.18090159},
//            {-0.038424317, 0.49863806, 0.13834822},
//            {0.12896655, 0.87876785, -0.051150344},
//            {0.24060263, 0.5114272, 0.19786319},
//            {0.18728295, 0.5113577, -0.2171151},
//            {0.116506316, 0.8777253, -0.08145551},
//            {0.39942312, 0.87779474, 0.31423768},
//            {0.2661556, 0.5123308, 0.17781466},
//            {0.23965387, 0.48960224, 0.12839034},
//            {0.21100168, 0.3742913, 0.20868409},
//            {-0.26754707, 0.8688285, 0.075016074},
//            {-0.15698627, 0.8688285, 0.6175549},
//            {0.0347242, 0.5625143, 0.48863345},
//            {0.25964084, 0.53026336, 0.1271954},
//            {0.08159237, 0.42384928, 0.41607365},
//            {0.23079889, 1.0156257, -0.0701699},
//            {-0.0027829956, 0.87869835, 0.27659687},
//            {0.17317823, 0.5634874, 0.54346824},
//            {0.48942775, 0.5634874, 0.28635558},
//            {-0.034344696, 0.61610365, 0.6104516},
//            {0.0916491, 0.5114272, -0.2656764},
//            {0.25964084, 0.5123308, 0.12719539},
//            {0.23193738, 0.5312365, 0.11036657},
//            {0.116000324, 0.7692955, 0.14087088},
//            {0.16407025, 0.48049688, -0.059714116},
//            {0.23193738, 0.5113577, 0.11036657},
//            {0.2229559, 0.42148608, 0.21787855},
//            {-0.0347242, 0.82761127, 0.1752256},
//            {0.18665047, 0.42134705, 0.13227391},
//            {0.047880176, 0.3729707, 0.6816173},
//            {-0.091712356, 0.9418796, 0.10090658},
//            {0.18728295, 0.5312365, 0.11481442},
//            {0.116506316, 0.8777253,  0.250474},
//            {0.18728295, 0.5113577, 0.11481442},
//            {0.39942312, 0.87779474, -0.017691841},
//            {0.2661556, 0.5123308, -0.15411487},
//            {-0.083173625, 0.5123308, 0.14827292},
//            {0.23851536, 0.5123308,  0.645437},
//            {0.11919444, 0.5212971, 0.14056551},
//            {-0.054837663, 0.97058564, -0.032163966},
//            {0.16315944, 0.5212971,  0.181592},
//            {0.05836701, 0.61671525, -0.016895207},
//            {-0.20024918, 0.9683615, 0.02509388},
//            {0.31182203, 0.72536755, 0.2953177},
//            {-0.030676205, 0.4073763, 0.53815734},
//            {0.14541152, 0.7223788, -0.03455386},
//            {0.32301724, 0.8688285, 0.4893637},
//            {   0.5029, 0.42134705, 0.36562037},
//            {0.14477903, 0.4795933, 0.15783249},
//            { 0.472793, 0.48960224, 0.45003006},
//            {0.22179843, 0.50197434, -0.17539157},
//            {0.1756513, 0.42609435, 0.21127315},
//            {0.17317823, 0.47910678, 0.2115387},
//            {0.08949861, 0.6355654, 0.6421509},
//            {0.01309273, 0.4259345,  0.603813},
//            {0.3416127, 0.41168568, -0.18402171},
//            {-0.1276383, 0.3564977, 0.09091549},
//            {-0.2580596, 0.4802189, -0.14883721},
//            {-0.07684863, 0.8688285,  0.609323},
//            {0.18285546, 0.48960224, 0.09566209},
//            {0.08545062, 0.37449983, 0.095794864},
//            {-0.08494462, 0.42996582, 0.09934651},
//            {0.21549243, 0.41231126, 0.5176441},
//            {-0.034344696, 0.42649055, 0.6104516},
//            {-0.23193738, 0.87876785, 0.5534925},
//            {0.19974318, 0.5301939, 0.08450925},
//            {-0.18728295, 0.87876785, 0.2171151},
//            {0.015432976, 0.5625143, 0.60351425},
//            {0.18285546, 0.5529919, 0.09566209},
//            {0.38683638, 0.83706415, -0.033558074},
//            {0.19974321, 0.5124003, 0.08450925},
//            {0.23965387, 0.5529919, -0.20353918},
//            {0.0154899005, 0.5212971, -0.15080555},
//            {0.15472192, 0.54451215, 0.08029375},
//            {0.24737035, 0.5212971, -0.15192415},
//            {0.15698627, 0.5212971, 0.04630417},
//            {-0.16116074, 0.06700407, 0.22046758},
//            {0.30081654, 0.91004574, 0.4375495},
//            {-0.13099055, 0.027176948, -0.017559057},
//            {-0.02536321, 0.7592171, 0.31387258},
//            {-0.11919444, 0.8688285,  0.191364},
//            {-0.01922797, 0.63215965, 0.6384664},
//            {0.29430178, 0.83713365, 0.2189739},
//            {0.15654351, 0.5529919, 0.047864236},
//            {0.090320855, 0.47980183, 0.46579668},
//            {0.054837663, 0.4195399,  0.696023},
//            {0.06343965,  0.636469, 0.47973776},
//            {0.1547219,  0.498082, 0.4122233},
//            {-0.24616863, 0.9005233, 0.15291992},
//            {0.15654351, 0.48960224, 0.047864236},
//            { 0.170838, 0.37484738, 0.13141091},
//            {0.010461533, 0.5212971, 0.57329875},
//            {0.14123704, 0.5312365, 0.035914775},
//            {-0.23851536, 0.8598622, 0.3503516},
//            {-0.17279872, 0.8688285, 0.08966412},
//            {0.15692301, 0.63591295, 0.51863986},
//            {0.17343123, 0.5301939, 0.036711406},
//            {-0.1756513, 0.77362573, 0.12065638},
//            {0.021947715, 0.5529919, 0.3849387},
//            {0.14123702, 0.5113577, 0.035914775},
//            {0.0027829956, 0.5114272, -0.27659687},
//            {0.17343123, 0.5124003, 0.036711406},
//            {0.011005484, 0.96599823, 0.5297927},
//            {0.012523481, 0.42002645, 0.11073169},
//            {0.1547219, 0.54451215, -0.2516358},
//            {0.0154899005, 0.5212971, 0.18112397},
//            {-0.05836701, 0.2691839, 0.016895207},
//            {0.060024157, 0.5212971, 0.14056551},
//            {-0.14123702, 0.8588891, 0.29601476},
//            {0.0050663175, 0.42555913, 0.21044332},
//            {0.016508223, 0.40668124, 0.57201415},
//            {0.4623188, 0.61677086, 0.009692338},
//            {0.07008089, 0.9005233, 0.013044839},
//            {0.03261165, 0.4259345, 0.1554426},
//            {0.0065779896, 0.4793848, 0.21147229},
//            {-0.09247135,  0.985182, 0.5672676},
//            {-0.011005484, 0.4241273, 0.46599588},
//            {-0.035166945, 0.5104541, 0.09118104},
//            {-0.08880919, 0.8688285, 0.35025534},
//            {-0.17343123, 0.1648689, 0.2952181},
//            {-0.07558363, 0.9005233, 0.27656367},
//            {0.0347242, 0.48007986, 0.15670392},
//            {0.090023585, 0.42582327, 0.13287139},
//            {0.024287963, 0.40772384, -0.14671284},
//            {0.11941582, 0.56237525, 0.6458353},
//            {0.090320855, 0.47980183, 0.13386717},
//            {0.15654351, 0.5529919, 0.37979376},
//            {-0.11919444, 0.8688285, 0.5232935},
//            {-0.01922797, 0.63215965, 0.3065369},
//            {0.09247135, 0.40494356, 0.096591495},
//            {0.024287963, 0.63487035, 0.5171462},
//            {0.31346652, 0.8589586, -0.11063211},
//            {0.111762576, 0.4053606, 0.09055037},
//            {-0.00594549, 0.5133734, 0.10495611},
//            {0.05836701, 0.61671525, 0.64696383},
//            {-0.23687088, 0.5212971, -0.07179636},
//            {-0.026628207, 0.8778643, 0.6306661},
//            {0.012523481, 0.6225678, 0.4426612},
//            {0.075646885, 0.87869835, 0.36382794},
//            {-0.016508223, 0.9834444, 0.0918449},
//            {0.09164911,  0.531167, 0.06625313},
//            {0.25964084, 0.5123308, -0.20473413},
//            {0.0916491, 0.5114272, 0.06625313},
//            {0.04802565, 0.5408284, 0.07355558},
//            {-0.0916491, 0.8589586, 0.59760594},
//            {0.2660291, 0.5529919, -0.07521523},
//            {0.105247825, 0.7218227, 0.37464887},
//            {0.24060263,  0.531167, 0.52979267},
//            {-0.01309273, 0.9641911, 0.39197555},
//            {0.23079889, 0.72203124, 0.26175964},
//            {0.22592865, 0.8273333, -0.03209759},
//            {0.0066412394, 0.63243765, 0.46128246},
//            {0.07684863, 0.5212971, 0.054536022},
//            {0.22377814,  0.985182, 0.26255625},
//            {-0.08949861, 0.9830968, 0.021708194},
//            {-0.14281829, 0.53019387, 0.12925336},
//            {0.021188717, 0.5212971, 0.051681425},
//            {-0.06343965, 0.28893763, 0.1841213},
//            {0.04802565, 0.50176585, 0.07355558},
//            {0.06887915, 0.8688285, 0.014040614},
//            {-0.159706, 0.5529919, 0.45003006},
//            {0.07558363, 0.5529919, 0.055365846},
//            {-0.06687412, 0.39451763, 0.033657648},
//            {-0.24009663, 0.4669432, -0.03797275},
//            {0.07552038, 0.8688285, 0.29322654},
//            {0.14123702, 0.5113577, 0.3678443},
//            {0.021947715, 0.5529919, 0.05300915},
//            {-0.12372312, 0.5212971, 0.3429164},
//            {0.07558363, 0.48960224, 0.055365846},
//            {0.08026412, 0.5303329, 0.03554965},
//            {-0.07058689, 0.48953274, -0.19952284},
//            {0.17501247, 0.87876785, 0.20187953},
//            {0.08026413, 0.5122613, 0.035549656},
//            {0.05818991, 0.4802189, 0.6467315},
//            {0.076152876,  0.575651, 0.53586704},
//            {-0.14477903, 0.91053224, 0.50602657},
//            {-0.021821218, 0.6452268, 0.06572205},
//            {0.021947715, 0.48960224, 0.05300915},
//            {0.30315676, 0.9641911, 0.10591872},
//            {0.0027829956,  0.531167, 0.05533265},
//            {0.026628207, 0.5303329, 0.033192955},
//            {-0.2753901, 0.18370508, 0.09237599},
//            {0.17343123, 0.5124003, -0.2952181},
//            {0.0027829956, 0.5114272, 0.055332657},
//            {-0.00594549, 0.5133734, -0.22697341},
//            {0.16315944, 0.5212971, 0.5135215},
//            {0.111762576, 0.4053606, -0.24137916},
//            {   0.5029, 0.42134705, 0.033690847},
//            {0.32301724, 0.8688285, 0.15743418},
//            {-0.20024918, 0.9683615, 0.35702342},
//            {-0.030676205, 0.4073763, 0.20622781},
//            {0.23712389, 0.5320706, 0.38221687},
//            {0.14756203, 0.7546296,  0.295882},
//            {-0.037254192, 0.41550857, 0.18754019},
//            {0.23465714, 0.77138066, -0.08182062},
//            {-0.0316882, 0.37512535, 0.20848492},
//            {0.2580596, 0.82775027, 0.48076674},
//            {0.024287963, 0.40772384, 0.1852167},
//            {0.090023585, 0.42582327, -0.19905813},
//            {0.11467207, 0.6363995, 0.48295745},
//            {0.06343965, 0.4061252, 0.14780822},
//            {0.0066412394, 0.41015655, 0.12935293},
//            {0.2246004, 0.8589586, -0.09971163},
//            {-0.04085943, 0.5113577, 0.25834075},
//            {0.037064448, 0.87779474, 0.43814698},
//            {0.1646395, 0.63605195, 0.5732755},
//            {-0.021947715, 0.83713365, 0.61084986},
//            {0.17501247, 0.8588891, 0.20187953},
//            {-0.01922797, 0.41043457, 0.3065369},
//            {0.1686875, 0.40709826, -0.20201232},
//            {0.010562734, 0.84394526, 0.33076778},
//            {0.3469257, 0.9827493, 0.37219256},
//            {0.11467207, 0.4061947, 0.15102793},
//            {0.07937863, 0.5212971, 0.23776111},
//            {0.15692301, 0.40668124, 0.18671036},
//            {0.1686875, 0.40709826, 0.12991722},
//            {-0.01922797, 0.41043457, -0.025392609},
//            {0.16407025, 0.56209725, 0.60414493},
//            {0.17501247, 0.8588891, -0.13004999},
//            {0.4709714, 0.54451215, 0.41760054},
//            {-0.116506316, 0.53019387, 0.41338503},
//            {-0.08144058, 0.54221845, 0.23158723},
//            {0.31182203, 1.0122894, -0.03661183},
//            {0.16116074, 0.41453546, 0.11146193},
//            {0.02536321, 0.41168568,  0.681916},
//            {0.24066588, 0.9005233, 0.2213306},
//            {-0.2661556, 0.8598622, 0.4860444},
//            {0.21024267, 0.40626422, 0.20519884},
//            {0.14477903, 0.56300086,  0.489762},
//            {0.18665047, 0.6212471, 0.46420342},
//            {-0.07158624, 0.5212971, 0.13297097},
//            {0.2896213, 0.8778643, 0.19915773},
//            {0.57589036, 0.5123308, 0.3706989},
//            {0.076152876, 0.4669432, -0.12799202},
//            {-0.07058689, 0.5530615, 0.1324067},
//            {-0.07058689, 0.48953274, 0.1324067},
//            {0.08026413, 0.5303329, -0.29637986},
//            {0.022833215, 0.8688285, 0.16948321},
//            {-0.08317362, 0.53026336, 0.14827292},
//            {0.1686875, 0.6354959, 0.46184677},
//            {0.084312126, 0.8588891, -0.055598192},
//            {-0.07912563, 0.5320706, 0.11567744},
//            {-0.07912563, 0.5105236, 0.115677446},
//            {0.19670719, 0.40501308, 0.2896085},
//            {0.24072912, 0.5212971, 0.45919132},
//            {-0.16116074, 0.76206684, 0.55239713},
//            {-0.15692301, 0.9834444, 0.14521916},
//            {-0.024983712, 0.5212971, 0.10698089},
//            {-0.28322673, 0.5212971, -0.13118851},
//            {0.16407025, 0.48049688, 0.60414493},
//            {0.2815253, 0.82761127, 0.3226687},
//            {-0.05818991, 0.9099067, 0.3490571},
//            {0.23480892, 0.50037575, 0.59823656},
//            {-0.026628207, 0.48960224, 0.10791028},
//            {0.40619087, 0.5689784, 0.27958423},
//            {0.28363785, 0.9641911, 0.32140738},
//            {0.03870894, 0.4061252, 0.2910026},
//            {-0.14281829, 0.5124003, 0.46118286},
//            {-0.005945491, 0.5292208, 0.10495612},
//            {-0.035166945, 0.5321401, 0.09118104},
//            {0.093293615, 0.76901746, 0.3838433},
//            {-0.026628207, 0.5529919, 0.1079103},
//            {0.2253594, 0.5529919,  0.629969},
//            {0.03116323, 0.7501187, 0.30169073},
//            {0.19705509, 0.8688285, 0.3065303},
//            {-0.15382376, 0.9772583, 0.31563178},
//            {0.01309273, 0.6166597,  0.603813},
//            {0.25280985, 0.75365657, 0.31377298},
//            {-0.011005483, 0.61846685, 0.13406634},
//            {0.27918506, 0.53026336, -0.059747316},
//            {-0.15926325, 0.5212971, 0.45159012},
//            {-0.037254192, 0.6270856, 0.18754019},
//            {0.116506316, 0.85993165,  0.250474},
//            {0.2896213, 0.8597927, -0.13277182},
//            {0.107714586, 0.47952378, 0.20934795},
//            {0.4732358, 0.5212971, 0.45159012},
//            {0.2660291, 0.48960224, -0.07521523},
//            {-0.076152876, 0.9231824, 0.45992154},
//            {-0.042187683, 0.62089956, 0.22046758},
//            {-0.01960747, 0.4078628, 0.28605688},
//            {0.33168247, 0.5625143, -0.10561997},
//            {0.08949861, 0.40702876, 0.31022134},
//            {0.119770005, 0.42613602, 0.64656556},
//            {-0.21100168, 0.32077146, 0.123245426},
//            {-0.124981806, 0.2729512, -0.031367324},
//            {-0.0067677395, 0.5212971, 0.32339895},
//            {0.13820103, 0.4073763, 0.28864592},
//            {0.3469257, 0.75490767, 0.37219256},
//            {0.21232992, 0.40855792, 0.27676284},
//            {-0.124981806,  0.769643, -0.031367324},
//            {0.08994135, 0.4736158, 0.21831004},
//            {0.10391958, 0.9815677, 0.44272763},
//            {0.2896213, 0.8778643, -0.13277182},
//            {-0.07158624, 0.5212971, -0.19895856},
//            {0.076152876, 0.4669432, 0.20393752},
//            {-0.322195, 0.5292208, 0.3929382},
//            {0.10771458, 0.56307036, 0.20934796},
//            {0.3479377, 1.0150002, 0.37444967},
//            {0.08994136, 0.5689784, 0.21831004},
//            {-0.1547219, 0.8920436, 0.5835653},
//            {-0.1547219, 0.84561336, 0.2516358},
//            {-0.14541152, 0.66774684, 0.03455386},
//            {0.2789953, 0.6270856, 0.3103541},
//            {-0.01960747, 0.6347313, 0.28605685},
//            {0.25964084, 0.5123308, 0.4591249},
//            {-0.034344696, 0.61610365, 0.27852207},
//            {0.050220422, 0.9005233, 0.42267907},
//            {-0.21024267, 0.75379556, 0.1267307},
//            {0.20024918, 0.62083006, 0.63876516},
//            {0.5285794, 0.63403624, 0.22113144},
//            {0.13099055, 0.6678859, 0.6814181},
//            {-0.29974127, 0.63591295, 0.25780967},
//            {0.124981806, 0.62048256, 0.69522643},
//            {0.3412332, 0.8688285, 0.27294564},
//            {0.035166945, 0.85798544,  0.572678},
//            {0.119770005, 0.6164582, 0.64656556},
//            {0.2839288, 0.6239579, 0.08902349},
//            {0.091712356, 0.5943482, 0.56295246},
//            {-0.021947715, 0.20546056, 0.27892038},
//            {0.045033928, 0.6375116, 0.6805551},
//            {0.15508877, 0.9755901, 0.2774267},
//            {0.02536321, 0.6309085,  0.681916},
//            {0.011005484, 0.96599823, 0.19786319},
//            {0.012523481, 0.42002645, 0.4426612},
//            {0.047880176, 0.6696235, 0.6816173},
//            {-0.21100168, 1.0158342, 0.123245426},
//            {-0.124981806,  0.968014, -0.031367324},
//            {0.054837663, 0.62305427,  0.696023},
//            {0.17343123, 0.53019387, -0.2952181},
//            {0.09089011, 0.8688285, 0.13366802},
//            {-0.2661556, 0.87779474, 0.4860444},
//            {0.15382376, 0.41286728, 0.016297754},
//            {0.03870894,  0.636469, 0.62293214},
//            {0.29196155, 0.9824017, 0.35118148},
//            {-0.22744031, 0.5212971, 0.18429059},
//            {0.27918506, 0.53026336, 0.60411173},
//            {0.2753901, 0.5312365, 0.5714831},
//            {0.21644118, 0.5212971, 0.5755326},
//            {0.26754707, 0.5212971,  0.588843},
//            {0.22622591, 0.77335465, -0.033093374},
//            {-0.35141644, 0.5321401, 0.40671325},
//            {0.2253594, 0.5212971, 0.63156235},
//            {0.2660291, 0.5529919, 0.58864385},
//            {0.4050587, 0.5212971, 0.18429059},
//            {0.2660291, 0.48960224, 0.58864385},
//            {-0.24737035, 0.8688285, 0.48385367},
//            {-0.22630815, 0.5689784, 0.27958423},
//            {0.2253594, 0.48960224,  0.629969},
//            {0.51067966, 0.6689979, 0.20370515},
//            {0.22453715, 0.79577744, 0.06505818},
//            {0.038424317, 0.8461695, 0.1935813},
//            {-0.06343965, 0.75365657, 0.1841213},
//            {0.20739642, 0.5312365, 0.6370059},
//            {0.23851536, 0.53026336,  0.645437},
//            {0.14345077, 0.8688285, 0.40823016},
//            {-0.0027829956, 0.8589586, 0.27659687},
//            {-0.037064448, 0.5123308, -0.106217444},
//            {0.09651935, 0.5212971, 0.59375554},
//            {0.32289076, 0.41015655, 0.03661183},
//            {0.27918506, 0.5123308, -0.059747316},
//            {0.14217313, 0.5212971, 0.62528884},
//            {-0.022833215, 0.66566163, 0.6127087},
//            {0.076595634, 0.9005233, 0.2943551},
//            {-0.03870894, 0.9840004, 0.37285647},
//            {0.16510755, 0.42601788,  0.604676},
//            {0.24009663, 0.8144746, 0.03797275},
//            {0.124981806, 0.42211163, 0.69522643},
//            {-0.025869206, 0.67539257,  0.404091},
//            {-0.010562734, 0.5461803, 0.66502076},
//            {0.026628207, 0.5303329, -0.29873657},
//            {0.13099055, 0.37470832, 0.6814181},
//            {0.2578825, 0.9642467, 0.14906956},
//            {0.15382376, 0.41286728, 0.68015677},
//            {0.2712156,  0.985043, 0.18266082},
//            {0.11941582, 0.4802189, 0.6458353},
//            {0.045033928, 0.40508258, 0.6805551},
//            {0.015432976, 0.48007986, 0.60351425},
//            {-0.20024918, 0.7692955, 0.35702342},
//            {-0.28363785, 0.4259345, 0.010522162},
//            {0.1646395, 0.40654224, 0.5732755},
//            {0.010562734, 0.14888246, -0.0011617483},
//            {0.21024267, 0.63632995, 0.5371283},
//            {0.1944302, 0.37359628, 0.62611866},
//            {-0.08026413, 0.16472986, 0.29637986},
//            {0.27406183, 0.4216946, 0.2774267},
//            {0.14606933, 0.61677086, 0.48820195},
//            {0.39942312, 0.8598622, 0.31423768},
//            {-0.091712356, 0.79577744, 0.43283612},
//            {-0.08159237, 0.9662763, 0.5797149},
//            {0.1756513, 0.6164998, 0.54320264},
//            {0.16116074, 0.6280587, 0.44339147},
//            {-0.060024157, 0.8688285, 0.5232935},
//            {0.1912677,  0.769643, 0.52926165},
//            { 0.170838, 0.66774684, 0.46334043},
//            {-0.090023585, 0.77335465, 0.5309877},
//            {0.03261165, 0.6166597, 0.48737213},
//            {0.026628207, 0.5122613, 0.3651225},
//            {0.05836701, 0.42587888, -0.016895207},
//            {0.40796188, 0.5943482, 0.26687133},
//            {0.090320855, 0.56279236, 0.46579668},
//            {-0.08026413, 0.8778643, 0.6283094},
//            {0.29088628, 0.7592171, 0.5159513},
//            {0.091712356, 0.44824603, 0.56295246},
//            {0.004427493,  0.664758, 0.46128246},
//            {-0.010461533, 0.8688285, 0.42248988},
//            {-0.16510755, 0.9641077, 0.39111257},
//            {0.09247135, 0.6376506, 0.42852104},
//            {0.4623188, 0.42582327, 0.009692338},
//            {0.111762576, 0.6372336, 0.42247987},
//            {-0.042187683, 0.4216946, 0.55239713},
//            {-0.14123702, 0.87876785, 0.6279443},
//            {0.07912563, 0.16299222, 0.21625209},
//            {0.37968916,  0.636469, 0.01815654},
//            {-0.038424317, 0.54395616, 0.47027776},
//            {0.22179843, 0.54061985, 0.48846748},
//            {0.106006846, 0.9838614, 0.3711636},
//            {-0.01309273, 0.7734659, 0.39197555},
//            {-0.23079889, 0.66809434, 0.0701699},
//            {-0.030676205, 0.63521785, 0.53815734},
//            {0.047880176, 0.6696235, 0.017758248},
//            {0.33168247, 0.48007986, 0.22630955},
//            {0.19252639, 0.5212971, 0.48690742},
//            {0.13402654, 0.40668124, 0.015799856},
//            {0.24737035, 0.5212971, 0.5119349},
//            {-0.15217926, 0.56209725, 0.22567888},
//            {0.15932651, 0.9834444, 0.35267514},
//            {-0.16510755, 0.2690449, 0.059183042},
//            {0.22179843, 0.50197434, 0.48846748},
//            {0.24616863, 0.5529919, 0.5109391},
//            {-0.0066412394, 0.7576879, 0.5345061},
//            {0.23965387, 0.5529919, 0.46031988},
//            {0.09445107, 0.8881512, -0.0094268},
//            {-0.038424317, 0.49863806, 0.47027776},
//            {-0.26141185, 0.4195399, -0.19812873},
//            {0.24060263, 0.5114272, 0.52979267},
//            {0.2661556, 0.5123308, 0.50974417},
//            {0.116000324, 0.7692955, 0.47280043},
//            {0.23193738, 0.5312365, 0.4422961},
//            {-0.107714586, 0.8270552, 0.12258158},
//            {-0.024983712, 0.5212971, -0.22494864},
//            {-0.28322673, 0.5212971, 0.20074102},
//            {0.23193738, 0.5113577, 0.4422961},
//            {-0.0067677395, 0.5212971, -0.00853058},
//            {0.2229559, 0.42148608, 0.5498081},
//            {-0.0347242, 0.82761127, 0.5071551},
//            {0.18665047, 0.42134705, 0.46420342},
//            {-0.091712356, 0.9418796, 0.43283612},
//            {0.18728295, 0.5312365, 0.44674394},
//            {0.18728295, 0.5113577, 0.44674394},
//            {0.14477903, 0.4795933,  0.489762},
//            {-0.15217926, 0.56209725, -0.106250644},
//            {-0.13086404, 0.8688285,  0.236367},
//            {0.15114197, 0.9641077, 0.43871126},
//            {0.1756513, 0.42609435, 0.54320264},
//            {-0.1007571, 0.41231126, -0.0197498},
//            {0.17317823, 0.47910678, 0.54346824},
//            {0.19670719, 0.40501308, -0.04232101},
//            {0.045033928, 0.40508258, 0.016696054},
//            {0.18443671, 0.5212971, 0.42782396},
//            {-0.30372605, 0.42002645, 0.05523307},
//            {-0.07912563, 0.5105236, -0.21625209},
//            {-0.116000324, 0.42176414, -0.14087088},
//            {0.18285546, 0.5529919, 0.42759162},
//            {0.19974321, 0.5124003, 0.41643876},
//            {-0.07583663, 0.3706075, 0.13987511},
//            {0.15698627, 0.5212971, 0.3782337},
//            {0.17343123, 0.5124003, 0.36864093},
//            {0.0154899005, 0.5212971, 0.5130535},
//            {0.060024157, 0.5212971, 0.47249505},
//            {-0.14123702, 0.8588891, 0.6279443},
//            {0.29974127, 0.9834444, 0.40604937},
//            {-0.30075958, 0.5212971, 0.31677032},
//            {0.42396408, 0.56307036, 0.28854635},
//            {0.0050663175, 0.42555913, 0.5423728},
//            {-0.07558363, 0.9005233, 0.6084932},
//            {-0.08880919, 0.8688285, 0.018325826},
//            {0.0347242, 0.48007986, 0.48863345},
//            {0.0916491, 0.5114272, 0.39818266},
//            {0.045033928, 0.6375116, 0.016696054},
//            {0.034344696, 0.7740219, 0.38533697},
//            {0.034471195, 0.5212971, 0.4319067},
//            { 0.322195, 0.86090475, -0.06100865},
//            {-0.02536321, 0.0641543, -0.018056955},
//            {-0.0316882, 0.66746885, 0.54041445},
//            {-0.08949861, 0.9830968, 0.35363773},
//            {-0.14281829, 0.53019387, 0.46118286},
//            {0.021188717, 0.5212971, 0.38361096},
//            {0.04802565, 0.50176585, 0.40548512},
//            {0.48031977, 0.56209725, -0.106250644},
//            {0.06887915, 0.8688285, 0.34597012},
//            {-0.07008089, 0.5529919, -0.013044839},
//            {-0.159706, 0.5529919, 0.118100524},
//            {0.07558363, 0.5529919, 0.38729537},
//            {-0.12372312, 0.5212971, 0.01098687},
//            {0.07558363, 0.48960224, 0.38729537},
//            {-0.2229559, 0.07395469, 0.11405098},
//            {-0.15382376, 0.06533589, -0.016297754},
//            {0.32131582, 0.42555913, 0.28745097},
//            {0.08026413, 0.5303329, 0.36747918},
//            {0.14281829, 0.85993165, 0.20267616},
//            {-0.024287963, 0.9824017, 0.47864237},
//            {-0.14477903, 0.91053224, 0.17409705},
//            {-0.021821218, 0.6452268, 0.39765155},
//            {0.021947715, 0.48960224, 0.3849387},
//            {0.2246004, 0.87869835, 0.2322179},
//            {-0.24060263, 0.8589586, 0.13406634},
//            {-0.20853493, 0.56307036, 0.28854635},
//            {0.0027829956,  0.531167, 0.3872622},
//            {0.3317394, 0.5212971, 0.31677032},
//            {0.30315676, 0.9641911, 0.43784824},
//            {0.026628207, 0.5303329, 0.3651225},
//            {0.0027829956, 0.5114272, 0.3872622},
//            {-0.037254192, 0.41550857, 0.5194697},
//            {-0.1686875, 0.28796455, 0.20201232},
//            {-0.1646395, 0.059010845, 0.090583555},
//            {0.2580596, 0.82775027, 0.14883721},
//            {-0.0316882, 0.37512535, 0.54041445},
//            {0.29702154, 0.63215965, 0.19135737},
//            {0.06343965, 0.4061252, 0.47973776},
//            {0.11467207, 0.4061947, 0.48295745},
//            {0.07937863, 0.5212971, 0.56969064},
//            {0.3469257, 0.9827493, 0.040263046},
//            {0.39537513,  0.879602, 0.2816422},
//            {-0.19670719, 0.7525445, 0.37425053},
//            {-0.04802565, 0.19329698, 0.25837395},
//            {0.15692301, 0.40668124, 0.51863986},
//            {-0.116506316, 0.53019387, 0.08145551},
//            {0.16116074, 0.41453546, 0.44339147},
//            {0.4709714, 0.54451215, 0.085671015},
//            { 0.159706, 0.9005233,  0.213829},
//            {-0.054837663, 0.27552286, -0.032163966},
//            {-0.07058689, 0.5530615, 0.4643362},
//            {-0.1740764, 0.5212971, 0.20453496},
//            {-0.07058689, 0.48953274, 0.4643362},
//            {0.45842263, 0.5212971, 0.20453496},
//            {-0.07912563, 0.5320706, 0.44760698},
//            {0.03870894, 0.4061252, 0.62293214},
//            {-0.14281829, 0.5124003, 0.12925336},
//            {-0.00594549, 0.5292208, 0.43688565},
//            {-0.035166945, 0.5321401, 0.42311054},
//            {-0.011005484, 0.61846685, 0.46599588},
//            {-0.15926325, 0.5212971, 0.11966059},
//            {-0.037254192, 0.6270856, 0.5194697},
//            {0.40796188, 0.44824603, 0.26687133},
//            {0.32671103, 0.5212971, 0.2565251},
//            {0.4732358, 0.5212971, 0.11966059},
//            {-0.30578798, 0.5212971, 0.2565251},
//            {-0.22453715, 0.44824603, 0.26687133},
//            {0.107714586, 0.47952378, 0.54127747},
//            {-0.076152876, 0.9231824, 0.12799202},
//            {-0.042187683, 0.62089956, 0.55239713},
//            {-0.01960747, 0.4078628, 0.6179864},
//            {0.12959905, 0.9687785, -0.033690847},
//            {0.08949861, 0.40702876, 0.6421509},
//            {-0.03870894, 0.75365657, 0.37285647},
//            {-0.010562734, 0.49641386, 0.66502076},
//            {0.13820103, 0.4073763, 0.6205754},
//            {0.21232992, 0.40855792, 0.60869235},
//            {0.21232992, 0.63403624, -0.055166684},
//            {-0.0154899005, 0.1737657, 0.15080555},
//            {0.3469257, 0.75490767, 0.040263046},
//            {0.08994135, 0.4736158, 0.55023956},
//            {0.10391958, 0.9815677, 0.11079808},
//            {0.13339405, 0.83713365, -0.07030268},
//            {0.076152876, 0.4669432, 0.53586704},
//            {0.30075958, 0.8688285, 0.34708872},
//            {-0.322195, 0.5292208, 0.06100865},
//            {0.107714586, 0.56307036, 0.54127747},
//            {0.43601954, 0.6164582, 0.18325828},
//            {0.08994135, 0.5689784, 0.55023956},
//            {0.3479377, 1.0150002, 0.042520165},
//            {0.056608666, 0.87779474, 0.29316014},
//            {-0.1547219, 0.84561336, 0.5835653},
//            {-0.1547219, 0.8920436, 0.2516358},
//            {-0.14541152, 0.66774684, 0.3664834},
//            {-0.01960747, 0.63473135, 0.6179864},
//            {0.19670719, 0.6375811, -0.04232101},
//            {0.10885309, 0.8588891, 0.13911167},
//            {0.02536321, 0.6309085, 0.018056955},
//            {-0.012523481, 0.76755786, 0.55312735},
//            {0.077734135, 0.8598622, 0.14754269},
//            {0.1944302, 0.6689979, -0.03774038},
//            {0.35350367, 0.76303995, 0.021575417},
//            {0.20024918, 0.62083006, -0.02509388},
//            {0.054837663, 0.62305427, 0.032163966},
//            {-0.022833215, 0.5212971, -0.16948321},
//            {0.16510755, 0.6165763, -0.059183042},
//            {0.13402654, 0.63591295, 0.015799856},
//            { 0.480889, 0.40654224, -0.0753812},
//            {0.3584372, 0.76922596, 0.38643235},
//            {-0.105247825, 0.3742913, -0.042719338},
//            {0.13099055, 0.6678859, 0.017559057},
//            {0.15382376, 0.62972695, 0.016297754},
//            {-0.08545062, 0.32056296, 0.23613468},
//            {0.2229559, 0.6211081, -0.11405098},
//            {-0.35467383, 0.54395616, 0.027616538},
//            {0.124981806, 0.62048256, 0.031367324},
//            {0.09089011, 0.9005233, 0.13207476},
//            {0.33908272, 0.72446394, 0.44674394},
//            {0.04085943, 0.87876785, 0.4055183},
//            {-0.15932651, 0.63591295, -0.020745596},
//            {0.13820103, 0.63521785, -0.043283608},
//            {0.2253594, 0.5529919, -0.03389001},
//            {0.01922797,  0.757966, 0.025392609},
//            {0.01309273, 0.6166597, -0.060046043},
//            {-0.27754056, 0.4061252, 0.20689169},
//            {0.01309273, 0.4259345, -0.060046043},
//            {0.08949861, 0.6355654, -0.021708194},
//            {-0.14477903, 0.13206193, 0.17409705},
//            {0.016508223, 0.63591295, -0.0918449},
//            { -0.15161, 0.40654224, 0.25654835},
//            {0.40170014, 0.37449983, 0.40209946},
//            {0.21549243, 0.41231126, -0.14621496},
//            {0.5272512, 0.3742913, 0.2892102},
//            {-0.034344696, 0.42649055, -0.05340746},
//            {-0.35350367, 0.41550857, 0.3103541},
//            {0.015432976, 0.5625143, -0.06034479},
//            {-0.045033928, 0.75261396, 0.31523347},
//            {0.2753901, 0.5312365, -0.09237599},
//            {0.13086404, 0.5212971, -0.236367},
//            {0.2753901, 0.5113577, -0.09237599},
//            {0.23480892, 0.54221845, -0.06562248},
//            {0.4412313, 0.62048256, -0.1973321},
//            {0.1964795, 0.77366745, 0.48060074},
//            {0.26754707, 0.5212971, -0.075016074},
//            {-0.21549243, 0.9778143, 0.14621496},
//            {0.5164987, 0.62083006, -0.14087088},
//            {0.2253594, 0.5212971, -0.03229674},
//            {0.23480892, 0.50037575, -0.06562248},
//            {-0.08994135, 0.8211472, 0.11361948},
//            {0.2253594, 0.48960224, -0.03389001},
//            {0.20739642, 0.5312365, -0.026853103},
//            {-0.170838, 1.0152782, 0.20051864},
//            {0.23851536, 0.53026336, -0.01842208},
//            {0.23851536, 0.5123308, -0.01842208},
//            {0.24616863, 0.48960224, -0.15291992},
//            {0.20024918, 0.42176414, -0.02509388},
//            {-0.26754707, 0.1737657, 0.075016074},
//            {0.09651935, 0.5212971, -0.07010352},
//            {0.17279872, 0.5212971, -0.08966412},
//            {0.30960828,  0.979969, -0.03661183},
//            {0.39309815, 0.5212971, 0.44335827},
//            {0.14217313, 0.5212971, -0.0385702},
//            {0.16510755, 0.42601788, -0.059183042},
//            {-0.17343123, 0.8777253, 0.6271477},
//            {0.119770005, 0.42613602, -0.01729353},
//            {0.11941582, 0.4802189, -0.018023778},
//            {0.076595634, 0.83713365, 0.2943551},
//            {0.08880919, 0.5212971, -0.018325826},
//            {-0.14123702, 0.1638263, 0.29601476},
//            {-0.02536321, 0.9784399, 0.31387258},
//            {0.2896213, 0.5529919,  0.389984},
//            {-0.037064448, 0.53026336, -0.106217444},
//            {0.004427493, 0.37783617, -0.2025766},
//            {0.033022773, 0.5212971, -0.03477626},
//            {0.02536321, 0.41168568, 0.018056955},
//            {0.3327577, 0.40668124, -0.07411986},
//            {-0.1547219, 0.1505506, 0.2516358},
//            {0.04870242, 0.8688285, 0.4228782},
//            {0.047880176, 0.3729707, 0.017758248},
//            {-0.105247825, 0.66830283, -0.042719338},
//            {-0.030676205, 0.63521785, -0.12570171},
//            {0.076152876,  0.575651, -0.12799202},
//            {0.05818991, 0.4802189, -0.017127557},
//            {0.010461533, 0.5212971, -0.09056034},
//            {0.1646395, 0.40654224, -0.090583555},
//            {-0.05818991, 0.21484388, 0.017127557},
//            {0.38683638, 0.90059286, -0.033558074},
//            {0.14541152, 1.0152782, 0.29737568},
//            {0.016508223, 0.40668124, -0.0918449},
//            {0.18665047, 0.6212471, -0.1996556},
//            {0.21024267, 0.63632995, -0.1267307},
//            {-0.1912677, 0.42211163, 0.13459744},
//            {0.21549243, 0.63028294, -0.14621496},
//            {-0.16407025, 0.2145659, 0.059714116},
//            {0.34287772, 0.9005233, -0.058054473},
//            {0.21100168, 0.66830283, -0.123245426},
//            {-0.20448694, 0.6372336, 0.40734392},
//            {0.28322673, 0.8688285, 0.13118851},
//            {0.14606933, 0.61677086, -0.17565711},
//            {0.30372605, 0.76755786, -0.05523307},
//            {0.04085943, 0.8588891, 0.4055183},
//            {0.042187683, 0.76922596, 0.11146193},
//            {0.1756513, 0.6164998, -0.12065638},
//            {0.3584372,  0.968431, 0.38643235},
//            {0.15692301, 0.63591295, -0.14521916},
//            {-0.045033928, 0.05755118, -0.016696054},
//            {0.14477903, 0.56300086, -0.17409705},
//            {0.1686875, 0.6354959, -0.20201232},
//            {0.33547747,  0.757966, 0.14057215},
//            {-0.050220422, 0.5529919, 0.24117999},
//            {0.17317823, 0.5634874, -0.12039083},
//            {0.16116074, 0.6280587, -0.22046758},
//            {-0.0050663175, 0.26950365, 0.1214862},
//            {0.35059422, 0.7740219, 0.1125573},
//            {0.10885309, 0.8588891, 0.4710412},
//            {0.0050663175, 0.61703503, -0.1214862},
//            {-0.13402654, 0.9834444, -0.015799856},
//            {0.026628207, 0.5122613, -0.29873657},
//            {0.03261165, 0.6166597, -0.17648691},
//            {0.024287963, 0.63487035, -0.14671284},
//            {0.14606933, 0.42582327, -0.17565711},
//            {0.0065779896, 0.56320935, -0.12045723},
//            {-0.15698627, 0.1737657, 0.28562537},
//            {0.090023585, 0.61677086, -0.19905813},
//            {-0.32681224, 0.49641386, -0.1671265},
//            {0.22744031, 0.8688285, 0.14763893},
//            {0.1547219,  0.498082, -0.2516358},
//            {0.06343965,  0.636469, -0.1841213},
//            {-0.07684863, 0.1737657, 0.27739352},
//            {0.0347242, 0.5625143, -0.1752256},
//            {-0.21024267, 0.28879857, 0.1267307},
//            {-0.024287963, 0.06019244, 0.14671284},
//            {0.16242574, 0.7603987, 0.18226251},
//            {0.0066412394, 0.63243765, -0.2025766},
//            {0.012523481, 0.6225678, -0.22119783},
//            {-0.04802565, 0.8883597, 0.25837395},
//            {0.004427493,  0.664758, -0.2025766},
//            {0.075646885, 0.8589586, 0.36382794},
//            {-0.07008089, 0.48960224, 0.31888467},
//            {-0.13086404, 0.1737657,  0.236367},
//            {0.09247135, 0.6376506, -0.23533802},
//            {0.111762576, 0.6372336, -0.24137916},
//            {0.08159237, 0.6187449, -0.24778539},
//            {-0.121819325, 0.37359628, 0.20370515},
//            {0.08545062, 0.66809434, -0.23613468},
//            {0.050093923, 0.87779474, 0.34377944},
//            {-0.011005484, 0.4241273, -0.19786319},
//            {-0.042187683, 0.4216946, -0.11146193},
//            {0.116000324, 0.9683615, 0.14087088},
//            {0.29664204, 0.63473135, -0.1200921},
//            {-0.038424317, 0.54395616, -0.1935813},
//            {0.22179843, 0.54061985, -0.17539157},
//            {0.19252639, 0.5212971, -0.17695165},
//            {0.1740764, 0.8688285, 0.12739457},
//            {-0.09971979, 0.64392006, 0.11859842},
//            {0.30967152, 0.91074073, 0.045507535},
//            {0.24616863, 0.5529919, -0.15291992},
//            {0.24060263,  0.531167, -0.13406634},
//            {0.2661556, 0.53026336, -0.15411487},
//            {0.24060263, 0.5114272, -0.13406634},
//            {0.12896655, 0.87876785, 0.28077918},
//            {-0.038424317, 0.49863806, -0.1935813},
//            {0.21100168, 0.3742913, -0.123245426},
//            {0.23965387, 0.48960224, -0.20353918},
//            {0.25964084, 0.53026336, -0.20473413},
//            {-0.11467207, 0.0586633, 0.18090159},
//            {0.037254192, 0.27955425, 0.14438935},
//            {0.23193738, 0.5312365, -0.22156295},
//            {0.18665047, 0.42134705, -0.1996556},
//            {-0.08949861, 0.7545602, 0.021708194},
//            {0.4078986, 0.5114272, 0.43164113},
//            {0.18728295, 0.5312365, -0.2171151},
//            {0.11919444, 0.5212971, -0.191364},
//            {-0.083173625, 0.5123308, -0.1836566},
//            {0.14477903, 0.4795933, -0.17409705},
//            {0.17147048, 0.91053224, 0.32379726},
//            {-0.111762576, 0.28970218, 0.24137916},
//            {0.17317823, 0.47910678, -0.12039083},
//            {-0.2580596, 0.4802189, 0.18309233},
//            {-0.24616863, 0.20546056, 0.15291992},
//            {0.08545062, 0.37449983, -0.23613468},
//            {0.18285546, 0.48960224, -0.23626743},
//            {0.19974321, 0.53019387, -0.24742028},
//            {-0.1756513, 0.07856296, 0.12065638},
//            {0.18285546, 0.5529919, -0.23626743},
//            {0.38683638, 0.83706415, 0.29837143},
//            {0.15654351, 0.5529919, -0.28406528},
//            {0.29430178, 0.83713365, -0.112955615},
//            { 0.170838, 0.37484738, -0.20051864},
//            {0.15654351, 0.48960224, -0.28406528},
//            {0.47007325, 0.41286728, 0.14966701},
//            {0.17147048, 0.8271248, 0.32379726},
//            {0.50353247, 0.5113577, 0.38307986},
//            {-0.14606933, 0.2692395, 0.17565711},
//            {0.14123702, 0.5312365, -0.29601476},
//            {0.07008089, 0.83713365, 0.013044839},
//            {-0.090320855, 0.21526095, 0.19806236},
//            {-0.23465714, 0.42384928, 0.41375014},
//            {0.14123702, 0.5113577, -0.29601476},
//            {0.012523481, 0.42002645, -0.22119783},
//            {-0.1944302, 0.32146654, 0.03774038},
//            {0.034344696, 0.26857227, 0.05340746},
//            {-0.02536321, 0.2833771, -0.018056955},
//            {0.060024157, 0.5212971, -0.191364},
//            {-0.107714586, 0.91060174, 0.4545111},
//            {0.0050663175, 0.42555913, -0.1214862},
//            {0.03261165, 0.4259345, -0.17648691},
//            {-0.23598538, 0.5122613, 0.13041511},
//            {-0.31182203, 0.37783617, 0.36854136},
//            {0.0065779896, 0.4793848, -0.12045723},
//            {0.050093923, 0.87779474, 0.011849892},
//            {-0.012523481, 0.07249506, 0.22119783},
//            {0.0347242, 0.48007986, -0.1752256},
//            {0.090320855, 0.47980183, -0.19806236},
//            {-0.18665047, 0.27371573, 0.1996556},
//            {0.5159927, 0.5124003, 0.41338503},
//            {0.08159237, 0.42384928, -0.24778539},
//            {0.31346652, 0.8589586, 0.22129741},
//            {0.09247135, 0.40494356, -0.23533802},
//            {0.0916491,  0.531167, -0.2656764},
//            {0.034471195, 0.5212971, -0.23195235},
//            {0.22592865, 0.8273333, 0.29983196},
//            {0.16242574, 0.7603987, 0.51419204},
//            {0.07684863, 0.5212971, -0.27739352},
//            {0.22377814,  0.985182, -0.069373265},
//            {-0.0065779896, 0.21567799, 0.12045723},
//            {0.27782518, 0.49863806, 0.027616538},
//            {0.021188717, 0.5212971, -0.2802481},
//            {0.26836935, 0.7205021, 0.51565254},
//            {0.04802565, 0.50176585, -0.25837395},
//            {-0.24009663, 0.4669432, 0.2939568},
//            {0.07552038, 0.8688285, -0.038702987},
//            {0.46102855, 0.56300086, 0.3400618},
//            {0.07558363, 0.5529919, -0.27656367},
//            {0.14281829, 0.8777253, 0.20267616},
//            {0.20853493, 0.91060174, 0.04338319},
//            {0.021947715, 0.5529919, -0.27892038},
//            {0.14307128, 0.8266382, 0.045573927},
//            {-0.22622591, 0.42582327, 0.3650229},
//            {0.07558363, 0.48960224, -0.27656367},
//            {0.43092155, 0.6363995, 0.34686637},
//            {0.08026413, 0.5122613, -0.29637986},
//            {0.17501247, 0.87876785, -0.13004999},
//            {0.021947715, 0.48960224, -0.27892038},
//            {0.0027829956,  0.531167, -0.27659687},
//            {-0.030676205, 0.4073763, -0.12570171},
//            {0.011005484, 0.7716587, 0.5297927},
//            {-0.037254192, 0.41550857, -0.14438935},
//            {0.55476487, 0.53026336, 0.18438683},
//            {0.23465714, 0.77138066, 0.2501089},
//            {0.2896213, 0.48960224,  0.389984},
//            {-0.0316882, 0.37512535, -0.123444594},
//            {-0.39537513, 0.5105236, 0.050287317},
//            {-0.14541152, 0.37484738, 0.3664834},
//            {0.06343965, 0.4061252, -0.1841213},
//            {-0.07937863, 0.1737657, 0.0941684},
//            {0.19683369, 0.9099067, 0.14794098},
//            {0.2246004, 0.8589586, 0.2322179},
//            {0.0066412394, 0.41015655, -0.2025766},
//            {0.17804848, 0.9827493, 0.12268115},
//            {0.010562734, 0.84394526, -0.0011617483},
//            {0.07937863, 0.5212971, -0.0941684},
//            {0.11467207, 0.4061947, -0.18090159},
//            {0.3509737, 0.5625143, 0.009260837},
//            {0.15692301, 0.40668124, -0.14521916},
//            {-0.08144058, 0.54221845, -0.10034229},
//            {0.31182203, 1.0122894, 0.2953177},
//            {0.16116074, 0.41453546, -0.22046758},
//            {0.12896655, 0.8588891, 0.28077918},
//            {0.21024267, 0.40626422, -0.1267307},
//            {-0.07058689, 0.5530615, -0.19952284},
//            {0.33547747,  0.757966, 0.47250167},
//            {0.084312126, 0.8588891, 0.27633134},
//            {-0.083173625, 0.53026336, -0.1836566},
//            {0.056608666, 0.8598622, 0.29316014},
//            {-0.07912563, 0.5320706, -0.21625209},
//            {0.40619087, 0.5689784, -0.052345287},
//            {-0.026628207, 0.48960224, -0.22401923},
//            {-0.00594549, 0.5292208, -0.22697341},
//            {0.03870894, 0.4061252, -0.040926915},
//            {0.23598538, 0.8778643, 0.20151441},
//            {-0.035166945, 0.5321401, -0.24074848},
//            {-0.026628207, 0.5529919, -0.22401923},
//            {0.35350367,  0.974617, 0.35350496},
//            {-0.035166945, 0.5104541, -0.24074848},
//            {0.25280985, 0.75365657, -0.01815654},
//            {-0.011005484, 0.61846685, -0.19786319},
//            {0.40574813, 0.40702876, -0.14425656},
//            {-0.0316882, 0.66746885, -0.123444594},
//            {-0.042187683, 0.62089956, -0.11146193},
//            {-0.01960747, 0.4078628, -0.045872655},
//            {-0.022833215, 0.37693256, -0.051150344},
//            {0.33168247, 0.5625143, 0.22630955},
//            {0.08949861, 0.40702876, -0.021708194},
//            {0.13820103, 0.4073763, -0.043283608},
//            {0.3479377, 0.7226568, 0.042520165},
//            {0.08994135, 0.4736158, -0.11361948},
//            {0.107714586, 0.56307036, -0.12258158},
//            { 0.328773, 0.6225678, 0.3871626},
//            {0.08994135, 0.5689784, -0.11361948},
//            {0.26836935, 1.0171549, 0.51565254},
//            {0.2789953, 0.6270856, -0.021575417},
//            {-0.01960747, 0.63473135, -0.045872655},
//            {-0.022833215, 0.66566163, -0.051150344},
//            {-0.01922797, 0.63215965, -0.025392609},
//            {-0.07286389, 0.40459603, 0.38414204},
//            {0.4574865, 0.5312365, 0.4619795},
//            {-0.034344696, 0.61610365, -0.05340746},
//            {-0.069258645, 0.39055577, 0.15457958},
//            {-0.3584372, 0.62089956, 0.2774267},
//            {0.28508627, 0.64000684, 0.13572599},
//            {-0.060277157, 0.3607376, 0.12839034},
//            {-0.08469162, 0.34475115, 0.14843889},
//            {-0.09971979, 0.39867413, 0.11859842},
//            {0.5624181, 0.48960224, 0.31888467},
//            {-0.1615276,  0.498082, 0.41760054},
//            {-0.10866333, 0.35302237, 0.09755409},
//            {-0.016508223, 0.7542126, 0.0918449},
//            {-0.096203096, 0.35135424, 0.08145551},
//            {-0.090023585, 0.2692395, 0.19905813},
//            {0.28905204, 0.59149843, 0.14146836},
//            {0.1007571, 0.75984263, 0.0197498},
//            {0.28190482, 0.61610365, 0.21937221},
//            {0.2944283, 0.6452268, 0.10024271},
//            { 0.314099, 0.6392492, 0.09815156},
//            {-0.04870242, 0.5212971, 0.24098083},
//            {0.2903803, 0.6753925, 0.09380328},
//            {0.37108716, 0.62305427, 0.13380079},
//            {0.24937539, 0.6480766, 0.13230711},
//            {0.34287772, 0.5303329, 0.46470135},
//            { 0.327255, 0.7716587, 0.3000311},
//            {0.030676205, 0.75490767, 0.45763123},
//            {0.3190325,  0.531167, 0.44256166},
//            {0.24509336, 0.6796324, 0.13449784},
//            {-0.26141185, 0.62305427, 0.13380079},
//            {-0.2896213, 0.5303329, 0.46470135},
//            {0.23978037, 0.63104755, 0.14681242},
//            {0.24338561, 0.63799816, 0.113752246},
//            { 0.335857, 0.7553942, 0.45202166},
//            {0.48968074, 0.53019387, 0.12925336},
//            {0.010562734, 0.8937117, -0.0011617483},
//            {0.20448694,  0.984765, 0.25651512},
//            {-0.13820103, 0.9827493, 0.043283608},
//            {-0.07552038, 0.5212971, 0.038702987},
//            {0.23079889, 1.0156257, 0.26175964},
//            {-0.25964084, 0.87779474, 0.20473413},
//            {0.2934163, 0.5212971,  0.335448},
//            {0.2817783, 0.8688285, -0.06598759},
//            {0.25280985, 0.9840004, 0.31377298},
//            {0.33908272, 0.8688285, 0.32841107},
//            {0.20157744, 0.9839309, 0.3169927},
//            {-0.1944302, 0.02606487, 0.03774038},
//            {-0.23480892, 0.19468708, 0.06562248},
//            {0.14756203, 0.98302734,  0.295882},
//            {0.38783577, 0.8688285, 0.2989357},
//            {0.23940088, 0.8688285, -0.111428745},
//            {0.39537513,  0.858055, 0.2816422},
//            {0.25280985, 0.9840004, -0.01815654},
//            {0.2817783, 0.8688285, 0.26594195},
//            {0.050093923, 0.8598622, 0.34377944},
//            {0.07008089, 0.83713365, 0.34497437},
//            {-0.19974321, 0.18266249, 0.24742028},
//            {0.09445107, 0.8495057, 0.32250273},
//            {-0.26141185, 0.4195399, 0.13380079},
//            {-0.14477903, 0.21546946, 0.17409705},
//            {0.09445107, 0.8881512, 0.32250273},
//            {0.22622591, 0.77335465, 0.29883614},
//            {0.5916396, 0.5312365, -0.07358878},
//            {0.20157744, 0.7537261, 0.3169927},
//            {-0.2253594, 0.83713365, 0.03389001},
//            {0.22377814, 0.75247496, 0.26255625},
//            {0.121819325, 1.0165293, 0.12822439},
//            {0.29088628, 0.9784399, 0.5159513},
//            {-0.1944302, 0.7211276, 0.3696699},
//            {-0.23480892, 0.8897499,  0.397552},
//            { 0.335857, 0.98226273, 0.45202166},
//            {0.09089011, 0.83713365, 0.13207476},
//            {0.33908272,  1.013193, 0.44674394},
//            {-0.26141185, 0.62305427, -0.19812873},
//            {0.33547747,  0.979691, 0.47250167},
//            {0.27754056, 0.9840004, 0.45696738},
//            {-0.2246004, 0.5114272, 0.09971163},
//            {0.18222298, 0.9834444, 0.51369417},
//            {0.042187683,  0.968431, 0.44339147},
//            {0.18525895, 1.0154172, 0.51545334},
//            {0.2267509, 0.9830968, 0.4761861},
//            {0.11954232, 0.9851125, 0.4555733},
//            {0.29088628, 0.9784399, 0.18402171},
//            {0.121819325, 1.0165293, 0.4601539},
//            {0.17804848, 0.9827493, 0.45461068},
//            {-0.20448694, 0.6372336, 0.07541439},
//            {0.28322673, 0.8688285, 0.46311805},
//            {-0.3479377, 0.66746885, 0.28940937},
//            {0.32681224, 0.84394526, 0.49905604},
//            {0.09445107, 0.8495057, -0.0094268},
//            {0.32681224, 0.8937117, 0.49905604},
//            {0.22744031, 0.8688285, 0.47956848},
//            {0.52364594, 0.5113577, 0.19281788},
//            {0.48942775, 0.47910678, -0.045573927},
//            {-0.20024918, 0.07423272, 0.02509388},
//            {-0.24009663,  0.575651, -0.03797275},
//            {0.2578825, 0.9642467, 0.48099908},
//            {-0.24009663,  0.575651, 0.2939568},
//            {0.1964795, 0.9639895, 0.48060074},
//            {0.2580596, 0.9099067, 0.48076674},
//            {-0.14541152, 0.37484738, 0.03455386},
//            {0.29702154, 0.63215965, -0.14057215},
//            {-0.39537513, 0.5105236, 0.38221687},
//            {-0.0066412394, 0.28490624, 0.2025766},
//            {0.19683369, 0.9099067, 0.47987053},
//            {0.33908272,  1.013193, 0.11481442},
//            {0.09089011, 0.83713365, 0.46400428},
//            {0.1318128, 0.8688285, 0.2618592},
//            {0.29506078, 0.8688285, -0.11428334},
//            {0.10885309, 0.87876785, 0.4710412},
//            {0.33168247, 0.48007986, -0.10561997},
//            {-0.170838, 0.027315984, 0.20051864},
//            {0.1740764, 0.8688285, 0.4593241},
//            {0.2712156, 0.75261396, 0.5145903},
//            {-0.25964084, 0.8598622, 0.53666365},
//            {-0.111762576, 0.05782921, 0.24137916},
//            {0.35059422, 0.7740219, 0.44448683},
//            {0.27754056, 0.75365657, 0.45696738},
//            {0.5272512, 0.66830283, -0.042719338},
//            {0.2578825, 0.77341026, 0.48099908},
//            {0.5129567, 0.6375811, -0.12364375},
//            {0.2267509, 0.7545602, 0.4761861},
//            {-0.35467383, 0.49863806, 0.027616538},
//            {0.010562734, 0.8937117, 0.33076778},
//            {0.20448694,  0.984765, -0.07541439},
//            {0.18222298, 0.7542126, 0.51369417},
//            {0.11954232, 0.7525445, 0.4555733},
//            {0.010562734, 0.19864893, -0.0011617483},
//            {-0.13820103, 0.2876865, 0.043283608},
//            {0.121819325, 0.7211276, 0.4601539},
//            {-0.19670719, 0.057481706, 0.04232101},
//            {-0.39537513, 0.5320706, 0.050287317},
//            {-0.18728295, 0.8588891, 0.2171151},
//            {0.10391958, 0.75608927, 0.44272763},
//            {0.0316882, 0.7226568, 0.123444594},
//            {0.17804848, 0.75490767, 0.45461068},
//            {0.1007571, 0.9778143, 0.35167933},
//            {-0.20739642, 0.8588891, 0.35878265},
//            {0.105247825, 1.0158342, 0.37464887},
//            {  0.15161, 0.98358333, 0.40731075},
//            { 0.322195, 0.86090475, 0.27092087},
//            {0.15309007, 0.8688285, 0.34755674},
//            {0.30315676, 0.7734659, 0.43784824},
//            {0.14059821, 0.9640312, 0.3772379},
//            {0.14307128, 0.9110188, 0.37750345},
//            {0.09980834, 0.8688285, 0.40956786},
//            {0.31030402, 0.5292208, 0.06100865},
//            {0.35059422, 0.9636351, 0.44448683},
//            {0.15217926, 0.9096287, 0.43818018},
//            {0.030676205, 0.9827493, 0.12570171},
//            {0.14059821, 0.77362573, 0.3772379},
//            {0.15114197, 0.77354926, 0.43871126},
//            {0.30948177, 0.5212971, -0.15743418},
//            {0.15932651, 0.7542126, 0.35267514},
//            {-0.16407025, 0.82802826, 0.059714116},
//            {0.14541152, 0.7223788, 0.29737568},
//            {0.31182203, 0.72536755, -0.03661183},
//            {0.011005484, 0.7716587, 0.19786319},
//            {  0.15161, 0.7540736, 0.40731075},
//            {0.14307128, 0.8266382, 0.37750345},
//            {-0.32681224, 0.5461803, 0.16480301},
//            {0.106006846, 0.75379556, 0.3711636},
//            {0.15217926, 0.82802826, 0.43818018},
//            {0.1007571, 0.75984263, 0.35167933},
//            {0.27754056, 0.9840004, 0.12503785},
//            {0.17018019, 0.77335465, 0.32223716},
//            {0.12372312, 0.8688285, 0.32094267},
//            {0.2934163, 0.5212971, 0.0035184533},
//            {0.17018019, 0.9643023, 0.32223716},
//            {-0.25964084, 0.87779474, 0.53666365},
//            {-0.30081654, 0.5625143, -0.10561997},
//            {0.08144058, 0.84790707, 0.4322718},
//            {-0.06343965, 0.9840004, 0.1841213},
//            {0.08144058, 0.8897499, 0.4322718},
//            {-0.08949861,  0.288034, 0.021708194},
//            {0.48031977, 0.48049688, -0.106250644},
//            {0.037064448, 0.8598622, 0.43814698},
//            {0.093293615, 0.96863955, 0.3838433},
//            {0.022833215,  1.013193, 0.38307986},
//            {0.30578798, 0.8688285, 0.40733397},
//            {-0.056608666, 0.53026336, 0.03876937},
//            {0.31118318, 0.96456647, 0.37640807},
//            {0.30967152, 0.91074073, 0.37743706},
//            {-0.3469257, 0.4073763, 0.29166648},
//            {0.24009663, 0.9231824, 0.03797275},
//            {0.1615276, 0.8920436, -0.085671015},
//            {0.31118318, 0.77309054, 0.37640807},
//            {-0.20157744, 0.4061947, 0.014936824},
//            {-0.16407025, 0.13296549, 0.059714116},
//            {0.29974127, 0.7542126, 0.40604937},
//            {0.30081654, 0.82761127, 0.4375495},
//            {-0.021188717, 0.8688285, 0.6121776},
//            {0.42801207, 0.4053606, 0.07541439},
//            {0.29196155, 0.7552552, 0.35118148},
//            {0.23079889, 0.72203124, -0.0701699},
//            {0.30967152, 0.8269162, 0.37743706},
//            {-0.03870894, 0.058593784, 0.040926915},
//            {-0.27918506, 0.8598622, 0.059747316},
//            {0.3479377, 0.7226568, 0.37444967},
//            {-0.17147048, 0.4795933, 0.3400618},
//            {0.18538547, 0.8688285, -0.07040225},
//            {0.28363785, 0.7734659, 0.32140738},
//            {0.25622535, 0.8688285, 0.3065303},
//            {0.037254192, 0.76303995, 0.47631887},
//            {0.2815253, 0.91004574, 0.3226687},
//            {0.077734135, 0.8598622, 0.4794722},
//            {0.39651364, 0.5303329, 0.13041511},
//            {0.077734135, 0.87779474, 0.4794722},
//            {-0.09247135, 0.057412185, 0.23533802},
//            {-0.032320693, 0.41863632, 0.07694126},
//            {-0.17501247, 0.5113577, 0.13004999},
//            {0.23687088, 0.8688285, 0.4037259},
//            {0.34886116, 0.4259345, 0.3424517},
//            {0.21973015, 0.8688285, 0.42779076},
//            {0.22630815, 0.9165098, 0.3842748},
//            {0.22453715, 0.9418796, 0.3969877},
//            {0.31118318, 0.77309054, 0.044478558},
//            {0.24009663, 0.9231824, 0.36990228},
//            {0.20853493, 0.91060174, 0.37531272},
//            {0.22630815, 0.8211472, 0.3842748},
//            {0.41276884, 0.5212971, 0.23606828},
//            {0.24009663, 0.8144746, 0.36990228},
//            {-0.29974127, 0.40668124, 0.25780967},
//            {0.20853493, 0.8270552, 0.37531272},
//            {-0.09089011, 0.5529919, 0.19985478},
//            {0.116000324, 0.9683615, 0.47280043},
//            {0.1912677,  0.968014, 0.52926165},
//            {-0.3584372, 0.62089956, -0.054502834},
//            {0.16242574, 0.9772583, 0.51419204},
//            {-0.34287772, 0.5529919, 0.058054473},
//            {0.26141185, 0.97058564, 0.53005826},
//            {0.15508877, 0.76206684, 0.2774267},
//            {0.30960828, 0.7576879, 0.2953177},
//            {-0.2660291, 0.9005233, 0.40714473},
//            {0.084312126, 0.87876785, 0.27633134},
//            {0.24066588, 0.83713365, -0.110598914},
//            {0.22622591, 0.9643023, 0.29883614},
//            {0.22592865, 0.91032374, 0.29983196},
//            {-0.111762576,  0.752892, 0.5733087},
//            {-0.0050663175, 0.96456647, 0.1214862},
//            {0.28363785, 0.7734659, -0.010522162},
//            {0.24566263, 0.48953274, 0.3654876},
//            {0.18538547, 0.8688285, 0.26152727},
//            {-0.056608666, 0.53026336, 0.3706989},
//            {0.1615276, 0.84561336, 0.24625853},
//            {0.31118318, 0.96456647, 0.044478558},
//            {-0.3469257, 0.4073763, -0.040263046},
//            {0.1615276, 0.8920436, 0.24625853},
//            {-0.19974321, 0.8777253, 0.5793498},
//            {0.15926325, 0.8688285, 0.21226895},
//            { 0.159706, 0.83713365,  0.213829},
//            {-0.2229559, 0.76901746, 0.4459805},
//            {0.35467383, 0.89148754, 0.30431297},
//            {0.34287772, 0.83713365, 0.27387506},
//            {0.35141644, 0.85798544, 0.2571458},
//            { 0.322195, 0.8767522, 0.27092087},
//            {-0.0347242, 0.21498291, 0.1752256},
//            {0.35141644, 0.8796715, 0.2571458},
//            {0.01922797, 0.28462824, 0.025392609},
//            {0.31346652, 0.87869835, 0.22129741},
//            {0.26822385, 0.8492972, 0.23952034},
//            {0.5087759, 0.5212971, 0.3429164},
//            {0.26822385, 0.8883597, 0.23952034},
//            {0.24066588, 0.83713365, 0.2213306},
//            {0.084312126, 0.87876785, -0.055598192},
//            {-0.1686875, 0.7546296, 0.53394186},
//            {0.38783577, 0.8688285, -0.032993797},
//            {0.23940088, 0.8688285, 0.2205008},
//            {0.29506078, 0.8688285, 0.2176462},
//            {0.1318128, 0.8688285, -0.07007033},
//            {0.12959905, 0.9687785, 0.29823866},
//            {0.23465714, 0.9662763, 0.2501089},
//            { 0.327255, 0.96599823, 0.3000311},
//            {-0.30960828, 0.41015655, 0.36854136},
//            {0.30372605, 0.9700991, 0.27669644},
//            {0.20157744, 0.9839309, -0.014936824},
//            {0.33908272, 0.8688285, -0.0035184533},
//            {0.14756203, 0.98302734, -0.036047544},
//            {0.042187683, 0.2733682, 0.11146193},
//            {0.075836636, 0.71813893, 0.1920544},
//            {0.15508877, 0.9755901, -0.054502834},
//            {0.04085943, 0.8588891, 0.07358878},
//            {-0.26822385, 0.5408284, 0.4243387},
//            {0.39942312, 0.8598622, -0.017691841},
//            {0.39537513,  0.858055, -0.050287317},
//            {0.55697864, 0.5212971, 0.3706325},
//            {-0.16315944, 0.8688285, 0.15033753},
//            {0.39537513,  0.879602, -0.050287317},
//            {-0.0050663175, 0.77309054, 0.45341572},
//            {0.075646885, 0.8589586, 0.031898417},
//            {0.050093923, 0.8598622, 0.011849892},
//            {0.20157744, 0.7537261, -0.014936824},
//            {0.5916396, 0.5312365, 0.25834075},
//            {0.22377814, 0.75247496, -0.069373265},
//            {-0.24066588, 0.5529919, 0.44252843},
//            {0.20448694,  0.752892, -0.07541439},
//            {-0.1944302, 0.7211276, 0.03774038},
//            {-0.23480892, 0.8897499, 0.06562248},
//            { 0.335857, 0.98226273, 0.1200921},
//            {0.33547747,  0.979691, 0.14057215},
//            {-0.2246004, 0.5114272, 0.43164113},
//            {0.18222298, 0.9834444, 0.18176462},
//            {0.042187683,  0.968431, 0.11146193},
//            {0.18525895, 1.0154172, 0.18352382},
//            {0.53269064, 0.5212971, 0.2542912},
//            {0.2267509, 0.9830968, 0.14425656},
//            {0.11954232, 0.9851125, 0.12364375},
//            {0.32681224, 0.84394526, 0.1671265},
//            {0.32681224, 0.8937117, 0.1671265},
//            {0.24509336, 0.3629618, 0.13449785},
//            {0.1964795, 0.9639895, 0.14867124},
//            {0.2580596, 0.9099067, 0.14883721},
//            {0.22592865, 0.91032374, -0.03209759},
//            {-0.060024157, 0.1737657,  0.191364},
//            {-0.30967152, 0.4793848, -0.045507535},
//            {0.10885309, 0.87876785, 0.13911167},
//            {0.5481869, 0.5312365, 0.3875277},
//            {0.2712156, 0.75261396, 0.18266082},
//            {-0.25964084, 0.8598622, 0.20473413},
//            { 0.335857, 0.7553942, 0.1200921},
//            {0.27754056, 0.75365657, 0.12503785},
//            {0.2578825, 0.77341026, 0.14906956},
//            {0.1964795, 0.77366745, 0.14867124},
//            {-0.0916491, 0.18363558, 0.2656764},
//            {0.2267509, 0.7545602, 0.14425656},
//            {0.19683369, 0.82775027, 0.14794098},
//            {0.18222298, 0.7542126, 0.18176462},
//            {0.18525895, 0.72223973, 0.18352382},
//            {0.11954232, 0.7525445, 0.12364375},
//            {0.121819325, 0.7211276, 0.12822439},
//            {-0.39537513, 0.5320706, 0.38221687},
//            {-0.18728295, 0.8588891, 0.5490446},
//            {0.10391958, 0.75608927, 0.11079808},
//            {0.0316882, 0.7226568, 0.45537412},
//            {0.17804848, 0.75490767, 0.12268115},
//            {-0.076152876, 0.1194118, 0.12799202},
//            {0.1007571, 0.9778143, 0.0197498},
//            {-0.20739642, 0.8588891, 0.026853103},
//            {0.105247825, 1.0158342, 0.042719338},
//            {  0.15161, 0.98358333, 0.0753812},
//            {0.15932651, 0.9834444, 0.020745596},
//            {0.15309007, 0.8688285, 0.015627235},
//            {0.14345077, 0.8688285, 0.07630064},
//            {0.14059821, 0.9640312, 0.04530838},
//            {0.30315676, 0.7734659, 0.10591872},
//            {0.14307128, 0.9110188, 0.045573927},
//            {-0.021188717, 0.1737657, 0.2802481},
//            {0.09980834, 0.8688285, 0.07763832},
//            {-0.13086404, 0.8688285, 0.56829655},
//            {0.15114197, 0.9641077, 0.10678172},
//            {0.15217926, 0.9096287, 0.106250644},
//            {0.030676205, 0.9827493, 0.45763123},
//            {-0.25964084, 0.18273202, 0.20473413},
//            {0.14059821, 0.77362573, 0.04530838},
//            {0.15114197, 0.77354926, 0.10678172},
//            {0.15932651, 0.7542126, 0.020745596},
//            {-0.16407025, 0.82802826, 0.39164364},
//            {0.5285794, 0.40855792, 0.22113144},
//            {  0.15161, 0.7540736, 0.0753812},
//            {0.106006846, 0.75379556, 0.039234072},
//            {-0.22179843, 0.19308847, 0.17539157},
//            {0.15217926, 0.82802826, 0.106250644},
//            {0.093293615, 0.76901746, 0.051913787},
//            {0.17018019, 0.77335465, -0.009692338},
//            {0.47007325, 0.41286728, -0.18226251},
//            {0.17147048, 0.8271248, -0.008132277},
//            {0.34287772, 0.5303329, 0.13277182},
//            {0.12372312, 0.8688285, -0.01098687},
//            {0.17018019, 0.9643023, -0.009692338},
//            {-0.14606933, 0.078291886, 0.17565711},
//            {-0.012523481, 0.9700991, 0.22119783},
//            {0.17147048, 0.91053224, -0.008132277},
//            {0.19705509, 0.8688285, -0.025399245},
//            {-0.21232992, 0.28650486, 0.055166684},
//            {-0.15382376, 0.9772583, -0.016297754},
//            {0.08144058, 0.84790707, 0.10034229},
//            {-0.1547219, 0.19698079, 0.2516358},
//            {0.04870242, 0.8688285, 0.090948686},
//            {-0.06343965, 0.9840004, 0.5160509},
//            {0.08144058, 0.8897499, 0.10034229},
//            {0.050220422, 0.83713365, 0.09074953},
//            {0.037064448, 0.8598622, 0.106217444},
//            {-0.17804848, 0.63521785, -0.12268115},
//            {0.04085943, 0.87876785, 0.07358878},
//            {0.093293615, 0.96863955, 0.051913787},
//            {0.022833215,  1.013193, 0.051150344},
//            {0.35350367,  0.974617, 0.021575417},
//            {0.29974127, 0.9834444, 0.07411986},
//            {0.29196155, 0.9824017, 0.019251922},
//            {0.30578798, 0.8688285, 0.07540443},
//            {0.30075958, 0.8688285, 0.015159213},
//            {0.30081654, 0.91004574, 0.10561997},
//            {0.29974127, 0.7542126, 0.07411986},
//            {0.30081654, 0.82761127, 0.10561997},
//            {-0.021188717, 0.8688285, 0.2802481},
//            {0.29196155, 0.7552552, 0.019251922},
//            {0.53920543, 0.42148608, 0.28001574},
//            {0.42801207, 0.4053606, 0.40734392},
//            {0.30967152, 0.8269162, 0.045507535},
//            {-0.27918506, 0.8598622, 0.39167684},
//            {0.3584372, 0.76922596, 0.054502834},
//            {-0.11941582, 0.13268751, 0.018023778},
//            {0.2815253, 0.82761127, -0.009260837},
//            {0.25622535, 0.8688285, -0.025399245},
//            {0.022833215, 0.1737657, 0.16948321},
//            {0.28363785, 0.9641911, -0.010522162},
//            {0.2815253, 0.91004574, -0.009260837},
//            {0.3584372,  0.968431, 0.054502834},
//            {0.26836935, 0.7205021,  0.183723},
//            {0.077734135, 0.87779474, 0.14754269},
//            {0.39651364, 0.5303329, 0.46234465},
//            {-0.032320693, 0.41863632, 0.4088708},
//            {-0.17501247, 0.5113577, 0.4619795},
//            {0.23687088, 0.8688285, 0.07179636},
//            {0.21973015, 0.8688285, 0.09586125},
//            {0.22630815, 0.9165098, 0.052345287},
//            {0.22453715, 0.9418796, 0.06505818},
//            {0.22630815, 0.8211472, 0.052345287},
//            {0.20853493, 0.8270552, 0.04338319},
//            {0.57589036, 0.53026336, 0.3706989},
//            {0.1912677,  0.968014, 0.1973321},
//            {-0.04085943, 0.5312365, -0.07358878},
//            {0.16242574, 0.9772583, 0.18226251},
//            {-0.34287772, 0.5529919,  0.389984},
//            {0.26141185, 0.97058564, 0.19812873},
//            {0.35059422, 0.9636351, 0.1125573},
//            {0.31030402, 0.5292208, 0.3929382},
//            {0.26836935, 1.0171549,  0.183723},
//            {0.037254192, 0.06797715, 0.14438935},
//            {0.14756203, 0.7546296, -0.036047544},
//            {-0.22179843, 0.8881512, 0.17539157},
//            {0.15508877, 0.76206684, -0.054502834},
//            {0.30960828, 0.7576879, -0.03661183},
//            { 0.327255, 0.7716587, -0.031898428},
//            {-0.09445107, 0.50197434, 0.3413563},
//            {-0.2661556, 0.16479939, 0.15411487},
//            {0.076595634, 0.9005233, -0.037574425},
//            {0.056608666, 0.8598622, -0.03876937},
//            {0.056608666, 0.87779474, -0.03876937},
//            {0.3642752, 0.50176585, 0.4243387},
//            {0.22622591, 0.9643023, -0.033093374},
//            {   0.5029, 0.6212471, 0.36562037},
//            {0.1615276, 0.84561336, -0.085671015},
//            {-0.09445107, 0.50197434, 0.0094268},
//            {0.13339405, 0.9005233, -0.07030268},
//            {0.12896655, 0.8588891, -0.051150344},
//            {0.15926325, 0.8688285, -0.11966059},
//            { 0.159706, 0.9005233, -0.118100524},
//            {0.22004642, 0.6912399, 0.08450925},
//            {0.37627366, 0.5212971, 0.35732877},
//            {0.14281829, 0.85993165, -0.12925336},
//            { 0.159706, 0.83713365, -0.118100524},
//            {0.49910495, 0.48960224, 0.4022322},
//            {0.14281829, 0.8777253, -0.12925336},
//            {0.3412332, 0.8688285, -0.05898388},
//            {0.35467383, 0.89148754, -0.027616538},
//            {-0.13820103, 0.059844896, 0.043283608},
//            {0.035166945, 0.1846087, 0.24074848},
//            {0.34287772, 0.83713365, -0.058054473},
//            {0.35141644, 0.85798544, -0.07478373},
//            { 0.322195, 0.8767522, -0.06100865},
//            {-0.24616863, 0.83713365, 0.15291992},
//            {0.35141644, 0.8796715, -0.07478373},
//            {0.31346652, 0.87869835, -0.11063211},
//            {-0.03870894, 0.28893763, 0.040926915},
//            {0.24066588, 0.9005233, -0.110598914},
//            { 0.320677,  0.664758, 0.03661183},
//            {0.26822385, 0.8492972, -0.09240918},
//            {0.4062731, 0.42582327, 0.3650229},
//            {-0.071156144, 0.3629618, 0.36339647},
//            {0.26822385, 0.8883597, -0.09240918},
//            {0.2246004, 0.87869835, -0.09971163},
//            {-0.35059422, 0.42649055, 0.21937221},
//            {-0.1646395, 0.28852054, 0.090583555},
//            {0.5006862, 0.5212971, 0.40199986},
//            {-0.30315676, 0.6166597, -0.10591872},
//            {0.23598538, 0.8597927, -0.13041511},
//            {0.23598538, 0.8778643, -0.13041511},
//            {-0.14059821, 0.42609435, 0.28662115},
//            {0.29430178, 0.9005233, -0.112955615},
//            {0.23465714, 0.9662763, -0.08182062},
//            { 0.327255, 0.96599823, -0.031898428},
//            {-0.20739642, 0.1638263, 0.026853103},
//            {0.30372605, 0.9700991, -0.05523307},
//            {0.027197465, 0.79862714, 0.30743313},
//            {-0.08545062, 0.72203124, 0.5680642},
//            {0.07058689, 0.83706415, 0.53145236},
//            {0.083173625, 0.8598622, 0.1836566},
//            {0.06687412, 0.74204904, 0.29827186},
//            {0.0021504953, 0.7508763, 0.26411632},
//            {0.025869206, 0.71473306, 0.25976804},
//            {0.40619087, 0.4736158, -0.052345287},
//            {0.032320693, 0.7661677, 0.25498825},
//            {0.071156144, 0.7104932, 0.3004626},
//            {0.07646914, 0.7590781, 0.3127772},
//            {0.07912563,  0.858055, 0.5481816},
//            {-0.14123702, 0.18370508, 0.29601476},
//            {0.39784187, 0.6187449, 0.41375014},
//            {-0.076595634, 0.5529919, 0.36950395},
//            {0.07286389, 0.7521274,  0.279717},
//            {-0.2712156, 0.40508258, -0.18266082},
//            {-0.08545062, 1.0156257, 0.5680642},
//            {-0.111762576,  0.984765, 0.5733087},
//            {0.43566534, 0.4802189, 0.18398854},
//            {-0.04802565, 0.8883597, 0.5903035},
//            {-0.04802565, 0.8492972, 0.5903035},
//            {-0.021947715, 0.9005233, 0.61084986},
//            {-0.07558363, 0.83713365, 0.6084932},
//            {-0.15654351, 0.83713365, 0.28406528},
//            {-0.2253594, 0.9005233, 0.03389001},
//            {-0.08026413, 0.8597927, 0.6283094},
//            {0.3381972, 0.48960224, 0.44488513},
//            {-0.026628207, 0.8597927, 0.6306661},
//            {-0.15654351, 0.9005233, 0.6159948},
//            {-0.15654351, 0.83713365, 0.6159948},
//            {-0.2253594, 0.9005233, 0.36581954},
//            {-0.07558363, 0.83713365, 0.27656367},
//            {-0.17343123, 0.85993165, 0.6271477},
//            {-0.08880919, 0.1737657, 0.018325826},
//            {-0.23193738, 0.8588891, 0.5534925},
//            {0.024983712, 0.8688285, 0.55687815},
//            {0.026628207, 0.83713365, 0.55594873},
//            {-0.27918506, 0.87779474, 0.39167684},
//            {-0.13820103, 0.75490767, 0.37521312},
//            {0.035166945, 0.8796715,  0.572678},
//            {0.026628207, 0.9005233, 0.55594873},
//            {-0.08159237, 0.77138066, 0.24778539},
//            {0.00594549, 0.86090475, 0.5589029},
//            {-0.19670719, 0.9851125, 0.37425053},
//            {-0.119770005, 0.77366745, 0.34922305},
//            {-0.08159237, 0.77138066, 0.5797149},
//            {0.026628207, 0.9005233, 0.22401923},
//            {-0.004427493, 1.0122894, 0.5345061},
//            {-0.1686875, 0.98302734, 0.53394186},
//            {-0.1646395, 0.7540736, 0.4225131},
//            {-0.16116074, 0.9755901, 0.55239713},
//            {-0.16510755, 0.77354926, 0.059183042},
//            {-0.170838, 1.0152782, 0.5324481},
//            {-0.18665047, 0.9687785, 0.53158516},
//            {0.060277157,  0.708269, 0.20353918},
//            {-0.11467207, 0.9839309, 0.5128311},
//            {-0.21024267, 0.9838614, 0.45866022},
//            {-0.024287963, 0.7552552, 0.47864237},
//            {-0.21549243, 0.9778143, 0.4781445},
//            {0.29702154, 0.41043457, 0.19135737},
//            {-0.21100168, 1.0158342, 0.45517495},
//            {-0.2229559, 0.96863955, 0.4459805},
//            {0.037254192,  0.974617, 0.47631887},
//            {0.54160887, 0.48960224, 0.19985478},
//            {-0.08994135, 0.9165098, 0.44554898},
//            {0.038424317, 0.89148754, 0.52551085},
//            {0.07158624, 0.8688285, 0.5308881},
//            {0.07058689, 0.90059286, 0.53145236},
//            {0.39183316, 0.48960224, 0.44252843},
//            {0.083173625, 0.87779474, 0.51558614},
//            {0.29664204, 0.4078628, -0.1200921},
//            {0.07912563,  0.879602, 0.5481816},
//            {0.083173625, 0.8598622, 0.51558614},
//            {-0.08545062, 0.72203124, 0.23613468},
//            {0.07058689, 0.83706415, 0.19952284},
//            {-0.0154899005, 0.8688285, 0.48273507},
//            {-0.054837663, 0.7670713, 0.29976556},
//            {-0.0050663175, 0.96456647, 0.45341572},
//            {-0.111762576,  0.752892, 0.24137916},
//            {-0.03261165, 0.9641911, 0.5084164},
//            {-0.21549243, 0.75984263, 0.4781445},
//            {-0.0065779896, 0.91074073, 0.45238677},
//            {-0.0347242, 0.91004574, 0.5071551},
//            {-0.090023585, 0.9643023, 0.5309877},
//            {-0.21973015, 0.5212971, -0.09586125},
//            {-0.22179843, 0.15444294, 0.17539157},
//            {-0.0050663175, 0.77309054, 0.1214862},
//            {-0.16315944, 0.8688285, 0.48226705},
//            {-0.14606933, 0.9643023, 0.50758666},
//            {-0.1756513, 0.9640312, 0.4525859},
//            {-0.11954232, 0.40501308, -0.12364375},
//            {-0.17317823, 0.9110188, 0.45232034},
//            {-0.08994135, 0.8211472, 0.44554898},
//            {-0.076152876, 0.8144746, 0.45992154},
//            {0.042187683, 0.76922596, 0.44339147},
//            {0.44711354, 0.5212971, 0.07040225},
//            {-0.03261165, 0.7734659, 0.5084164},
//            {-0.0065779896, 0.8269162, 0.45238677},
//            {-0.090320855, 0.8273333, 0.52999187},
//            {0.32289076, 0.63243765, 0.36854136},
//            {0.30524403, 0.61846685, 0.36382794},
//            {-0.004427493, 0.72536755, 0.5345061},
//            {-0.14606933, 0.77335465, 0.50758666},
//            {-0.15692301, 0.7542126, 0.47714868},
//            {-0.14477903, 0.8271248, 0.50602657},
//            {-0.17018019, 0.42582327, 0.34162188},
//            {-0.015432976, 0.82761127, 0.06034479},
//            {-0.17317823, 0.8266382, 0.45232034},
//            {-0.170838, 0.7223788, 0.5324481},
//            {0.48968074, 0.53019387, 0.46118286},
//            {-0.13820103, 0.9827493, 0.37521312},
//            {-0.1646395, 0.98358333, 0.4225131},
//            {-0.21232992, 0.9815677, 0.3870962},
//            {0.01960747, 0.98226273, 0.3778022},
//            {0.01922797,  0.979691, 0.35732213},
//            {-0.07646914, 0.6310475, 0.35108188},
//            {0.0067677395, 0.8688285, 0.3404601},
//            {0.5559034, 0.5529919, 0.36950395},
//            {-0.23465714, 0.6187449, 0.41375014},
//            {-0.033022773, 0.8688285, 0.3667058},
//            {-0.23480892, 0.84790707,  0.397552},
//            {-0.05836701, 0.9642467, 0.3488247},
//            {-0.015432976, 0.91004574, 0.3922743},
//            {-0.045033928,  0.985043, 0.31523347},
//            {-0.14217313, 0.8688285, 0.3704997},
//            {-0.16407025, 0.9096287, 0.39164364},
//            {-0.119770005, 0.9639895, 0.34922305},
//            {-0.11941582, 0.9099067, 0.34995332},
//            {-0.05836701, 0.77341026, 0.3488247},
//            {0.01960747, 0.7553942, 0.3778022},
//            {0.022833215, 0.72446394, 0.38307986},
//            {-0.21100168, 0.7218227, 0.45517495},
//            {0.01922797,  0.757966, 0.35732213},
//            {-0.08949861, 0.7545602, 0.35363773},
//            {-0.05818991, 0.82775027, 0.3490571},
//            {-0.1615276,  0.498082, 0.085671015},
//            {-0.016508223, 0.7542126, 0.42377445},
//            {-0.17317823, 0.8266382, 0.12039083},
//            {-0.17018019, 0.42582327, 0.009692338},
//            {-0.015432976, 0.82761127, 0.3922743},
//            {-0.16510755, 0.77354926, 0.39111257},
//            {-0.16116074, 0.9755901, 0.22046758},
//            {-0.13402654, 0.7542126, 0.31612965},
//            {-0.21232992, 0.75608927, 0.3870962},
//            {-0.21644118, 0.8688285, 0.42025596},
//            {-0.2660291, 0.83713365, 0.40714473},
//            {0.27406183, 0.62089956, -0.054502834},
//            {-0.2253594, 0.83713365, 0.36581954},
//            {0.55476487, 0.5123308, -0.14754269},
//            {-0.20739642, 0.87876785, 0.35878265},
//            {-0.23851536, 0.87779474, 0.3503516},
//            {-0.2753901, 0.87876785, 0.4243055},
//            {-0.2753901, 0.8588891, 0.4243055},
//            {-0.39942312, 0.53026336, 0.017691841},
//            {0.58379656, 0.5212971, -0.090948686},
//            {-0.09651935, 0.8688285, 0.40203303},
//            {-0.18665047, 0.76887846, 0.53158516},
//            {-0.107714586, 0.8270552, 0.4545111},
//            {-0.22179843, 0.8881512, 0.50732106},
//            {-0.19252639, 0.8688285, 0.50888115},
//            {-0.24072912, 0.8688285, 0.5365973},
//            {-0.16315944, 0.1737657, 0.15033753},
//            {-0.22179843, 0.8495057, 0.50732106},
//            {0.49910495, 0.5529919, 0.07030268},
//            {-0.23965387, 0.9005233, 0.53546876},
//            {-0.24616863, 0.83713365, 0.48484945},
//            {-0.24060263, 0.87869835, 0.46599588},
//            {-0.24060263, 0.8589586, 0.46599588},
//            {-0.23965387, 0.83713365, 0.53546876},
//            {-0.012523481, 0.9700991, 0.55312735},
//            {-0.1007571, 0.63028294, 0.3121797},
//            {-0.13099055, 0.72223973, 0.31437045},
//            {0.39240238,  0.575651, 0.2939568},
//            {-0.15382376, 0.7603987, 0.31563178},
//            {-0.047880176, 0.7205021, 0.31417128},
//            {-0.18443671, 0.8688285, 0.5679646},
//            {-0.18285546, 0.83713365, 0.56819695},
//            {0.47741026, 0.6280587, 0.38643235},
//            {-0.18285546, 0.9005233, 0.56819695},
//            {-0.19974321, 0.85993165, 0.5793498},
//            {-0.0916491, 0.87869835, 0.59760594},
//            {-0.13099055, 1.0154172, 0.31437045},
//            {-0.17343123, 0.18266249, 0.2952181},
//            {-0.124981806,  0.968014, 0.3005622},
//            {-0.105247825, 0.3742913, 0.2892102},
//            {-0.09247135, 0.2901192, 0.23533802},
//            {-0.034471195, 0.1737657, 0.23195235},
//            {0.00594549, 0.18168941, 0.22697341},
//            {-0.2578825, 0.42587888, -0.14906956},
//            {-0.047880176, 0.025439294, -0.017758248},
//            {-0.15309007, 0.5212971, 0.3163023},
//            {-0.04802565, 0.15423442, 0.25837395},
//            {-0.07558363, 0.20546056, 0.27656367},
//            {-0.17343123, 0.85993165, 0.2952181},
//            {-0.07558363, 0.14207084, 0.27656367},
//            {-0.1964795, 0.6164582, 0.18325828},
//            {-0.021947715, 0.14207084, 0.27892038},
//            {-0.0027829956, 0.18363558, 0.27659687},
//            {-0.026628207, 0.18280151, 0.29873657},
//            {-0.0027829956, 0.16389583, 0.27659687},
//            {-0.026628207, 0.16472986, 0.29873657},
//            {-0.15654351, 0.20546056, 0.28406528},
//            {-0.15654351, 0.14207084, 0.28406528},
//            {-0.2253594, 0.20546056, 0.03389001},
//            {-0.23193738, 0.18370508, 0.22156295},
//            {-0.19683369, 0.4802189, 0.18398854},
//            {-0.23193738, 0.1638263, 0.22156295},
//            {0.024983712, 0.1737657, 0.22494864},
//            { 0.480889, 0.63605195, 0.25654835},
//            {0.026628207, 0.14207084, 0.22401923},
//            {-0.27918506, 0.18273202, 0.059747316},
//            {0.026628207, 0.20546056, 0.22401923},
//            {0.00594549,  0.165842, 0.22697341},
//            {-0.19670719, 0.29004967, 0.04232101},
//            {-0.119770005, 0.07860463, 0.01729353},
//            {-0.31118318, 0.61703503, -0.044478558},
//            {0.035166945, 0.16292271, 0.24074848},
//            {-0.08545062, 0.026968436, 0.23613468},
//            {0.07058689, 0.14200132, 0.19952284},
//            {0.34053746, 0.63487035, 0.31267762},
//            {-0.08159237, 0.07631788, 0.24778539},
//            {-0.015432976, 0.91004574, 0.06034479},
//            {-0.004427493, 0.31722665, 0.2025766},
//            {-0.024287963, 0.28733897, 0.14671284},
//            {-0.16116074, 0.28052732, 0.22046758},
//            {0.40872085, 0.40494356, 0.40130278},
//            {-0.15692301, 0.28838158, 0.14521916},
//            {-0.17501247, 0.5312365, 0.4619795},
//            {-0.11467207, 0.2888681, 0.18090159},
//            {-0.21549243, 0.2827515, 0.14621496},
//            {-0.2229559, 0.2735767, 0.11405098},
//            {0.030676205, 0.2876865, 0.12570171},
//            {0.59543455, 0.5123308, 0.22571209},
//            {0.39183316, 0.5529919, 0.44252843},
//            {-0.35141644, 0.5104541, 0.07478373},
//            {0.0316882, 0.3199374, 0.123444594},
//            {-0.2253594, 0.1737657, 0.03229674},
//            {0.34287772, 0.5122613, 0.13277182},
//            {-0.08994135,  0.221447, 0.11361948},
//            {-0.07286389, 0.40459603, 0.052212514},
//            {0.4574865, 0.5312365, 0.13004999},
//            {-0.076152876, 0.2281196, 0.12799202},
//            {-0.091712356, 0.24681678, 0.10090658},
//            {0.038424317, 0.19642474, 0.1935813},
//            {0.07158624, 0.1737657, 0.19895856},
//            {0.038424317, 0.15110666, 0.1935813},
//            {-0.06343965, 0.058593784, 0.1841213},
//            {-0.2712156, 0.6375116, 0.1492687},
//            {0.07058689, 0.20553008, 0.19952284},
//            {0.083173625, 0.18273202, 0.1836566},
//            {-0.14217313, 0.8688285, 0.0385702},
//            {0.07912563, 0.18453918, 0.21625209},
//            {0.083173625, 0.16479939, 0.1836566},
//            {-0.03261165, 0.26912832, 0.17648691},
//            {-0.21549243, 0.06477987, 0.14621496},
//            {-0.11919444, 0.1737657,  0.191364},
//            {-0.1756513, 0.26896843, 0.12065638},
//            {0.35495847,  0.636469, 0.20689169},
//            {-0.091712356, 0.10071462, 0.10090658},
//            {-0.08159237, 0.2712135, 0.24778539},
//            {0.030676205, 0.059844896, 0.12570171},
//            {-0.18728295, 0.1638263, 0.2171151},
//            {0.30568677, 0.5461803, -0.1671265},
//            {0.0316882, 0.027593974, 0.123444594},
//            {0.35495847,  0.636469, -0.12503785},
//            {0.042187683, 0.074163206, 0.11146193},
//            {-0.0050663175, 0.07802773, 0.1214862},
//            {-0.22179843, 0.8495057, 0.17539157},
//            {-0.03261165, 0.078403085, 0.17648691},
//            {-0.0065779896, 0.13185342, 0.12045723},
//            {-0.2896213, 0.5122613, 0.46470135},
//            {-0.090023585, 0.078291886, 0.19905813},
//            {-0.0347242, 0.13254847, 0.1752256},
//            {0.011005484, 0.07659591, 0.19786319},
//            {-0.16242574, 0.62972695, -0.18226251},
//            {-0.0066412394, 0.06262515, 0.2025766},
//            {-0.090320855, 0.13227044, 0.19806236},
//            {-0.004427493, 0.030304754, 0.2025766},
//            {-0.15692301, 0.059149843, 0.14521916},
//            {0.5380479, 0.54061985, 0.3413563},
//            {-0.1686875, 0.059566863, 0.20201232},
//            {-0.13402654, 0.7542126, -0.015799856},
//            {-0.17317823, 0.13157539, 0.12039083},
//            {0.23978037, 0.4115467, 0.14681242},
//            {-0.016508223, 0.28838158, 0.0918449},
//            {-0.20024918, 0.27329868, 0.02509388},
//            {0.01960747, 0.28719994, 0.045872655},
//            {0.55476487, 0.5123308, 0.18438683},
//            {0.28108257, 0.5104541, 0.40671325},
//            {0.022833215, 0.31813025, 0.051150344},
//            {-0.25622535, 0.5212971, 0.35732877},
//            {0.0067677395, 0.1737657, 0.00853058},
//            {-0.033022773, 0.1737657, 0.03477626},
//            {-0.23480892, 0.15284432, 0.06562248},
//            {-0.010461533, 0.1737657, 0.09056034},
//            {0.3509737, 0.5625143, 0.34119037},
//            {-0.01309273, 0.26912832, 0.060046043},
//            {-0.3584372, 0.4216946, -0.054502834},
//            {-0.045033928, 0.2899802, -0.016696054},
//            {-0.17279872, 0.1737657, 0.08966412},
//            {-0.14217313, 0.1737657, 0.0385702},
//            {0.083173625, 0.87779474, 0.1836566},
//            {-0.119770005, 0.26892674, 0.01729353},
//            {-0.11941582, 0.21484388, 0.018023778},
//            {-0.13402654, 0.28838158, -0.015799856},
//            {-0.26822385, 0.50176585, 0.09240918},
//            {-0.05836701, 0.07834748, 0.016895207},
//            {0.01960747, 0.060331438, 0.045872655},
//            {0.022833215, 0.029401146, 0.051150344},
//            {-0.21100168, 0.026759924, 0.123245426},
//            {0.01922797, 0.06290318, 0.025392609},
//            {0.034344696, 0.078959145, 0.05340746},
//            {-0.01309273, 0.078403085, 0.060046043},
//            {0.36412966, 0.3729707, -0.183723},
//            {-0.08949861, 0.059497386, 0.021708194},
//            {-0.06887915, 0.5212971, -0.014040614},
//            {-0.05818991, 0.13268751, 0.017127557},
//            {-0.016508223, 0.059149843, 0.0918449},
//            {-0.015432976, 0.13254847, 0.06034479},
//            {0.5264922, 0.40626422, -0.039234072},
//            {0.40872085, 0.40494356, 0.069373265},
//            {-0.16510755, 0.07848648, 0.059183042},
//            {-0.13402654, 0.059149843, -0.015799856},
//            {-0.1686875, 0.7546296, 0.20201232},
//            {0.36412966, 0.6696235, -0.183723},
//            {-0.1964795, 0.42613602, 0.18325828},
//            {-0.21232992, 0.061026532, 0.055166684},
//            {-0.21644118, 0.1737657, 0.08832644},
//            {-0.2660291, 0.20546056, 0.07521523},
//            {0.5624181, 0.5529919, 0.31888467},
//            {-0.17147048, 0.56300086, 0.008132277},
//            {-0.2660291, 0.14207084, 0.07521523},
//            {0.2896213, 0.48960224, 0.058054473},
//            {-0.2253594, 0.14207084, 0.03389001},
//            {-0.20739642, 0.18370508, 0.026853103},
//            {-0.29506078, 0.5212971, 0.4462129},
//            {-0.23851536, 0.18273202, 0.01842208},
//            {-0.23851536, 0.16479939, 0.01842208},
//            {-0.2753901, 0.1638263, 0.09237599},
//            {-0.27918506, 0.16479939, 0.059747316},
//            {-0.03870894, 0.75365657, 0.040926915},
//            {-0.09651935, 0.1737657, 0.07010352},
//            {-0.21024267, 0.05873282, 0.1267307},
//            {-0.18665047, 0.07381565, 0.1996556},
//            {-0.107714586, 0.13199241, 0.12258158},
//            {-0.19252639, 0.1737657, 0.17695165},
//            {-0.116000324, 0.42176414, 0.19105865},
//            {-0.24737035, 0.1737657, 0.15192415},
//            {-0.24072912, 0.1737657, 0.20466773},
//            {-0.23965387, 0.20546056, 0.20353918},
//            {-0.24616863, 0.14207084, 0.15291992},
//            {-0.24060263, 0.18363558, 0.13406634},
//            {0.40657037, 0.47980183, 0.3640271},
//            {-0.17317823, 0.9110188, 0.12039083},
//            {-0.14756203, 0.40709826, 0.36797708},
//            {-0.2661556, 0.18273202, 0.15411487},
//            {-0.24060263, 0.16389583, 0.13406634},
//            {0.32131582, 0.61703503, -0.044478558},
//            {-0.23965387, 0.14207084, 0.20353918},
//            {-0.25964084, 0.16479939, 0.20473413},
//            {-0.012523481, 0.2750363, 0.22119783},
//            {0.3642752, 0.5408284, 0.09240918},
//            {0.011005484, 0.27093548, 0.19786319},
//            {-0.124981806, 0.07458023, -0.031367324},
//            {-0.18443671, 0.1737657, 0.2360351},
//            {-0.1964795, 0.6164582, -0.14867124},
//            {0.1276383, 0.7040291, 0.24101402},
//            {-0.18285546, 0.14207084, 0.23626743},
//            {-0.18285546, 0.20546056, 0.23626743},
//            {-0.19974321, 0.1648689, 0.24742028},
//            {-0.18728295, 0.18370508, 0.2171151},
//            {-0.0916491, 0.16389583, 0.2656764},
//            {-0.20448694, 0.4053606, 0.07541439},
//            {-0.13099055, 0.32035443, -0.017559057},
//            {-0.15382376, 0.2821955, -0.016297754},
//            {-0.21232992, 0.9815677, 0.055166684},
//            {-0.047880176, 0.32209212, -0.017758248},
//            {-0.08545062, 1.0156257, 0.23613468},
//            {0.35495847, 0.4061252, -0.12503785},
//            {-0.111762576,  0.984765, 0.24137916},
//            {-0.034471195, 0.8688285, 0.23195235},
//            {0.00594549, 0.8767522, 0.22697341},
//            {-0.047880176, 0.7205021, -0.017758248},
//            {-0.04802565, 0.8492972, 0.25837395},
//            {-0.021947715, 0.9005233, 0.27892038},
//            {-0.08026413, 0.8597927, 0.29637986},
//            {-0.026628207, 0.8597927, 0.29873657},
//            {-0.15654351, 0.9005233, 0.28406528},
//            {-0.17343123, 0.8777253, 0.2952181},
//            {-0.23193738, 0.8588891, 0.22156295},
//            {0.024983712, 0.8688285, 0.22494864},
//            {0.026628207, 0.83713365, 0.22401923},
//            {-0.27918506, 0.87779474, 0.059747316},
//            {-0.13820103, 0.75490767, 0.043283608},
//            {0.035166945, 0.8796715, 0.24074848},
//            {0.00594549, 0.86090475, 0.22697341},
//            {-0.19670719, 0.9851125, 0.04232101},
//            {-0.119770005, 0.77366745, 0.01729353},
//            {0.5822786, 0.48960224, 0.24117999},
//            {0.035166945, 0.85798544, 0.24074848},
//            {-0.09247135, 0.75247496, 0.23533802},
//            {-0.024287963, 0.9824017, 0.14671284},
//            {-0.1686875, 0.98302734, 0.20201232},
//            {-0.1646395, 0.7540736, 0.090583555},
//            {-0.33547747, 0.63215965, 0.19135737},
//            {-0.18665047, 0.9687785, 0.1996556},
//            {-0.11467207, 0.9839309, 0.18090159},
//            {-0.21024267, 0.9838614, 0.1267307},
//            {-0.024287963, 0.7552552, 0.14671284},
//            {0.43601954, 0.42613602, -0.14867124},
//            {-0.26836935, 0.6696235, 0.14820652},
//            {-0.2229559, 0.96863955, 0.11405098},
//            {0.037254192,  0.974617, 0.14438935},
//            {-0.08994135, 0.9165098, 0.11361948},
//            {0.5317419, 0.63028294, 0.3121797},
//            {0.038424317, 0.89148754, 0.1935813},
//            {0.07158624, 0.8688285, 0.19895856},
//            {0.39183316, 0.48960224, 0.110598914},
//            {0.07058689, 0.90059286, 0.19952284},
//            {0.07912563,  0.879602, 0.21625209},
//            {0.07912563,  0.858055, 0.21625209},
//            {0.28108257, 0.5321401, 0.40671325},
//            {-0.0154899005, 0.8688285, 0.15080555},
//            {-0.03261165, 0.9641911, 0.17648691},
//            {-0.21549243, 0.75984263, 0.14621496},
//            {-0.0065779896, 0.91074073, 0.12045723},
//            {-0.0347242, 0.91004574, 0.1752256},
//            {-0.090023585, 0.9643023, 0.19905813},
//            {-0.14606933, 0.9643023, 0.17565711},
//            {-0.1756513, 0.9640312, 0.12065638},
//            {-0.076152876, 0.8144746, 0.12799202},
//            {-0.2712156, 0.6375116, -0.18266082},
//            {0.3190325,  0.531167, 0.11063211},
//            {0.030676205, 0.75490767, 0.12570171},
//            {0.037254192, 0.76303995, 0.14438935},
//            {-0.29088628, 0.6309085, 0.14790781},
//            {-0.03261165, 0.7734659, 0.17648691},
//            {0.44711354, 0.5212971, 0.40233177},
//            {-0.0065779896, 0.8269162, 0.12045723},
//            {0.32289076, 0.63243765, 0.03661183},
//            {-0.090320855, 0.8273333, 0.19806236},
//            {-0.004427493, 0.72536755, 0.2025766},
//            {-0.14606933, 0.77335465, 0.17565711},
//            {-0.14477903, 0.8271248, 0.17409705},
//            {-0.170838, 0.7223788, 0.20051864},
//            {-0.1646395, 0.98358333, 0.090583555},
//            {0.01960747, 0.98226273, 0.045872655},
//            {0.01922797,  0.979691, 0.025392609},
//            {-0.07646914, 0.6310475, 0.019152336},
//            {0.0067677395, 0.8688285, 0.00853058},
//            {-0.23465714, 0.6187449, 0.08182062},
//            {-0.033022773, 0.8688285, 0.03477626},
//            {-0.23480892, 0.84790707, 0.06562248},
//            {-0.05836701, 0.9642467, 0.016895207},
//            {-0.045033928,  0.985043, -0.016696054},
//            {-0.16407025, 0.9096287, 0.059714116},
//            {-0.119770005, 0.9639895, 0.01729353},
//            {-0.11941582, 0.9099067, 0.018023778},
//            {-0.05836701, 0.77341026, 0.016895207},
//            {0.2934163, 0.37693256, 0.2171151},
//            {0.01960747, 0.7553942, 0.045872655},
//            {0.022833215, 0.72446394, 0.051150344},
//            {-0.21100168, 0.7218227, 0.123245426},
//            {-0.05818991, 0.82775027, 0.017127557},
//            {-0.11941582, 0.82775027, 0.018023778},
//            {-0.21973015, 0.5212971, 0.23606828},
//            {-0.19670719, 0.7525445, 0.04232101},
//            {0.34053746, 0.40772384, 0.31267762},
//            {-0.21232992, 0.75608927, 0.055166684},
//            {-0.21644118, 0.8688285, 0.08832644},
//            {-0.2660291, 0.9005233, 0.07521523},
//            {-0.2660291, 0.83713365, 0.07521523},
//            {-0.23851536, 0.87779474, 0.01842208},
//            {-0.2753901, 0.87876785, 0.09237599},
//            {-0.2753901, 0.8588891, 0.09237599},
//            {0.24466328, 0.5212971, 0.3649233},
//            {-0.39942312, 0.53026336, 0.34962136},
//            {-0.09651935, 0.8688285, 0.07010352},
//            {-0.2229559, 0.76901746, 0.11405098},
//            {-0.15382376, 0.7603987, -0.016297754},
//            {-0.18665047, 0.76887846, 0.1996556},
//            {-0.19252639, 0.8688285, 0.17695165},
//            {-0.24072912, 0.8688285, 0.20466773},
//            {-0.23965387, 0.9005233, 0.20353918},
//            {0.49910495, 0.5529919, 0.4022322},
//            {-0.24060263, 0.87869835, 0.13406634},
//            {-0.2661556, 0.8598622, 0.15411487},
//            {-0.23965387, 0.83713365, 0.20353918},
//            {-0.30075958, 0.5212971, -0.015159213},
//            {0.42396408, 0.56307036, -0.04338319},
//            {-0.107714586, 0.91060174, 0.12258158},
//            {-0.1007571, 0.63028294, -0.0197498},
//            {-0.13099055, 0.72223973, -0.017559057},
//            {-0.02536321, 0.7592171, -0.018056955},
//            {-0.054837663, 0.7670713, -0.032163966},
//            {-0.18443671, 0.8688285, 0.2360351},
//            {-0.19974321, 0.8777253, 0.24742028},
//            {-0.18285546, 0.9005233, 0.23626743},
//            {-0.19974321, 0.85993165, 0.24742028},
//            {0.40872085, 0.6376506, 0.40130278},
//            {-0.0916491, 0.87869835, 0.2656764},
//            {-0.13099055, 1.0154172, -0.017559057},
//            {0.27782518, 0.54395616, 0.35954607},
//            {0.108663335, 0.70055383, 0.23437543},
//            {0.096203096, 0.6988856,  0.250474},
//            {-0.30967152, 0.4793848, 0.28642198},
//            {0.08494462, 0.77749723, 0.23258303},
//            {0.09971979, 0.7462055, 0.2133311},
//            {-0.2580596, 0.56237525, -0.14883721},
//            {0.37461653, 0.61671525, -0.14906956},
//            {0.08469162, 0.69228256, 0.18349063},
//            {-0.10885309, 0.5312365, -0.13911167},
//            { 0.320677, 0.37783617, 0.03661183},
//            {0.10714534, 0.6732378, 0.22332218},
//            {0.39651364, 0.5122613, 0.46234465},
//            {0.32934222, 0.6166597, 0.2260108},
//            {0.43544397, 0.5212971, 0.35732877},
//            {-0.14345077, 0.5212971, 0.25562888},
//            {0.5510585, 0.50037575, 0.23158723},
//            {0.58379656, 0.5212971, 0.24098083},
//            {0.5510585, 0.54221845, 0.23158723},
//            {0.5159927, 0.53019387, 0.41338503},
//            {0.5822786, 0.5529919, 0.24117999},
//            {0.5916396, 0.5113577, 0.25834075},
//            {0.59543455, 0.53026336, 0.22571209},
//            {0.54160887, 0.5212971, 0.1982615},
//            {0.54160887, 0.5529919, 0.19985478},
//            {0.39562812, 0.5212971, 0.26013318},
//            {0.40619087, 0.4736158, 0.27958423},
//            {-0.17147048, 0.56300086, 0.3400618},
//            {0.39240238, 0.4669432, 0.2939568},
//            {-0.106006846, 0.63632995, 0.29269546},
//            {0.42396408, 0.47952378, 0.28854635},
//            { 0.472793, 0.5529919, 0.45003006},
//            {0.4574865, 0.5113577, 0.4619795},
//            {-0.032320693, 0.6239579, 0.4088708},
//            {0.48968074, 0.5124003, 0.46118286},
//            {0.33743823, 0.5212971, 0.4462129},
//            {0.3381972, 0.5529919, 0.44488513},
//            {0.34287772, 0.5122613, 0.46470135},
//            {0.3190325, 0.5114272, 0.44256166},
//            {0.48493698, 0.6354959, 0.36797708},
//            {0.4870875, 0.66774684, 0.3664834},
//            {0.4062731, 0.61677086, 0.3650229},
//            {0.40657037, 0.56279236, 0.3640271},
//            {0.42801207, 0.6372336, 0.40734392},
//            {0.5272512, 0.66830283, 0.2892102},
//            {0.40170014, 0.66809434, 0.40209946},
//            {0.34886116, 0.6166597, 0.3424517},
//            { 0.320677,  0.664758, 0.36854136},
//            {0.5481869, 0.5113577, 0.3875277},
//            {0.5559034, 0.48960224, 0.36950395},
//            {0.4709714,  0.498082, 0.41760054},
//            {0.50353247, 0.5312365, 0.38307986},
//            {0.5380479, 0.50197434, 0.3413563},
//            {0.46102855, 0.4795933, 0.3400618},
//            {0.4078986,  0.531167, 0.43164113},
//            {0.39784187, 0.42384928, 0.41375014},
//            { 0.328773, 0.42002645, 0.3871626},
//            {0.30524403, 0.4241273, 0.36382794},
//            {0.27782518, 0.49863806, 0.35954607},
//            {0.2912658, 0.5212971, 0.39091343},
//            {0.31030402, 0.5133734, 0.3929382},
//            {0.3509737, 0.48007986, 0.34119037},
//            {0.32289076, 0.41015655, 0.36854136},
//            {0.39651364, 0.5122613, 0.13041511},
//            { 0.320677, 0.37783617, 0.36854136},
//            {-0.22453715, 0.5943482, -0.06505818},
//            {0.37968916, 0.4061252, 0.35008606},
//            {0.24566263, 0.5530615, 0.3654876},
//            {0.45445055, 0.4073763, -0.12268115},
//            {0.5824051, 0.5123308, 0.32007965},
//            {0.23307589, 0.5123308, 0.34962136},
//            {-0.33547747, 0.63215965, -0.14057215},
//            {0.4870875, 0.37484738, 0.03455386},
//            {0.23712389, 0.5105236, 0.38221687},
//            {0.23307589, 0.53026336, 0.34962136},
//            {0.47940895, 0.5212971, 0.3163023},
//            {0.43092155, 0.4061947, 0.34686637},
//            {0.3507207, 0.5212971, 0.39791712},
//            {0.49190077, 0.42609435, 0.28662115},
//            {0.56361985, 0.5212971, 0.31788892},
//            {0.48493698, 0.40709826, 0.36797708},
//            {0.47741026, 0.41453546, 0.38643235},
//            {0.23712389, 0.5105236, 0.050287317},
//            {0.4870875, 0.37484738, 0.3664834},
//            {0.32131582, 0.61703503, 0.28745097},
//            {-0.15932651, 0.40668124, -0.020745596},
//            {0.49190077, 0.6164998, 0.28662115},
//            {0.4731725, 0.63591295, 0.31118393},
//            {-0.15114197, 0.42601788, 0.22514781},
//            {0.32934222, 0.4259345, 0.2260108},
//            {0.32282752, 0.4793848, 0.28642198},
//            {0.47007325, 0.62972695, 0.14966701},
//            {0.45445055, 0.63521785, 0.20924836},
//            {-0.1964795, 0.42613602, -0.14867124},
//            {0.36412966, 0.6696235, 0.14820652},
//            {-0.2580596, 0.56237525, 0.18309233},
//            {0.37461653, 0.61671525, 0.18285997},
//            {0.40574813, 0.6355654, 0.18767294},
//            {-0.2578825, 0.61671525, 0.18285997},
//            {0.37443942, 0.56237525, 0.18309233},
//            {0.43566534, 0.56237525, 0.18398854},
//            {0.45027605, 0.63591295, 0.1501649},
//            {0.2934163, 0.66566163, -0.11481442},
//            {0.44724005, 0.6678859, 0.1484057},
//            {0.4412313, 0.42211163, 0.13459744},
//            {0.37108716, 0.4195399, 0.13380079},
//            {0.37461653, 0.42587888, 0.18285997},
//            {-0.26836935, 0.6696235, -0.183723},
//            {0.43601954, 0.42613602, 0.18325828},
//            {0.37443942, 0.4802189, 0.18309233},
//            {0.35495847, 0.4061252, 0.20689169},
//            {0.36412966, 0.3729707, 0.14820652},
//            {0.36128342, 0.40508258, 0.1492687},
//            {0.3416127, 0.41168568, 0.14790781},
//            {0.45027605, 0.40668124, 0.1501649},
//            {0.48135704, 0.6165763, -0.10678172},
//            {0.44724005, 0.37470832, 0.1484057},
//            {0.5824051, 0.5123308, -0.011849892},
//            {0.45445055, 0.4073763, 0.20924836},
//            {0.34927228, 0.5212971, 0.20074102},
//            {0.30948177, 0.5212971, 0.17449534},
//            {0.30568677, 0.49641386, 0.16480301},
//            {0.30568677, 0.5461803, 0.16480301},
//            {0.36128342, 0.6375116, 0.1492687},
//            {-0.18538547, 0.5212971, 0.40233177},
//            {0.3416127, 0.6309085, 0.14790781},
//            {-0.0021504953, 0.40334493, 0.067813195},
//            {-0.1318128, 0.5212971, 0.40199986},
//            {0.28190482, 0.42649055, 0.21937221},
//            {0.29664204, 0.4078628, 0.21183743},
//            {0.2855733, 0.4073763, 0.29166648},
//            {0.2789953, 0.41550857, 0.3103541},
//            {0.2845613, 0.37512535, 0.28940937},
//            {-0.08144058, 0.50037575, 0.23158723},
//            {0.48904824, 0.5212971, 0.25562888},
//            {-0.30315676, 0.4259345, 0.2260108},
//            {0.48135704, 0.42601788, 0.22514781},
//            {0.48031977, 0.48049688, 0.22567888},
//            {0.5164987, 0.42176414, 0.19105865},
//            {0.5129567, 0.40501308, 0.20828578},
//            {0.51067966, 0.37359628, 0.20370515},
//            {-0.23079889, 0.37449983, 0.40209946},
//            { 0.480889, 0.40654224, 0.25654835},
//            {0.5264922, 0.40626422, 0.29269546},
//            {0.5317419, 0.41231126, 0.3121797},
//            {0.24699087, 0.6520384, 0.011385182},
//            {-0.093293615, 0.6211081, -0.051913787},
//            {0.55685216, 0.5114272, 0.30003113},
//            {0.55685216,  0.531167, 0.30003113},
//            {-0.025869206, 0.36720166,  0.404091},
//            {0.5824051, 0.53026336, 0.32007965},
//            {0.3327577, 0.63591295, 0.25780967},
//            {0.32282752, 0.56320935, 0.28642198},
//            {0.2855733, 0.63521785, 0.29166648},
//            {0.2845613, 0.66746885, 0.28940937},
//            {0.27406183, 0.62089956, 0.2774267},
//            {0.44724005, 0.6678859, -0.18352382},
//            {0.2934163, 0.66566163, 0.2171151},
//            {-0.06687412, 0.64807653, 0.36558717},
//            {0.44724005, 0.37470832, -0.18352382},
//            {0.48135704, 0.6165763, 0.22514781},
//            {-0.20853493, 0.47952378, 0.28854635},
//            {0.5264922, 0.63632995, 0.29269546},
//            {0.48031977, 0.56209725, 0.22567888},
//            {-0.075646885, 0.5114272, -0.031898417},
//            {0.53920543, 0.6211081, 0.28001574},
//            {0.5129567, 0.6375811, 0.20828578},
//            {-0.22453715, 0.5943482, 0.26687133},
//            {-0.050220422, 0.48960224, 0.24117999},
//            {-0.037064448, 0.5123308, 0.22571209},
//            {-0.04085943, 0.5312365, 0.25834075},
//            {-0.35141644, 0.5104541, 0.40671325},
//            {-0.077734135, 0.5123308, 0.18438683},
//            {-0.09089011, 0.5212971, 0.1982615},
//            {-0.10885309, 0.5113577, 0.19281788},
//            {-0.09089011, 0.48960224, 0.19985478},
//            {-0.10885309, 0.5312365, 0.19281788},
//            {-0.20448694, 0.4053606, 0.40734392},
//            {-0.093293615, 0.42148608, 0.28001574},
//            {-0.23687088, 0.5212971, 0.26013318},
//            {-0.22630815, 0.4736158, 0.27958423},
//            {-0.3584372, 0.4216946, 0.2774267},
//            {-0.159706, 0.48960224, 0.45003006},
//            {-0.23940088, 0.5212971, 0.44335827},
//            {-0.31182203, 0.37783617, 0.03661183},
//            {0.52364594, 0.5312365, -0.13911167},
//            {-0.23598538, 0.5122613, 0.46234465},
//            {-0.23598538, 0.5303329, 0.46234465},
//            {-0.29430178, 0.5529919, 0.44488513},
//            {-0.29430178, 0.48960224, 0.44488513},
//            {-0.31346652, 0.5114272, 0.44256166},
//            {-0.31346652,  0.531167, 0.44256166},
//            {-0.17018019, 0.61677086, 0.34162188},
//            {-0.14756203, 0.6354959, 0.36797708},
//            {-0.15508877, 0.6280587, 0.38643235},
//            {-0.12959905, 0.6212471, 0.36562037},
//            {-0.22622591, 0.61677086, 0.3650229},
//            {-0.335857, 0.63473135, 0.21183743},
//            {-0.25280985,  0.636469, 0.35008606},
//            {-0.20157744, 0.6363995, 0.34686637},
//            {-0.22592865, 0.56279236, 0.3640271},
//            {-0.22377814, 0.6376506, 0.40130278},
//            {-0.2815253, 0.5625143, 0.34119037},
//            {-0.30960828, 0.63243765, 0.36854136},
//            {-0.327255, 0.61846685, 0.36382794},
//            {-0.30372605, 0.6225678, 0.3871626},
//            {-0.31182203,  0.664758, 0.36854136},
//            {-0.07552038, 0.5212971, 0.3706325},
//            {-0.027197465, 0.45109576, 0.3564259},
//            {-0.084312126, 0.5113577, 0.3875277},
//            {-0.056608666, 0.5123308, 0.3706989},
//            {-0.076595634, 0.48960224, 0.36950395},
//            {-0.084312126, 0.5312365, 0.3875277},
//            {-0.1615276, 0.54451215, 0.41760054},
//            {-0.13339405, 0.48960224, 0.4022322},
//            {-0.13339405, 0.5529919, 0.4022322},
//            {-0.16242574, 0.41286728, 0.14966701},
//            {-0.12896655, 0.5113577, 0.38307986},
//            {-0.116506316, 0.5124003, 0.41338503},
//            {-0.12896655, 0.5312365, 0.38307986},
//            {-0.12959905, 0.42134705, 0.36562037},
//            {-0.09445107, 0.54061985, 0.3413563},
//            {0.49190077, 0.6164998, -0.04530838},
//            {-0.15932651, 0.40668124, 0.31118393},
//            {-0.26822385, 0.50176585, 0.4243387},
//            {-0.24066588, 0.48960224, 0.44252843},
//            {-0.2246004,  0.531167, 0.43164113},
//            {-0.19705509, 0.5212971, 0.35732877},
//            {-0.22592865, 0.47980183, 0.3640271},
//            {-0.30372605, 0.42002645, 0.3871626},
//            {-0.327255, 0.4241273, 0.36382794},
//            {-0.35467383, 0.49863806, 0.35954607},
//            {-0.33908272, 0.5212971,  0.335448},
//            {-0.3412332, 0.5212971, 0.39091343},
//            {-0.34287772, 0.48960224,  0.389984},
//            {-0.322195, 0.5133734, 0.3929382},
//            {-0.2815253, 0.48007986, 0.34119037},
//            {-0.29196155, 0.40772384, 0.31267762},
//            {0.5087759, 0.5212971, 0.01098687},
//            {-0.31118318, 0.42555913, 0.28745097},
//            {-0.38783577, 0.5212971, 0.3649233},
//            {-0.38683638, 0.5530615, 0.3654876},
//            {-0.17804848, 0.4073763, -0.12268115},
//            {-0.050093923, 0.5123308, 0.32007965},
//            {-0.39942312, 0.5123308, 0.34962136},
//            {-0.38683638, 0.48953274, 0.3654876},
//            {-0.20157744, 0.4061947, 0.34686637},
//            {-0.22377814, 0.40494356, 0.40130278},
//            {-0.06887915, 0.5212971, 0.31788892},
//            {-0.07008089, 0.5529919, 0.31888467},
//            {-0.15508877, 0.41453546, 0.38643235},
//            {-0.31118318, 0.61703503, 0.28745097},
//            {-0.29196155, 0.63487035, 0.31267762},
//            {0.4731725, 0.40668124, -0.020745596},
//            {-0.14059821, 0.6164998, 0.28662115},
//            {-0.30081654, 0.48007986, 0.22630955},
//            {-0.16242574, 0.62972695, 0.14966701},
//            {-0.17804848, 0.63521785, 0.20924836},
//            {-0.2267509, 0.6355654, 0.18767294},
//            {-0.19683369, 0.56237525, 0.18398854},
//            {-0.18222298, 0.63591295, 0.1501649},
//            {-0.33908272, 0.66566163, -0.11481442},
//            {-0.18525895, 0.6678859, 0.1484057},
//            {-0.1912677, 0.62048256, 0.13459744},
//            {-0.2578825, 0.42587888, 0.18285997},
//            {-0.26836935, 0.3729707, 0.14820652},
//            {-0.2712156, 0.40508258, 0.1492687},
//            {-0.29088628, 0.41168568, 0.14790781},
//            {-0.2267509, 0.40702876, 0.18767294},
//            {-0.18222298, 0.40668124, 0.1501649},
//            {-0.15114197, 0.6165763, -0.10678172},
//            {-0.18525895, 0.37470832, 0.1484057},
//            {-0.32301724, 0.5212971, 0.17449534},
//            {-0.32681224, 0.49641386, 0.16480301},
//            {-0.27754056,  0.636469, 0.20689169},
//            {-0.335857, 0.4078628, 0.21183743},
//            {-0.33908272, 0.37693256, 0.2171151},
//            {-0.33547747, 0.41043457, 0.19135737},
//            {-0.3479377, 0.37512535, 0.28940937},
//            {-0.14307128, 0.47910678, 0.28635558},
//            {-0.09980834, 0.5212971, 0.2542912},
//            {-0.15217926, 0.48049688, 0.22567888},
//            {-0.11954232, 0.40501308, 0.20828578},
//            {-0.10391958, 0.40855792, 0.22113144},
//            {-0.106006846, 0.40626422, 0.29269546},
//            {-0.1007571, 0.41231126, 0.3121797},
//            {0.53920543, 0.6211081, -0.051913787},
//            {-0.075646885, 0.5114272, 0.30003113},
//            {-0.075646885,  0.531167, 0.30003113},
//            {-0.050093923, 0.53026336, 0.32007965},
//            {-0.30081654, 0.5625143, 0.22630955},
//            {-0.3469257, 0.63521785, 0.29166648},
//            {0.4870875, 0.66774684, 0.03455386},
//            {-0.35350367, 0.6270856, 0.3103541},
//            {-0.18525895, 0.6678859, -0.18352382},
//            {-0.33908272, 0.66566163, 0.2171151},
//            {-0.35059422, 0.61610365, 0.21937221},
//            {-0.18525895, 0.37470832, -0.18352382},
//            {-0.15114197, 0.6165763, 0.22514781},
//            { -0.15161, 0.63605195, 0.25654835},
//            {-0.14307128, 0.5634874, 0.28635558},
//            {0.55685216, 0.5114272, -0.031898417},
//            {-0.093293615, 0.6211081, 0.28001574},
//            {-0.11954232, 0.6375811, 0.20828578},
//            {-0.121819325, 0.6689979, 0.20370515},
//            {-0.10391958, 0.63403624, 0.22113144},
//            {0.28108257, 0.5321401, 0.07478373},
//            {-0.116000324, 0.62083006, 0.19105865},
//            {0.39240238,  0.575651, -0.03797275},
//            {0.32934222, 0.6166597, -0.10591872},
//            {-0.14345077, 0.5212971, -0.07630064},
//            {0.5510585, 0.50037575, -0.10034229},
//            {0.5510585, 0.54221845, -0.10034229},
//            {0.5822786, 0.48960224, -0.09074953},
//            {0.5822786, 0.5529919, -0.09074953},
//            {0.23307589, 0.5123308, 0.017691841},
//            {0.5916396, 0.5113577, -0.07358878},
//            {0.59543455, 0.5123308, -0.106217444},
//            {0.59543455, 0.53026336, -0.106217444},
//            {0.55476487, 0.53026336, -0.14754269},
//            {0.54160887, 0.5212971, -0.13366802},
//            {0.54160887, 0.5529919, -0.13207476},
//            {0.52364594, 0.5113577, -0.13911167},
//            {0.54160887, 0.48960224, -0.13207476},
//            {0.53920543, 0.42148608, -0.051913787},
//            {0.39562812, 0.5212971, -0.07179636},
//            {0.41276884, 0.5212971, -0.09586125},
//            {0.32671103, 0.5212971, -0.07540443},
//            {0.40796188, 0.44824603, -0.06505818},
//            {0.39240238, 0.4669432, -0.03797275},
//            {-0.106006846, 0.63632995, -0.039234072},
//            {0.42396408, 0.47952378, -0.04338319},
//            {0.27406183, 0.4216946, -0.054502834},
//            {0.5624181, 0.5529919, -0.013044839},
//            { 0.472793, 0.5529919, 0.118100524},
//            {0.4574865, 0.5113577, 0.13004999},
//            {-0.032320693, 0.6239579, 0.07694126},
//            {0.48968074, 0.5124003, 0.12925336},
//            { 0.472793, 0.48960224, 0.118100524},
//            {0.39309815, 0.5212971, 0.111428745},
//            {-0.108663335, 0.6895718, 0.097554095},
//            {0.33743823, 0.5212971, 0.11428334},
//            {0.3381972, 0.5529919, 0.112955615},
//            {0.3381972, 0.48960224, 0.112955615},
//            {0.3190325, 0.5114272, 0.11063211},
//            {0.46102855, 0.56300086, 0.008132277},
//            {0.48493698, 0.6354959, 0.036047544},
//            {0.47741026, 0.6280587, 0.054502834},
//            {   0.5029, 0.6212471, 0.033690847},
//            {0.4062731, 0.61677086, 0.033093374},
//            {-0.09980834, 0.5212971, -0.07763832},
//            {0.43092155, 0.6363995, 0.014936824},
//            {0.40657037, 0.56279236, 0.03209759},
//            {0.40872085, 0.6376506, 0.069373265},
//            {0.39784187, 0.6187449, 0.08182062},
//            {0.42801207, 0.6372336, 0.07541439},
//            {0.40170014, 0.66809434, 0.0701699},
//            {0.34886116, 0.6166597, 0.010522162},
//            {0.30524403, 0.61846685, 0.031898428},
//            { 0.328773, 0.6225678, 0.05523307},
//            {0.55697864, 0.5212971, 0.038702987},
//            {0.5559034, 0.5529919, 0.037574425},
//            {0.5481869, 0.5113577, 0.055598192},
//            {0.57589036, 0.5123308, 0.03876937},
//            {0.5559034, 0.48960224, 0.037574425},
//            {0.5481869, 0.5312365, 0.055598192},
//            {0.57589036, 0.53026336, 0.03876937},
//            {0.4709714,  0.498082, 0.085671015},
//            {0.5006862, 0.5212971, 0.07007033},
//            {0.49910495, 0.48960224, 0.07030268},
//            {0.50353247, 0.5113577, 0.051150344},
//            {0.5159927, 0.5124003, 0.08145551},
//            {0.50353247, 0.5312365, 0.051150344},
//            {0.5159927, 0.53019387, 0.08145551},
//            {0.5380479, 0.50197434, 0.0094268},
//            {0.5380479, 0.54061985, 0.0094268},
//            {0.46102855, 0.4795933, 0.008132277},
//            {0.4078986, 0.5114272, 0.09971163},
//            {0.39183316, 0.5529919, 0.110598914},
//            {0.3642752, 0.50176585, 0.09240918},
//            {0.4078986,  0.531167, 0.09971163},
//            {0.39784187, 0.42384928, 0.08182062},
//            {0.37627366, 0.5212971, 0.025399245},
//            {0.43544397, 0.5212971, 0.025399245},
//            {-0.31118318, 0.42555913, -0.044478558},
//            {0.4062731, 0.42582327, 0.033093374},
//            {0.18861121, 0.6860965, 0.075049266},
//            {0.40657037, 0.47980183, 0.03209759},
//            { 0.328773, 0.42002645, 0.05523307},
//            {0.30524403, 0.4241273, 0.031898428},
//            {0.2912658, 0.5212971, 0.05898388},
//            {0.27782518, 0.54395616, 0.027616538},
//            {0.2896213, 0.5529919, 0.058054473},
//            {0.31030402, 0.5133734, 0.06100865},
//            {0.28108257, 0.5104541, 0.07478373},
//            {0.34886116, 0.4259345, 0.010522162},
//            {0.3509737, 0.48007986, 0.009260837},
//            {0.34053746, 0.40772384, -0.019251922},
//            {-0.20853493, 0.56307036, -0.04338319},
//            {0.3317394, 0.5212971, -0.015159213},
//            {-0.22622591, 0.42582327, 0.033093374},
//            {0.32131582, 0.42555913, -0.044478558},
//            {0.37968916, 0.4061252, 0.01815654},
//            {0.24466328, 0.5212971, 0.032993797},
//            {0.24566263, 0.5530615, 0.033558074},
//            {0.24566263, 0.48953274, 0.033558074},
//            {-0.04870242, 0.5212971, -0.090948686},
//            {0.23307589, 0.53026336, 0.017691841},
//            {0.23712389, 0.5320706, 0.050287317},
//            {0.47940895, 0.5212971, -0.015627235},
//            {0.43092155, 0.4061947, 0.014936824},
//            {0.40170014, 0.37449983, 0.0701699},
//            {0.3507207, 0.5212971, 0.06598759},
//            {0.49190077, 0.42609435, -0.04530838},
//            {0.56361985, 0.5212971, -0.014040614},
//            {0.5624181, 0.48960224, -0.013044839},
//            {0.48493698, 0.40709826, 0.036047544},
//            {0.47741026, 0.41453546, 0.054502834},
//            {0.34053746, 0.63487035, -0.019251922},
//            {0.4731725, 0.63591295, -0.020745596},
//            {-0.15114197, 0.42601788, -0.10678172},
//            {0.32934222, 0.4259345, -0.10591872},
//            {0.32282752, 0.4793848, -0.045507535},
//            {0.47007325, 0.62972695, -0.18226251},
//            {0.45445055, 0.63521785, -0.12268115},
//            {0.37108716, 0.62305427, -0.19812873},
//            {0.43601954, 0.6164582, -0.14867124},
//            {0.40574813, 0.6355654, -0.14425656},
//            {-0.2578825, 0.61671525, -0.14906956},
//            {0.37443942, 0.56237525, -0.14883721},
//            {0.43566534, 0.56237525, -0.14794098},
//            {0.45027605, 0.63591295, -0.18176462},
//            {0.45842263, 0.5212971, -0.12739457},
//            {0.4412313, 0.42211163, -0.1973321},
//            {0.37108716, 0.4195399, -0.19812873},
//            {0.4050587, 0.5212971, -0.14763893},
//            {0.37461653, 0.42587888, -0.14906956},
//            {0.37443942, 0.4802189, -0.14883721},
//            {0.43566534, 0.4802189, -0.14794098},
//            {0.36128342, 0.40508258, -0.18266082},
//            {0.45027605, 0.40668124, -0.18176462},
//            {0.34927228, 0.5212971, -0.13118851},
//            {0.30568677, 0.49641386, -0.1671265},
//            {-0.31346652,  0.531167, 0.11063211},
//            {0.36128342, 0.6375116, -0.18266082},
//            {0.3416127, 0.6309085, -0.18402171},
//            {0.28190482, 0.42649055, -0.1125573},
//            {0.2934163, 0.37693256, -0.11481442},
//            {0.29702154, 0.41043457, -0.14057215},
//            {0.2855733, 0.4073763, -0.040263046},
//            {0.2789953, 0.41550857, -0.021575417},
//            {0.2845613, 0.37512535, -0.042520165},
//            {-0.08144058, 0.50037575, -0.10034229},
//            {0.48904824, 0.5212971, -0.07630064},
//            {-0.30315676, 0.4259345, -0.10591872},
//            {0.48135704, 0.42601788, -0.10678172},
//            {-0.20157744, 0.6363995, 0.014936824},
//            {0.53269064, 0.5212971, -0.07763832},
//            {0.5164987, 0.42176414, -0.14087088},
//            {0.5129567, 0.40501308, -0.12364375},
//            {0.51067966, 0.37359628, -0.12822439},
//            {0.5285794, 0.40855792, -0.11079808},
//            {0.5317419, 0.41231126, -0.0197498},
//            {0.55685216,  0.531167, -0.031898417},
//            {0.5824051, 0.53026336, -0.011849892},
//            {0.3327577, 0.63591295, -0.07411986},
//            {0.32282752, 0.56320935, -0.045507535},
//            {0.2855733, 0.63521785, -0.040263046},
//            {0.2845613, 0.66746885, -0.042520165},
//            {0.28190482, 0.61610365, -0.1125573},
//            { 0.480889, 0.63605195, -0.0753812},
//            {0.48942775, 0.5634874, -0.045573927},
//            {-0.20853493, 0.47952378, -0.04338319},
//            {0.5264922, 0.63632995, -0.039234072},
//            {0.5317419, 0.63028294, -0.0197498},
//            {0.51067966, 0.6689979, -0.12822439},
//            {0.5285794, 0.63403624, -0.11079808},
//            {-0.22630815, 0.5689784, -0.052345287},
//            {-0.050220422, 0.48960224, -0.09074953},
//            {-0.39942312, 0.5123308, 0.017691841},
//            {-0.04085943, 0.5113577, -0.07358878},
//            {-0.077734135, 0.5123308, -0.14754269},
//            {-0.077734135, 0.53026336, -0.14754269},
//            {-0.09089011, 0.5212971, -0.13366802},
//            {-0.09089011, 0.5529919, -0.13207476},
//            {-0.10885309, 0.5113577, -0.13911167},
//            {-0.09089011, 0.48960224, -0.13207476},
//            {-0.093293615, 0.42148608, -0.051913787},
//            {-0.22630815, 0.4736158, -0.052345287},
//            {-0.30578798, 0.5212971, -0.07540443},
//            {-0.22453715, 0.44824603, -0.06505818},
//            {-0.159706, 0.48960224, 0.118100524},
//            {-0.17501247, 0.5312365, 0.13004999},
//            {-0.23940088, 0.5212971, 0.111428745},
//            {-0.23598538, 0.5303329, 0.13041511},
//            {-0.29506078, 0.5212971, 0.11428334},
//            {-0.29430178, 0.5529919, 0.112955615},
//            {-0.2896213, 0.5122613, 0.13277182},
//            {-0.29430178, 0.48960224, 0.112955615},
//            {-0.2896213, 0.5303329, 0.13277182},
//            {-0.31346652, 0.5114272, 0.11063211},
//            {-0.17018019, 0.61677086, 0.009692338},
//            {-0.14756203, 0.6354959, 0.036047544},
//            {-0.15508877, 0.6280587, 0.054502834},
//            {-0.12959905, 0.6212471, 0.033690847},
//            {-0.22622591, 0.61677086, 0.033093374},
//            {-0.25280985,  0.636469, 0.01815654},
//            {-0.22592865, 0.56279236, 0.03209759},
//            {-0.22377814, 0.6376506, 0.069373265},
//            {-0.2815253, 0.5625143, 0.009260837},
//            {-0.30960828, 0.63243765, 0.03661183},
//            {-0.327255, 0.61846685, 0.031898428},
//            {-0.30372605, 0.6225678, 0.05523307},
//            {-0.31182203,  0.664758, 0.03661183},
//            {-0.076595634, 0.5529919, 0.037574425},
//            {-0.027197465, 0.45109576, 0.024496399},
//            {-0.084312126, 0.5113577, 0.055598192},
//            {-0.056608666, 0.5123308, 0.03876937},
//            {-0.076595634, 0.48960224, 0.037574425},
//            {-0.084312126, 0.5312365, 0.055598192},
//            {-0.0021504953, 0.40334493, 0.3997427},
//            {-0.1318128, 0.5212971, 0.07007033},
//            {-0.1615276, 0.54451215, 0.085671015},
//            {-0.13339405, 0.48960224, 0.07030268},
//            {-0.13339405, 0.5529919, 0.07030268},
//            {-0.12896655, 0.5113577, 0.051150344},
//            {-0.116506316, 0.5124003, 0.08145551},
//            {-0.12896655, 0.5312365, 0.051150344},
//            {-0.12959905, 0.42134705, 0.033690847},
//            {-0.09445107, 0.54061985, 0.0094268},
//            {-0.17147048, 0.4795933, 0.008132277},
//            {-0.24066588, 0.5529919, 0.110598914},
//            {-0.26822385, 0.5408284, 0.09240918},
//            {-0.24066588, 0.48960224, 0.110598914},
//            {-0.2246004,  0.531167, 0.09971163},
//            {-0.23465714, 0.42384928, 0.08182062},
//            {-0.25622535, 0.5212971, 0.025399245},
//            {-0.19705509, 0.5212971, 0.025399245},
//            {-0.22592865, 0.47980183, 0.03209759},
//            {-0.18538547, 0.5212971, 0.07040225},
//            {-0.327255, 0.4241273, 0.031898428},
//            {-0.33908272, 0.5212971, 0.0035184533},
//            {0.28508627, 0.40258732, 0.13572598},
//            {-0.3412332, 0.5212971, 0.05898388},
//            {-0.34287772, 0.48960224, 0.058054473},
//            {-0.322195, 0.5133734, 0.06100865},
//            {-0.2815253, 0.48007986, 0.009260837},
//            {-0.29196155, 0.40772384, -0.019251922},
//            {-0.30960828, 0.41015655, 0.03661183},
//            {-0.38783577, 0.5212971, 0.032993797},
//            {-0.38683638, 0.5530615, 0.033558074},
//            {-0.38683638, 0.48953274, 0.033558074},
//            {-0.15309007, 0.5212971, -0.015627235},
//            {-0.106006846, 0.40626422, -0.039234072},
//            {-0.22377814, 0.40494356, 0.069373265},
//            {-0.23079889, 0.37449983, 0.0701699},
//            {-0.2817783, 0.5212971, 0.06598759},
//            {-0.14756203, 0.40709826, 0.036047544},
//            {-0.15508877, 0.41453546, 0.054502834},
//            {-0.29196155, 0.63487035, -0.019251922},
//            {-0.30081654, 0.48007986, -0.10561997},
//            {-0.2267509, 0.6355654, -0.14425656},
//            {-0.19683369, 0.56237525, -0.14794098},
//            {-0.18222298, 0.63591295, -0.18176462},
//            {-0.1912677, 0.62048256, -0.1973321},
//            {-0.1740764, 0.5212971, -0.12739457},
//            {-0.22744031, 0.5212971, -0.14763893},
//            {-0.19683369, 0.4802189, -0.14794098},
//            {-0.26836935, 0.3729707, -0.183723},
//            {-0.29088628, 0.41168568, -0.18402171},
//            {-0.2267509, 0.40702876, -0.14425656},
//            {-0.16242574, 0.41286728, -0.18226251},
//            {-0.18222298, 0.40668124, -0.18176462},
//            {-0.32301724, 0.5212971, -0.15743418},
//            {-0.32681224, 0.5461803, -0.1671265},
//            {-0.27754056,  0.636469, -0.12503785},
//            {-0.29088628, 0.6309085, -0.18402171},
//            {-0.35059422, 0.42649055, -0.1125573},
//            {-0.335857, 0.4078628, -0.1200921},
//            {-0.33908272, 0.37693256, -0.11481442},
//            {-0.33547747, 0.41043457, -0.14057215},
//            {-0.29974127, 0.40668124, -0.07411986},
//            {-0.3479377, 0.37512535, -0.042520165},
//            {-0.14307128, 0.47910678, -0.045573927},
//            {-0.15217926, 0.48049688, -0.106250644},
//            {-0.10391958, 0.40855792, -0.11079808},
//            { -0.15161, 0.40654224, -0.0753812},
//            {-0.075646885,  0.531167, -0.031898417},
//            {-0.050093923, 0.53026336, -0.011849892},
//            {-0.29974127, 0.63591295, -0.07411986},
//            {-0.30967152, 0.56320935, -0.045507535},
//            {-0.3469257, 0.63521785, -0.040263046},
//            {-0.35350367, 0.6270856, -0.021575417},
//            {-0.3479377, 0.66746885, -0.042520165},
//            {-0.335857, 0.63473135, -0.1200921},
//            {-0.35059422, 0.61610365, -0.1125573},
//            { -0.15161, 0.63605195, -0.0753812},
//            {-0.14307128, 0.5634874, -0.045573927},
//            {-0.11954232, 0.6375811, -0.12364375},
//            {-0.121819325, 0.6689979, -0.12822439},
//            {-0.10391958, 0.63403624, -0.11079808},
//            {-0.116000324, 0.62083006, -0.14087088},
//            {0.21652973, 0.64392006, 0.04736634},
//            {0.2313049, 0.61262834, 0.06661826},
//            {0.20910417, 0.7168878, 0.057357423},
//            {0.20758617, 0.6895718, 0.06841067},
//            {0.23155789, 0.6978431, 0.01752588},
//            {0.24041288, 0.6719867, 0.026089657},
//            {0.25597236, 0.6818566, 0.037574425},
//            {-0.03116323, 0.40258732, 0.030238783},
//            {-0.021821218, 0.3973674, 0.06572205},
//            {-0.025869206, 0.36720166, 0.07216148},
//            {-0.071156144, 0.3629618, 0.031466916},
//            {-0.07646914, 0.4115467, 0.019152336},
//            {-0.03116323, 0.40258732, 0.3621683},
//            {-0.021821218, 0.3973674, 0.39765155},
//            {-0.06687412, 0.39451763, 0.36558717},
//            {-0.07646914, 0.4115467, 0.35108188},
//            {-0.03116323, 0.64000684, 0.030238783},
//            {-0.027197465, 0.59149843, 0.024496399},
//            {-0.06687412, 0.64807653, 0.033657648},
//            {-0.0021504953, 0.63924927, 0.067813195},
//            {-0.025869206, 0.67539257, 0.07216148},
//            {-0.071156144, 0.6796324, 0.031466916},
//            {-0.07286389, 0.63799816, 0.052212514},
//            {-0.03116323, 0.64000684, 0.3621683},
//            {-0.027197465, 0.59149843, 0.3564259},
//            {-0.0021504953, 0.63924927, 0.3997427},
//            {-0.071156144, 0.6796324, 0.36339647},
//            {-0.07286389, 0.63799816, 0.38414204},
//            {-0.08469162, 0.6978431, 0.14843889},
//            {-0.075836636, 0.6719867, 0.13987511},
//            {-0.10714534, 0.7168878, 0.10860734},
//            {-0.060277157, 0.6818566, 0.12839034},
//            {-0.1276383, 0.6860965, 0.09091549},
//            {-0.096203096, 0.6912399, 0.0814555},
//            {-0.08494462, 0.61262834, 0.0993465},
//            {0.2944283, 0.3973674, 0.10024271},
//            {0.2903803, 0.36720166, 0.09380328},
//            {0.24937539, 0.39451763, 0.13230711}
//    }};
//    voronoi::VoronoiDiagram::BoolVector bools{false, true, false, true,
//                                                    false, true, false, true,
//                                                    false, false, false, true,
//                                                    false, false, false, true,
//                                                    false, false, true, true,
//                                                    false, false, false, false,
//                                                    false, false, true, false,
//                                                    true, true, true, true,
//                                                    false, false, false, true,
//                                                    false, false, false, false,
//                                                    false, true, true, false,
//                                                    false, true, false, false,
//                                                    true, false, true, false,
//                                                    false, false, false, true,
//                                                    false, false, false, false,
//                                                    true, false, false, true,
//                                                    false, false, false, false,
//                                                    false, false, false, true,
//                                                    false, false, false, true,
//                                                    true, false, false, true,
//                                                    false, false, false, true,
//                                                    false, true, true, false,
//                                                    true, true, false, false,
//                                                    true, false, false, false,
//                                                    true, true, true, false,
//                                                    false, false, true, false,
//                                                    false, true, false, true,
//                                                    false, true, true, false,
//                                                    false, false, false, false,
//                                                    false, true, false, false,
//                                                    false, true, true, false,
//                                                    false, true, false, false,
//                                                    false, true, false, true,
//                                                    false, false, true, false,
//                                                    true, false, false, true,
//                                                    false, false, false, true,
//                                                    false, false, false, true,
//                                                    true, false, false, false,
//                                                    false, false, false, false,
//                                                    true, false, true, false,
//                                                    false, false, true, true,
//                                                    true, false, false, false,
//                                                    false, true, true, false,
//                                                    false, false, false, true,
//                                                    false, true, false, true,
//                                                    false, false, false, true,
//                                                    false, true, false, false,
//                                                    true, false, false, false,
//                                                    true, true, false, false,
//                                                    true, false, true, true,
//                                                    false, true, false, true,
//                                                    false, false, true, false,
//                                                    false, false, false, true,
//                                                    false, true, true, false,
//                                                    false, false, true, false,
//                                                    false, true, true, false,
//                                                    false, false, false, false,
//                                                    true, true, false, false,
//                                                    false, true, false, false,
//                                                    true, false, false, false,
//                                                    false, true, false, false,
//                                                    false, false, true, false,
//                                                    false, false, false, false,
//                                                    true, true, false, false,
//                                                    false, false, false, false,
//                                                    true, true, false, false,
//                                                    false, false, false, true,
//                                                    false, false, false, false,
//                                                    false, false, true, false,
//                                                    true, false, false, true,
//                                                    true, false, false, false,
//                                                    false, false, true, true,
//                                                    false, true, false, false,
//                                                    true, false, false, false,
//                                                    true, true, false, false,
//                                                    true, false, false, true,
//                                                    false, false, true, true,
//                                                    true, false, false, false,
//                                                    true, true, false, false,
//                                                    false, false, true, true,
//                                                    true, false, false, false,
//                                                    true, false, false, false,
//                                                    false, false, true, false,
//                                                    false, false, true, false,
//                                                    false, false, false, false,
//                                                    false, true, false, true,
//                                                    false, false, true, false,
//                                                    false, false, false, true,
//                                                    true, false, false, false,
//                                                    true, false, false, true,
//                                                    false, false, true, false,
//                                                    false, true, false, false,
//                                                    false, true, false, false,
//                                                    true, false, true, false,
//                                                    false, false, false, true,
//                                                    false, true, false, false,
//                                                    true, false, false, true,
//                                                    false, false, false, false,
//                                                    true, false, true, false,
//                                                    true, false, false, false,
//                                                    true, true, true, false,
//                                                    false, false, true, false,
//                                                    false, false, false, false,
//                                                    false, false, true, true,
//                                                    false, false, true, true,
//                                                    false, true, false, false,
//                                                    true, false, true, false,
//                                                    false, true, false, true,
//                                                    false, true, false, false,
//                                                    false, false, false, false,
//                                                    false, true, false, false,
//                                                    true, true, false, false,
//                                                    false, false, false, false,
//                                                    true, true, false, false,
//                                                    false, false, true, false,
//                                                    false, true, false, true,
//                                                    false, false, true, false,
//                                                    true, false, false, false,
//                                                    false, false, false, false,
//                                                    true, false, false, false,
//                                                    false, false, true, true,
//                                                    false, true, false, false,
//                                                    false, true, false, false,
//                                                    true, false, true, false,
//                                                    true, false, true, false,
//                                                    true, false, true, false,
//                                                    false, false, true, true,
//                                                    false, false, true, false,
//                                                    false, false, true, true,
//                                                    false, false, true, false,
//                                                    false, true, true, false,
//                                                    false, true, true, false,
//                                                    false, false, false, false,
//                                                    false, true, false, true,
//                                                    true, false, false, false,
//                                                    false, false, false, false,
//                                                    false, true, false, false,
//                                                    false, true, false, true,
//                                                    false, false, true, false,
//                                                    false, false, false, true,
//                                                    false, true, true, false,
//                                                    false, true, false, false,
//                                                    false, false, true, false,
//                                                    true, true, false, false,
//                                                    true, false, false, false,
//                                                    false, false, false, true,
//                                                    false, false, true, false,
//                                                    true, false, true, false,
//                                                    false, true, true, false,
//                                                    false, false, false, false,
//                                                    false, true, false, false,
//                                                    false, true, true, true,
//                                                    true, false, false, false,
//                                                    false, false, false, false,
//                                                    true, false, false, false,
//                                                    true, false, false, true,
//                                                    false, false, false, true,
//                                                    true, false, false, true,
//                                                    false, false, true, true,
//                                                    true, false, false, false,
//                                                    true, false, false, false,
//                                                    false, false, true, false,
//                                                    false, true, false, true,
//                                                    true, false, true, false,
//                                                    false, false, false, false,
//                                                    false, true, false, false,
//                                                    true, false, false, true,
//                                                    false, false, false, true,
//                                                    true, false, true, false,
//                                                    false, false, true, true,
//                                                    false, true, false, true,
//                                                    false, false, false, true,
//                                                    false, true, false, true,
//                                                    false, false, false, false,
//                                                    true, false, true, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false, false, false,
//                                                    false, false};
//
//    std::vector<geometry::Symm<FloatingPointType>> symmvec;
//    symmvec.reserve(symms.size());
//    for (int i = 0; i < symms.size(); i++)
//    {
//        symmvec.emplace_back(symms[i]);
//    }
//
//    cpplib::cluster_detail::UnitCellBuilder ucb(symmvec);
//    auto buildresult = ucb.build(data, std::vector<AtomTypeBase>(data.size(), AtomTypeBase(1)));
//
//    cpplib::geometry::SpatialGrid<FloatingPointType> space;
//    space.build(buildresult.atoms.points, cell, 15.0);
//    auto bonds = space.get_bonds();
//
//    voronoi::VoronoiDiagram vd(buildresult.atoms.points, bonds, cell, bools);
//    auto cells = vd.extractCells();
//}

// ==================== ENTRY POINT ====================
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}