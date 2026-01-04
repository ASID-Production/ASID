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
#include "../src/Classes/Geometry.h"
#include <vector>
#include <array>

using namespace std;
using namespace cpplib;
using namespace cpplib::geometry;

using namespace cpplib::basic_types;

constexpr Point<FloatingPointType> a(0.01423f, 0.27322f, 0.01346f);
constexpr Matrix<FloatingPointType> m({ 10.4804f, -5.2402f, 0.f, 0.f, 9.07629264f, 0.f, 0.f, 0.f, 31.8116f });

TEST(PointTest, OperatorMult_Point_Matrix) {
	Point<FloatingPointType> res;
	ASSERT_NO_THROW({ res = m * a; });
	EXPECT_NEAR(res[0], -1.28259146, 0.00001);
	EXPECT_NEAR(res[1], 2.47982478, 0.00001);
	EXPECT_NEAR(res[2], 0.428184122, 0.00001);
}

TEST(PointTest, CreationNothrow) {
	ASSERT_NO_THROW({
		auto pf = Point<FloatingPointType>();
		auto pd = Point<FloatingPointType>();
					});
	ASSERT_NO_THROW({
		auto pf = Point<FloatingPointType>(0.01423f, 0.27322f, 0.01346f);
		auto pi = Point<FloatingPointType>(0.0, 1.0, 2.0);
					});
}
TEST(PointTest, MemberFunction_r) {
	Point<FloatingPointType>::value_type res;
	ASSERT_NO_THROW({ res = a.r(); });
	EXPECT_NEAR(res, 0.273921, 0.00001);
}
TEST(PointTest, reverseTorsion) {
	Point<FloatingPointType> a1(0.865301423f, 0.21727322f, 0.032461346f);
	Point<FloatingPointType> a2(0.65630135423f, 0.23467322f, 0.835801346f);
	Point<FloatingPointType> a3(0.35601423f, 0.346346227322f, 0.601346f);
	Point<FloatingPointType> a4(0.65301423f, 0.2373334622f, 0.568501346f);
	EXPECT_NEAR(cpplib::geometry::Point<FloatingPointType>::torsionRad(a1, a2, a3, a4), cpplib::geometry::Point<FloatingPointType>::torsionRad(a4, a3, a2, a1), 0.00001);
	EXPECT_NEAR(cpplib::geometry::Point<FloatingPointType>::torsionRad(a2, a1, a3, a4), cpplib::geometry::Point<FloatingPointType>::torsionRad(a4, a3, a1, a2), 0.00001);
}


std::vector<const char*> test_data_symm{
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


TEST(SymmTest, Creation) {
	for (const auto str : test_data_symm) {
		Symm<FloatingPointType> symm(str);
		EXPECT_TRUE(true);
	}
}


TEST(PolygonTest, CreationFromPoints) {
    using PointType = Point<FloatingPointType>;
    std::vector<PointType> points = {
        PointType(0,0,0),
        PointType(1,0,0),
        PointType(1,1,0),
        PointType(0,1,0)
    };
    Polygon<FloatingPointType> poly(points);
    ASSERT_EQ(poly.size(), 4);
    ASSERT_TRUE(poly.isConvex());
}
TEST(PolygonTest, CreationFromPlaneAndCenter) {
    using PointType = Point<FloatingPointType>;
    PointType center(0.5, 0.5, 0);
    Plane<FloatingPointType> plane(center, PointType(0, 0, 1));
    Polygon<FloatingPointType> poly(plane, center, 1.0);
    ASSERT_EQ(poly.size(), 4);
    for (size_t i = 0; i < poly.size(); ++i) {
        ASSERT_NEAR(plane.distance(poly[i]), 0.0, 1e-10);
    }
}
TEST(PolygonTest, ClipByPlane) {
    using PointType = Point<FloatingPointType>;
    std::vector<PointType> points = {
        PointType(0,0,0),
        PointType(1,0,0),
        PointType(1,1,0),
        PointType(0,1,0)
    };
    Polygon<FloatingPointType> poly(points);
    Plane<FloatingPointType> clipping_plane(PointType(0.5, 0, 0), PointType(1, 0, 0));
    poly.clipByPlane(clipping_plane);
    ASSERT_TRUE(poly.isConvex());
    for (size_t i = 0; i < poly.size(); ++i) {
        ASSERT_GE(clipping_plane.side(poly[i]), 0.0);
    }
}
TEST(PolygonTest, PolygonClipPreservesOrder) {
    std::vector<Point<FloatingPointType>> square = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(2, 0, 0),
        Point<FloatingPointType>(2, 2, 0),
        Point<FloatingPointType>(0, 2, 0)
    };

    Polygon<FloatingPointType> poly(square);

    Plane<FloatingPointType> clipping_plane(Point<FloatingPointType>(1, 1, 0), Point<FloatingPointType>(1, 1, 0));
    poly.clipByPlane(clipping_plane);
    EXPECT_EQ(poly.size(), 3);
    EXPECT_TRUE(poly.isConvex());

    for (size_t i = 0; i < poly.size(); ++i) {
        Point<FloatingPointType> current = poly[i];
        Point<FloatingPointType> next = poly[(i + 1) % poly.size()];
        Point<FloatingPointType> next_next = poly[(i + 2) % poly.size()];

        Point<FloatingPointType> edge1 = next - current;
        Point<FloatingPointType> edge2 = next_next - next;
        Point<FloatingPointType> cross = Point<FloatingPointType>::Vector(edge1, edge2);

        EXPECT_GT(cross[2], -1e-10);
    }
}
TEST(PolygonTest, PolygonConvexityCorrectness) {
    std::vector<Point<FloatingPointType>> convex_square = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(1, 0, 0),
        Point<FloatingPointType>(1, 1, 0),
        Point<FloatingPointType>(0, 1, 0)
    };
    Polygon<FloatingPointType> convex_poly(convex_square);
    EXPECT_TRUE(convex_poly.isConvex());

    std::vector<Point<FloatingPointType>> concave_points = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(2, 0, 0),
        Point<FloatingPointType>(1, 1, 0),
        Point<FloatingPointType>(2, 2, 0),
        Point<FloatingPointType>(0, 2, 0)
    };
    Polygon<FloatingPointType> concave_poly(concave_points);
    EXPECT_FALSE(concave_poly.isConvex());

    std::vector<Point<FloatingPointType>> triangle = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(1, 0, 0),
        Point<FloatingPointType>(0, 1, 0)
    };
    Polygon<FloatingPointType> triangle_poly(triangle);
    EXPECT_TRUE(triangle_poly.isConvex());
}


class GeometryTest : public ::testing::Test {
protected:
    void SetUp() override {
        square_points = {
            Point<FloatingPointType>(0, 0, 0),
            Point<FloatingPointType>(1, 0, 0),
            Point<FloatingPointType>(1, 1, 0),
            Point<FloatingPointType>(0, 1, 0)
        };

        triangle_points = {
            Point<FloatingPointType>(0, 0, 0),
            Point<FloatingPointType>(1, 0, 0),
            Point<FloatingPointType>(0, 1, 0)
        };
    }

    std::vector<Point<FloatingPointType>> square_points;
    std::vector<Point<FloatingPointType>> triangle_points;
};

TEST_F(GeometryTest, PolygonCreation) {
    Polygon<FloatingPointType> poly(square_points);
    EXPECT_EQ(poly.size(), 4);
    EXPECT_TRUE(poly.isConvex());
}

TEST_F(GeometryTest, PolygonConvexity) {
    Polygon<FloatingPointType> convex_poly(square_points);
    EXPECT_TRUE(convex_poly.isConvex());

    std::vector<Point<FloatingPointType>> concave_points = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(2, 0, 0),
        Point<FloatingPointType>(1, 1, 0),
        Point<FloatingPointType>(2, 2, 0),
        Point<FloatingPointType>(0, 2, 0)
    };
    Polygon<FloatingPointType> concave_poly(concave_points);
    EXPECT_FALSE(concave_poly.isConvex());
}

TEST_F(GeometryTest, PolygonPlaneCreation) {
    Point<FloatingPointType> center(0.5, 0.5, 0);
    Plane<FloatingPointType> plane(Point<FloatingPointType>(0, 0, 0), Point<FloatingPointType>(1, 0, 0), Point<FloatingPointType>(0, 1, 0));

    Polygon<FloatingPointType> poly(plane, center, 1.0);
    EXPECT_EQ(poly.size(), 4);
    EXPECT_TRUE(poly.isConvex());
}

TEST_F(GeometryTest, PolygonAccessOperator) {
    Polygon<FloatingPointType> poly(square_points);

    EXPECT_EQ(poly[0], Point<FloatingPointType>(0, 0, 0));
    EXPECT_EQ(poly[1], Point<FloatingPointType>(1, 0, 0));
    EXPECT_EQ(poly[2], Point<FloatingPointType>(1, 1, 0));
    EXPECT_EQ(poly[3], Point<FloatingPointType>(0, 1, 0));
}

TEST_F(GeometryTest, PolygonClipByPlane) {
    Polygon<FloatingPointType> poly(square_points);
    Plane<FloatingPointType> clipping_plane(Point<FloatingPointType>(0.5, 0, 0), Point<FloatingPointType>(1, 0, 0));

    poly.clipByPlane(clipping_plane);

    for (size_t i = 0; i < poly.size(); ++i) {
        EXPECT_GE(poly[i][0], 0.5);
    }
    EXPECT_TRUE(poly.isConvex());
}

TEST_F(GeometryTest, PolygonEdgeCases) {
    std::vector<Point<FloatingPointType>> insufficient_points = {
        Point<FloatingPointType>(0, 0, 0),
        Point<FloatingPointType>(1, 0, 0)
    };

    Polygon<FloatingPointType> small_poly(insufficient_points);
    EXPECT_FALSE(small_poly.isConvex());

    Polygon<FloatingPointType> empty_poly;
    EXPECT_EQ(empty_poly.size(), 0);
}

// VoronoiCell
TEST_F(GeometryTest, VoronoiCellDefaultConstructor) {
    VoronoiCell<FloatingPointType> cell;

    Point<FloatingPointType> seed = cell.getSeed();
    EXPECT_EQ(seed, Point<FloatingPointType>(0, 0, 0));

    auto faces = cell.getFaces();
    EXPECT_FALSE(faces.empty());
}

TEST_F(GeometryTest, VoronoiCellSeedConstructor) {
    Point<FloatingPointType> seed(1.0, 2.0, 3.0);
    VoronoiCell<FloatingPointType> cell(seed);

    EXPECT_EQ(cell.getSeed(), seed);

    auto faces = cell.getFaces();
    EXPECT_FALSE(faces.empty());
}

TEST_F(GeometryTest, VoronoiCellInteraction) {
    VoronoiCell<FloatingPointType> cell1(Point<FloatingPointType>(0.1, 0.1, 0.1));
    VoronoiCell<FloatingPointType> cell2(Point<FloatingPointType>(0.4, 0.4, 0.4));

    int result = VoronoiCell<FloatingPointType>::interact(cell1, cell2);

    EXPECT_EQ(result, 0);

    auto faces1 = cell1.getFaces();
    auto faces2 = cell2.getFaces();
    EXPECT_FALSE(faces1.empty());
    EXPECT_FALSE(faces2.empty());
}

TEST_F(GeometryTest, VoronoiCellCloseSeeds) {
    VoronoiCell<FloatingPointType> cell1(Point<FloatingPointType>(0.1, 0.1, 0.1));
    VoronoiCell<FloatingPointType> cell2(Point<FloatingPointType>(0.1000001, 0.1000001, 0.1000001));

    int result = VoronoiCell<FloatingPointType>::interact(cell1, cell2);

    EXPECT_EQ(result, 1);
}

// VoronoiDiagram
TEST_F(GeometryTest, VoronoiDiagramConstructor) {
    std::vector<Point<FloatingPointType>> points = {
        Point<FloatingPointType>(0.1, 0.1, 0.1),
        Point<FloatingPointType>(0.4, 0.4, 0.4),
        Point<FloatingPointType>(0.254, 0.4, 0.364),
        Point<FloatingPointType>(0.954, 0.866, 0.23),
        Point<FloatingPointType>(0.7, 0.7, 0.7)
    };

    VoronoiDiagram<FloatingPointType> vd(points);

    auto cells = vd.extractCells();
    EXPECT_EQ(cells.size(), points.size());


    for (const auto& cell : cells) {
        const auto seed = cell.getSeed();
        for (const auto& face : cell.getFaces()) {
            EXPECT_TRUE(face.isConvex());
            VoronoiDiagram<FloatingPointType>::VoronCell::Face::PlaneType plane(face[0], face[1], face[2]);
            bool res = plane.side(seed) >= 0;
            EXPECT_TRUE(res);
        }
    }
}

TEST_F(GeometryTest, VoronoiDiagramWithFlags) {
    std::vector<Point<double>> points = {
        Point<double>(0.1, 0.1, 0.1),
        Point<double>(0.4, 0.4, 0.4),
        Point<double>(0.254, 0.4, 0.364),
        Point<double>(0.954, 0.866, 0.23),
        Point<double>(0.7, 0.7, 0.7)
    };

    std::vector<bool> flags = { true, false, true , true, true};

    VoronoiDiagram<double> vd(points, flags);

    std::vector<std::pair<int, int>> bondlist = {
        {0,1}, {0,2}, {0,3}, {0,4},
        {1,2}, {1,3}, {1,4},
        {2,3}, {2,4},
        {3,4}
    };
    vd.calculateFaces<int>(bondlist);

    auto cells = vd.extractCells();
    EXPECT_EQ(cells.size(), points.size());

    for (const auto& cell : cells) {
        const auto seed = cell.getSeed();
        for (const auto& face : cell.getFaces()) {
            EXPECT_TRUE(face.isConvex());
            VoronoiDiagram<double>::VoronCell::Face::PlaneType plane(face[0], face[1], face[2]);
            EXPECT_TRUE(plane.side(seed) >= 0);
        }
    }

    VoronoiFused<double> VF;
    VF.AddCells(cells);
    for (size_t i = 0; i < VF.vertexes.size(); i++)
    {
        for (size_t j = i+1; j < VF.vertexes.size(); j++)
        {
            if ((VF.vertexes[i] - VF.vertexes[j]).r() < 0.0001)
                break;
        }
    }
    VF.centers.size();
}

TEST_F(GeometryTest, VoronoiDiagramAddPoints) {
    VoronoiDiagram<FloatingPointType> vd;

    std::vector<Point<FloatingPointType>> points = {
        Point<FloatingPointType>(0.2, 0.2, 0.2),
        Point<FloatingPointType>(0.5, 0.5, 0.5)
    };

    std::vector<bool> flags(points.size(), true);
    vd.addPoints(points, flags);

    auto cells = vd.extractCells();
    EXPECT_EQ(cells.size(), points.size());
}

TEST_F(GeometryTest, VoronoiDiagramCalculateLongestDiagonal) {
    std::vector<Point<FloatingPointType>> points = {
        Point<FloatingPointType>(0.1, 0.1, 0.1),
        Point<FloatingPointType>(0.9, 0.9, 0.9)
    };

    VoronoiDiagram<FloatingPointType> vd(points);

    Cell<FloatingPointType> cell(3, 3, 3, 90, 90, 120, true);
    FloatingPointType diagonal = vd.calculateLongestDiagonal(cell.fracToCart());

    EXPECT_GT(diagonal, 0.0);
}

TEST_F(GeometryTest, VoronoiDiagramExtractCells) {
    std::vector<Point<FloatingPointType>> points = {
        Point<FloatingPointType>(0.1, 0.1, 0.1),
        Point<FloatingPointType>(0.4, 0.4, 0.4)
    };

    VoronoiDiagram<FloatingPointType> vd(points);

    auto cells = vd.extractCells();
    EXPECT_EQ(cells.size(), points.size());

    auto empty_cells = vd.extractCells();
    EXPECT_TRUE(empty_cells.empty());
}

TEST_F(GeometryTest, VoronoiCellBoundaryConditions) {
    VoronoiCell<FloatingPointType> cell1(Point<FloatingPointType>(0.0, 0.0, 0.0));
    VoronoiCell<FloatingPointType> cell2(Point<FloatingPointType>(0.999, 0.999, 0.999));

    int result = VoronoiCell<FloatingPointType>::interact(cell1, cell2);

    EXPECT_EQ(result, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

