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
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <future>
#include "../src/Functions/Functions.h"
#include "../src/Functions/AllInOneAndCurrent.h"
#include "../src/Classes/FindGeometry.h"

#pragma warning( disable : 4305 )

using namespace cpplib::currents;

DistancesType testdistances("BondLength.ini");
struct FMIC_TS {
	std::array<float, 6> cell;
	std::vector<const char*> symm;
	std::vector<int> types;
	std::vector<float> xyz;
	FMIC_TS(std::array<float, 6>&& cell_in, std::vector<const char*>&& symm_in,
					  std::vector<int>&& types_in, std::vector<float>&& xyz_in) :
		cell(std::move(cell_in)),
		symm (std::move(symm_in)),
		types(std::move(types_in)),
		xyz(std::move(xyz_in))
		{}
	inline auto call() {
		return FindMoleculesInCell(cell, symm, types, xyz);
	}
};

TEST(SearchMainTest, Coordf) {
	const char search[] {"1 2 0 6 0 2 4 6 0 2 4 "};
	std::vector<const char*> dat(1, "1 2 1 6 0 6 1 1 2");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(search, std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 0);
}
TEST(SearchMainTest, 107403t) {
	const char search[]{ "1 11 11 17 0 6 0 6 0 6 0 6 0 6 0 7 0 6 0 6 0 6 0 6 0 1 2 2 3 2 10 3 4 3 5 5 6 5 7 7 8 8 9 8 10 10 11" };
	std::vector<const char*> dat(1,"107403 18 20 17 0 6 0 6 1 6 1 6 0 6 1 7 0 6 0 6 0 6 0 6 0 6 1 17 0 6 1 6 1 6 0 8 0 6 3 1 2 2 3 2 4 3 5 4 6 5 7 5 8 6 8 7 9 8 10 9 11 9 12 10 13 10 11 11 14 12 15 14 16 15 16 16 17 17 18");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 1);
}
TEST(SearchMainTest, Tricycle92807t) {
	const char search[] {"1 2 1 6 0 6 0 1 2"};
	std::vector<const char*> dat(1, "92807 35 63 6 1 6 0 17 0 6 0 6 0 17 1 17 1 6 1 17 0 6 0 17 0 6 1 7 2 6 3 17 0 7 1 6 0 6 2 6 1 15 0 6 2 6 0 8 0 6 2 6 1 6 1 6 1 8 0 6 2 7 1 6 0 6 1 6 0 6 1 6 1 1 21 1 31 1 2 1 6 1 17 1 34 2 31 2 29 2 6 2 8 2 34 2 17 3 10 4 23 4 30 4 10 5 33 5 11 5 12 6 34 7 25 7 32 7 22 8 29 8 18 8 34 8 26 9 10 10 15 12 19 13 20 13 14 13 21 14 21 14 31 14 32 14 17 16 20 16 24 17 21 17 31 17 32 17 22 18 29 18 22 18 25 18 26 19 27 20 28 20 30 21 31 22 31 22 25 22 32 22 26 24 33 25 32 25 26 26 29 27 35 29 34 31 32 33 35");
	
	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 1);
}
TEST(SearchMainTest, Tricycle52403f) {
	const char search[] {"1 14 16 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 3 7 4 5 4 8 5 6 7 10 7 11 8 9 9 10 10 14 11 12 12 13 13 14"};
	std::vector<const char*> dat(1,"52403 95 107 6 1 15 0 6 0 6 0 6 1 6 1 6 0 6 3 6 0 6 0 6 1 6 0 6 3 6 0 6 0 6 1 6 0 6 1 6 0 6 3 6 3 6 2 6 0 6 0 6 1 6 0 6 1 6 1 6 0 6 3 6 0 6 2 6 1 6 0 6 3 6 1 6 0 6 1 6 3 6 0 6 1 6 1 6 0 6 1 6 1 6 1 6 1 6 1 6 0 6 0 9 0 9 0 9 0 9 0 5 0 6 0 6 0 6 0 9 0 6 0 6 0 9 0 9 0 6 0 6 0 9 0 9 0 6 0 6 0 9 0 6 0 9 0 6 0 9 0 6 0 6 0 6 0 6 0 6 0 9 0 9 0 6 0 6 0 9 0 6 0 9 0 6 0 6 0 9 0 6 0 9 0 9 0 9 0 6 0 6 0 1 19 1 48 2 14 2 24 3 49 3 23 3 25 4 26 4 14 4 23 5 37 5 31 6 48 6 45 7 42 7 49 7 43 8 12 9 15 9 14 9 50 10 15 10 17 10 27 11 29 11 37 12 26 12 32 12 13 15 22 15 40 16 38 16 43 17 19 17 28 18 49 18 36 19 40 20 31 21 29 22 34 23 32 23 43 24 29 24 31 25 41 26 33 27 46 28 47 30 34 33 50 34 50 34 39 35 37 36 41 38 44 40 45 42 44 46 47 51 57 52 78 53 60 54 90 55 61 55 77 55 65 55 64 56 84 56 85 56 64 57 71 57 95 58 80 58 90 58 75 59 68 60 79 60 76 61 90 61 78 62 75 63 85 64 88 65 79 65 82 66 76 67 94 68 78 68 75 69 81 69 88 69 83 70 83 71 72 71 94 73 86 73 77 73 95 74 87 76 87 77 94 79 91 82 89 82 87 83 85 88 93 92 95");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 0);
}
TEST(SearchMainTest, 85443t) {
	const char search[] {"1 11 11 17 0 6 0 6 0 6 0 6 0 6 0 7 0 6 0 6 0 6 0 6 0 1 2 2 3 2 10 3 4 3 5 5 6 5 7 7 8 8 9 8 10 10 11"};
	std::vector<const char*> dat(1, "85443 19 21 17 0 6 0 6 0 6 0 6 0 6 1 6 0 6 1 7 0 6 1 6 2 6 3 6 1 6 1 6 1 6 2 6 1 6 3 6 3 1 2 2 3 2 4 3 5 3 6 4 7 4 8 5 9 5 10 6 11 6 12 7 9 7 13 8 14 10 15 10 16 11 16 13 17 14 17 15 18 15 19");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 1);
}
TEST(SearchMainTest, Tricycle52245t) {
	const char search[] {"1 14 16 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 3 7 4 5 4 8 5 6 7 10 7 11 8 9 9 10 10 14 11 12 12 13 13 14"};
	std::vector<const char*> dat(1, "52245 24 27 8 1 6 1 6 2 6 1 6 2 6 0 6 0 6 2 6 1 6 3 6 2 6 1 6 2 6 1 6 2 6 0 6 2 6 1 6 3 6 2 6 1 6 3 6 2 8 1 1 2 2 3 2 4 3 5 4 6 5 7 6 8 6 7 7 9 7 10 8 11 9 12 9 13 11 12 12 14 13 15 14 16 14 17 15 16 16 18 16 19 17 20 18 20 18 21 21 22 21 23 23 24");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 1);
}
TEST(SearchMainTest, Tricycle897641f) {
	const char search[] {"1 14 16 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 3 7 4 5 4 8 5 6 7 10 7 11 8 9 9 10 10 14 11 12 12 13 13 14"};
	std::vector<const char*> dat(1, "897641 157 186 48 0 16 0 16 0 16 0 16 0 48 0 48 0 48 0 48 0 6 0 48 0 6 0 48 0 6 0 16 0 16 0 16 0 16 0 16 0 16 0 16 0 16 0 16 0 16 0 6 1 6 1 16 0 6 1 6 1 16 0 6 1 6 1 6 0 6 0 48 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 1 6 1 6 0 6 1 6 1 6 0 6 1 6 1 6 1 6 1 6 1 6 1 17 0 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 16 0 8 0 6 1 16 0 8 0 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 3 6 3 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 44 0 7 0 7 0 7 0 7 0 7 0 7 0 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 0 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 6 1 8 0 8 0 6 0 1 2 1 3 1 4 1 5 13 15 13 5 13 22 13 30 8 2 8 21 8 22 8 23 11 16 11 19 11 4 11 27 7 2 7 18 7 19 7 20 35 56 35 17 35 18 35 23 9 3 9 21 9 20 9 24 6 2 6 15 6 16 6 17 23 42 22 41 20 39 5 14 21 40 4 12 19 38 3 10 18 37 17 36 16 34 34 54 34 55 54 82 82 103 103 83 83 55 15 33 27 46 74 46 30 49 77 49 78 101 78 49 75 100 75 46 10 25 10 26 25 44 44 73 73 45 45 26 14 31 14 32 31 50 50 79 79 51 51 32 12 28 12 29 28 47 47 76 76 48 48 29 39 63 39 64 63 90 90 107 107 91 91 64 37 59 37 60 59 86 86 105 105 87 87 60 41 67 41 68 67 94 94 109 109 95 95 68 38 61 38 62 61 88 88 106 106 89 89 62 40 65 40 66 65 92 92 108 108 93 93 66 42 69 42 70 69 96 96 110 110 97 97 70 36 57 36 58 57 84 84 104 104 85 85 58 33 52 33 53 52 80 80 102 102 81 81 53 24 43 43 71 43 72 71 98 98 111 111 99 99 72 112 113 112 114 112 115 112 116 112 117 112 118 130 118 130 142 142 153 153 141 141 129 141 154 129 118 129 127 154 151 151 139 139 127 139 152 127 117 152 140 140 128 128 117 126 116 126 138 138 149 149 137 137 125 137 150 125 116 125 123 150 147 147 135 135 123 135 148 123 115 148 136 136 124 124 115 122 114 122 134 134 145 145 133 133 121 133 146 121 114 121 119 146 143 143 131 131 119 131 144 119 113 144 132 132 120 120 113 156 157");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 0);
}
TEST(SearchMainTest, NoBonds) {
	const char search[] {"1 1 0 6 0"};
	std::vector<const char*> dat(1, "1 1 0 6 0");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);});
	EXPECT_EQ(res.size(), 1);
}
TEST(SearchMainTest, Tricycle484021t) {
	const char search[] {"1 14 16 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 3 7 4 5 4 8 5 6 7 10 7 11 8 9 9 10 10 14 11 12 12 13 13 14"};
	std::vector<const char*> dat(1, "484021 42 47 8 0 6 0 7 0 6 3 7 0 6 1 6 0 6 2 6 0 6 0 6 0 6 1 6 1 6 1 6 0 6 1 6 1 6 1 6 1 6 0 6 0 6 1 6 1 6 0 6 0 6 1 6 1 6 1 6 1 6 1 6 0 6 0 6 0 9 0 6 0 9 0 6 0 9 0 6 0 9 0 53 0 53 0 1 2 2 3 2 4 3 5 3 6 5 7 6 8 6 9 7 8 7 10 9 11 9 12 10 13 10 14 11 15 11 16 12 17 13 18 14 19 15 20 15 21 16 22 17 21 18 23 19 23 20 24 20 25 21 26 22 24 24 27 25 28 25 29 26 28 27 30 29 30 32 31 32 33 32 34 33 35 33 36 34 37 34 38 36 39 36 40 38 41 38 40 40 42");

	std::vector<int> res;
	ASSERT_NO_THROW({res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false); });
	EXPECT_EQ(res.size(), 1);
}
TEST(CompareGraphTest, Tricycle484021t) {
	const char search[] {"1 14 16 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 2 3 3 4 3 7 4 5 4 8 5 6 7 10 7 11 8 9 9 10 10 14 11 12 12 13 13 14"};
	const char data[] {"484021 42 47 8 0 6 0 7 0 6 3 7 0 6 1 6 0 6 2 6 0 6 0 6 0 6 1 6 1 6 1 6 0 6 1 6 1 6 1 6 1 6 0 6 0 6 1 6 1 6 0 6 0 6 1 6 1 6 1 6 1 6 1 6 0 6 0 6 0 9 0 6 0 9 0 6 0 9 0 6 0 9 0 53 0 53 0 1 2 2 3 2 4 3 5 3 6 5 7 6 8 6 9 7 8 7 10 9 11 9 12 10 13 10 14 11 15 11 16 12 17 13 18 14 19 15 20 15 21 16 22 17 21 18 23 19 23 20 24 20 25 21 26 22 24 24 27 25 28 25 29 26 28 27 30 29 30 32 31 32 33 32 34 33 35 33 36 34 37 34 38 36 39 36 40 38 41 38 40 40 42"};
	bool res = false;
	ASSERT_NO_THROW({res = CompareGraph(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), data, false); });
	EXPECT_TRUE(res);
}
TEST(CompareGraphTest, s161834) {
	const char search[] {"1 36 38 6 1 6 1 6 1 6 1 6 1 6 0 8 0 15 0 8 0 8 0 7 1 6 0 6 0 6 0 8 0 6 0 6 0 6 0 6 0 6 0 6 0 8 0 6 0 6 0 9 0 6 0 8 1 7 0 6 1 6 0 7 1 6 1 6 0 8 0 8 0 8 0 1 2 1 3 2 6 3 4 4 5 5 6 6 7 7 8 8 9 8 10 8 11 10 19 11 12 12 13 12 14 14 15 14 36 15 16 16 17 16 18 19 20 20 21 20 22 21 24 21 27 22 23 23 24 23 28 24 25 24 26 28 29 28 30 29 32 30 31 30 34 31 33 32 33 33 35"};
	const char data[] {"161834 36 38 15 0 8 0 8 0 8 0 7 1 6 2 6 0 6 1 6 1 6 1 6 1 6 3 6 0 8 0 6 1 6 1 6 1 8 0 8 0 6 1 8 1 6 0 6 1 6 1 7 0 9 0 6 3 6 3 6 3 6 0 6 1 8 0 7 1 6 1 6 0 8 0 1 2 1 3 1 4 1 5 3 7 7 10 7 11 10 16 16 23 23 17 17 11 2 6 6 9 9 14 9 15 14 20 20 25 20 22 25 30 25 31 30 32 30 33 33 35 35 36 35 34 34 31 22 26 22 27 22 15 15 21 5 8 8 12 8 13 13 18 13 19 18 24 24 28 24 29"};
	bool res = false;
	ASSERT_NO_THROW({res = CompareGraph(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), data, false); });
	EXPECT_TRUE(res);
}
TEST(FindMoleculesInCellTest, AADRIB) {
	p_distances = &testdistances;
	std::array<float,6> cell { 17.125, 11.453, 8.590, 90, 109.18, 90 };
	std::vector<const char*> symm  { "x,y,z", "-x,1/2+y,-z" };

	std::vector<int> types   { 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<float> xyz   { 0.0601, 0.3274, -0.2849, 0.1674, 0.3811, 0.0296, 0.2086, 0.4422, -0.244, 0.3034, 0.3104, -0.3598, 0.132, 0.1547, -0.2861, -0.0553, 0.2279, -0.3087, 0.1863, 0.2567, 0.238, 0.3156, 0.5543, -0.1163, 0.3753, 0.158, -0.3742, 0.1108, 0.2362, -0.1845, 0.1873, 0.2905, -0.0669, 0.246, 0.3419, -0.1495, 0.2608, 0.2519, -0.2644, 0.1804, 0.2034, -0.3793, -0.0223, 0.309, -0.3479, -0.0614, 0.3988, -0.4668, 0.1713, 0.354, 0.1839, 0.1534, 0.4568, 0.2719, 0.249, 0.5444, -0.2178, 0.2046, 0.6355, -0.3282, 0.3614, 0.258, -0.3995, 0.3981, 0.3291, -0.4998, 0.0733, 0.1923, -0.1348, 0.2228, 0.2307, 0.0088, 0.2946, 0.3598, -0.0664, 0.286, 0.1844, -0.2174, 0.1884, 0.141, -0.4395, 0.1463, 0.2648, -0.4441, -0.0359, 0.4636, -0.4712, -0.1127, 0.4157, -0.454, -0.0735, 0.3718, -0.5755, 0.1159, 0.5044, 0.2119, 0.1931, 0.5033, 0.2955, 0.1415, 0.4303, 0.376, 0.1481, 0.6371, -0.3571, 0.2128, 0.6213, -0.4436, 0.2364, 0.7093, -0.2728, 0.3981, 0.3042, -0.586, 0.4504, 0.3364, -0.4507, 0.374, 0.4104, -0.5151, -0.3808, 0.4875, -0.194, -0.3549, 0.4158, -0.4705, -0.3097, 0.2543, -0.212, -0.1438, 0.2407, -0.0133, -0.2411, 0.5294, -0.0646, -0.3998, 0.6757, -0.1433, -0.2746, 0.4124, -0.6314, -0.3341, 0.1447, -0.4391, -0.0215, 0.3286, 0.082, -0.3045, 0.5212, -0.2168, -0.2827, 0.4356, -0.3301, -0.2468, 0.3201, -0.2498, -0.18, 0.3487, -0.0918, -0.2121, 0.4205, 0.0176, -0.4255, 0.5758, -0.1602, -0.5053, 0.5367, -0.1486, -0.341, 0.4033, -0.6179, -0.4194, 0.3738, -0.7515, -0.3496, 0.1698, -0.3171, -0.4141, 0.1149, -0.2647, -0.0626, 0.2413, 0.0683, -0.0299, 0.1261, 0.1316, -0.3063, 0.5925, -0.2679, -0.2322, 0.4856, -0.3654, -0.2283, 0.2705, -0.3187, -0.1437, 0.4001, -0.0962, -0.1667, 0.451, 0.1056, -0.2583, 0.3862, 0.0412, -0.509, 0.5762, -0.0578, -0.5635, 0.5276, -0.2214, -0.501, 0.4652, -0.1164, -0.4447, 0.4515, -0.8113, -0.4083, 0.3244, -0.8162, -0.4584, 0.3306, -0.7256, -0.4382, 0.1708, -0.2086, -0.3887, 0.0685, -0.1766, -0.4566, 0.0726, -0.3519, 0.0207, 0.1297, 0.2064, -0.0623, 0.082, 0.1704, -0.0166, 0.0772, 0.0486 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
	EXPECT_TRUE(res.find(";") == std::string::npos);
}
TEST(FindMoleculesInCellTest, AAXTHP) {
	p_distances = &testdistances;
	std::array<float, 6> cell { 8.332, 13.644, 16.345, 90, 90, 90 };
	std::vector<const char*> symm { "x,y,z", "1/2+x,1/2-y,-z", "-x,1/2+y,1/2-z", "1/2-x,-y,1/2+z" };

	std::vector<int>  types { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
	std::vector<float> xyz    { 0.0257, 0.107, 0.3146, -0.17, -0.0132, 0.5453, 0.0349, 0.1655, 0.5245, 0.1099, 0.2586, 0.4978, 0.5946, -0.0048, 0.326, 0.6707, -0.0944, 0.2916, 0.1601, 0.1209, 0.0516, 0.0744, 0.1843, -0.0075, 0.1216, 0.0954, 0.3919, 0.2802, 0.0901, 0.3917, 0.3745, 0.1057, 0.3149, 0.2682, 0.1392, 0.2454, 0.3438, 0.1234, 0.1624, -0.1584, 0.2268, 0.2628, -0.1901, 0.3341, 0.2611, -0.0479, -0.0033, 0.4783, -0.07, 0.068, 0.314, -0.196, -0.078, 0.526, -0.129, -0.005, 0.602, -0.253, 0.04, 0.545, 0.2, 0.257, 0.478, 0.068, 0.289, 0.442, 0.104, 0.296, 0.53, 0.77, -0.12, 0.332, 0.604, -0.142, 0.286, 0.707, -0.075, 0.247, 0.001, 0.148, -0.034, 0.152, 0.226, -0.039, 0.01, 0.231, 0.018, 0.34, 0.074, 0.439, 0.451, 0.156, 0.322, 0.247, 0.214, 0.259, 0.346, 0.052, 0.15, 0.452, 0.149, 0.162, -0.243, 0.343, 0.209, -0.102, 0.365, 0.257, -0.245, 0.355, 0.3, 0.0292, 0.0886, 0.4666, 0.1173, 0.0876, 0.2451, -0.0317, 0.2067, 0.3134, -0.2286, 0.1658, 0.2255, -0.0155, -0.0684, 0.4322, -0.0242, 0.1573, 0.5911, 0.4524, 0.0148, 0.2902, 0.651, 0.0451, 0.3791, 0.2537, 0.1753, 0.101, 0.1491, 0.034, 0.0553 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
	EXPECT_TRUE(res.find(";") == std::string::npos);
}
TEST(FindMoleculesInCellTest, AABHTZ) {
	p_distances = &testdistances;
	std::array<float, 6> cell{ 11.372 ,  10.272 ,  7.359 ,  108.75 ,  71.07 ,  96.16 };
	std::vector<const char*> symm  {"x,y,z", "-x,-y,-z" };

	std::vector<int>  types   { 17, 17, 6, 6, 6, 6, 6, 6, 6, 7, 7, 6, 7, 7, 6, 7, 7, 6, 6, 8, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<float> xyz    { -0.3355, 0.0998, 0.1061, -0.6407, -0.3084, 0.327, -0.4775, 0.0388, 0.2307, -0.5727, 0.1337, 0.3424, -0.6845, 0.0916, 0.4483, -0.7021, -0.0445, 0.4447, -0.6069, -0.1387, 0.3295, -0.4904, -0.1007, 0.2181, -0.3844, -0.1938, 0.0886, -0.3668, -0.2977, 0.1351, -0.2629, -0.38, 0.0099, -0.1857, -0.3653, -0.1752, -0.2196, -0.3859, -0.3347, -0.1148, -0.3604, -0.4843, -0.0254, -0.3256, -0.406, -0.0642, -0.3268, -0.2104, 0.0023, -0.2837, -0.0741, -0.0305, -0.1606, 0.0744, 0.0405, -0.1244, 0.2215, -0.1102, -0.0883, 0.0793, -0.2394, -0.486, 0.0718, -0.3307, -0.5093, 0.2548, -0.1454, -0.5542, -0.0296, -0.558, 0.232, 0.347, -0.752, 0.157, 0.531, -0.784, -0.078, 0.53, -0.326, -0.175, -0.032, 0.057, -0.296, -0.469, 0.046, -0.34, -0.057, 0.081, -0.036, 0.217, 0.105, -0.198, 0.189, -0.006, -0.107, 0.333, -0.313, -0.598, 0.275, -0.329, -0.451, 0.374, -0.413, -0.525, 0.243 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
	EXPECT_TRUE(res.find(";") == std::string::npos);
}

TEST(FindMoleculesInCellTest, ABABAH) {
	p_distances = &testdistances;
	std::array<float, 6> cell    { 12.3717, 7.772, 32.821, 90.0, 90.0, 90.0 };
	std::vector<const char*> symm{ "x,y,z", "1/2-x,-y,1/2+z", "-x,1/2+y,1/2-z", "1/2+x,1/2-y,-z", "-x,-y,-z", "-1/2+x,y,-1/2-z", "x,-1/2-y,-1/2+z", "-1/2-x,-1/2+y,z" };
	std::vector<int>  types      { 8, 8, 8, 6, 6, 1, 6, 6, 1, 6, 1, 1, 1, 6, 1, 1, 1, 6, 1, 1, 6, 1, 6, 1, 1, 1, 6, 1, 1, 1, 6, 6, 1, 6, 1, 6, 6, 1, 6, 1, 6, 1, 1, 1 };
	std::vector<float> xyz       { 0.23473, 0.3421, 0.41042, 0.30561, 0.2589, 0.4695, 0.12786, 0.2198, 0.22618, 0.22729, 0.3035, 0.4506, 0.11298, 0.3248, 0.4643, 0.1042, 0.4477, 0.4726, 0.04731, 0.3, 0.42445, 0.12831, 0.385, 0.39434, 0.1191, 0.5126, 0.3963, 0.0865, 0.2155, 0.50145, 0.0096, 0.2267, 0.5079, 0.1298, 0.2546, 0.5247, 0.1033, 0.0947, 0.4956, -0.05842, 0.4008, 0.42562, -0.1089, 0.3447, 0.4444, -0.0901, 0.4041, 0.3983, -0.0442, 0.5185, 0.4349, 0.03096, 0.1075, 0.41599, 0.004, 0.0542, 0.4414, 0.103, 0.0571, 0.4105, -0.0451, 0.0508, 0.38116, -0.0428, 0.1396, 0.3591, -0.1618, 0.0309, 0.39557, -0.2065, -0.0117, 0.3731, -0.1894, 0.1427, 0.4047, -0.1646, -0.0511, 0.4182, -0.0053, -0.1192, 0.36384, -0.0519, -0.154, 0.3411, -0.0075, -0.2075, 0.3851, 0.069, -0.1058, 0.3541, 0.12725, 0.3367, 0.34992, 0.06433, 0.4256, 0.32196, 0.0202, 0.518, 0.331, 0.06475, 0.382, 0.28113, 0.0204, 0.4433, 0.2625, 0.12955, 0.2495, 0.26744, 0.1927, 0.1605, 0.2946, 0.2376, 0.0694, 0.2854, 0.1907, 0.204, 0.33554, 0.2341, 0.1408, 0.3541, 0.2029, 0.0987, 0.21071, 0.192, 0.0854, 0.1813, 0.1919, -0.0124, 0.2242, 0.2766, 0.1393, 0.2159 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
	EXPECT_TRUE(res.find(";") == std::string::npos);
}

TEST(FindMoleculesInCellTest, AZIVIO) {
	p_distances = &testdistances;

	std::array<float, 6> cell    { 20.479, 7.9523, 22.452, 90.0, 114.108, 90.0 };
	std::vector<const char*> symm{ "x,y,z", "x,-y,1/2+z", "1/2+x,1/2+y,z", "1/2+x,1/2-y,1/2+z", "-x,-y,-z", "-x,y,-1/2-z", "-1/2-x,-1/2-y,-z", "-1/2-x,-1/2+y,-1/2-z" };
	std::vector<int>  types      { 46, 46, 17, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 46, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 17, 1, 1 };
	std::vector<float> xyz       { 0.56986, 0.527, 0.04919, 0.4273, 0.5362, -0.0289, 0.913, 0.8558, 0.2095, 0.9391, 0.4981, 0.1186, 1.0269, 0.4049, 0.0946, 0.9633, 0.225, 0.1281, 0.9148, 0.343, 0.0291, 0.6274, 0.679, 0.0048, 0.5515, 0.717, -0.0252, 0.558, 0.613, -0.056, 0.5011, 0.361, 0.0813, 0.491, 0.305, 0.057, 0.5726, 0.3099, 0.1069, 0.6289, 0.403, 0.1496, 0.625, 0.271, 0.103, 0.7047, 0.356, 0.1712, 0.7315, 0.254, 0.1378, 0.8028, 0.217, 0.1615, 0.8487, 0.283, 0.2201, 0.8232, 0.387, 0.2542, 0.751, 0.429, 0.2286, 0.6632, 0.557, -0.0169, 0.7349, 0.525, 0.0201, 0.7706, 0.607, 0.0797, 0.7349, 0.718, 0.1023, 0.6642, 0.754, 0.0658, 0.9637, 0.677, 0.2503, 0.9614, 0.364, 0.0925, 0.5313, 0.8163, -0.0126, 0.4848, 0.4431, 0.1061, 0.5835, 0.201, 0.0937, 0.6215, 0.4716, 0.1843, 0.7, 0.2072, 0.0972, 0.8206, 0.1465, 0.1374, 0.8983, 0.2547, 0.2374, 0.85, 0.4258, 0.2972, 0.7373, 0.514, 0.2525, 0.6379, 0.4978, -0.0566, 0.7599, 0.4472, 0.005, 0.8217, 0.5521, 0.0989, 0.7629, 0.7628, 0.1488, 0.6399, 0.8316, 0.0821, 0.938, 0.577, 0.238, 0.969, 0.678, 0.295, 0.43014, 0.473, -0.04919, 0.3726, 0.321, -0.0048, 0.4485, 0.283, 0.0252, 0.4989, 0.639, -0.0813, 0.4687, 0.1837, 0.0126, 0.4274, 0.6901, -0.1069, 0.3711, 0.597, -0.1496, 0.4165, 0.799, -0.0937, 0.2953, 0.644, -0.1712, 0.2685, 0.746, -0.1378, 0.1972, 0.783, -0.1615, 0.1513, 0.717, -0.2201, 0.1768, 0.613, -0.2542, 0.249, 0.571, -0.2286, 0.2627, 0.486, -0.2525, 0.15, 0.5742, -0.2972, 0.1017, 0.7453, -0.2374, 0.1794, 0.8535, -0.1374, 0.3, 0.7928, -0.0972, 0.3785, 0.5284, -0.1843, 0.3358, 0.246, -0.0658, 0.3368, 0.443, 0.0169, 0.2651, 0.475, -0.0201, 0.2294, 0.393, -0.0797, 0.2651, 0.282, -0.1023, 0.3601, 0.1684, -0.0821, 0.2371, 0.2372, -0.1488, 0.1783, 0.4479, -0.0989, 0.2401, 0.5528, -0.005, 0.3621, 0.5022, 0.0566, 0.5152, 0.5569, -0.1061, 1.0363, 0.677, 0.2497, 1.087, 0.8558, 0.2905, 1.062, 0.577, 0.262, 1.031, 0.678, 0.205 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
}

TEST(FindMoleculesInCellTest, EKESIW) {
	p_distances = &testdistances;

	std::array<float, 6> cell    { 6.83, 11.339, 14.627, 94.58, 95.92, 102.54 };
	std::vector<const char*> symm{ "x,y,z", "-x,-y,-z" };
	std::vector<int>  types      { 29, 19, 19, 7, 6, 6, 7, 6, 1, 6, 1, 6, 6, 8, 8, 8, 8, 7, 6, 6, 7, 6, 1, 6, 1, 6, 6, 8, 8, 8, 8, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 8, 8, 8, 8, 8, 8, 19, 19, 19, 19, 19, 19, 1, 1, 1, 1 };
	std::vector<float> xyz       { -0.27059, -0.2371, 0.00097, -0.28596, -0.08164, 0.57449, -0.21667, -0.40913, 0.42584, -0.2366, -0.1124, 0.10899, -0.2214, -0.1596, 0.1874, -0.2208, -0.0927, 0.2724, -0.2343, 0.0234, 0.27433, -0.2488, 0.0702, 0.1945, -0.257, 0.151, 0.195, -0.2521, 0.0034, 0.111, -0.265, 0.039, 0.056, -0.2054, -0.2902, 0.1788, -0.2086, -0.1443, 0.3643, -0.2252, -0.3425, 0.0965, -0.1747, -0.3422, 0.24852, -0.3715, -0.1982, 0.38902, -0.0396, -0.1253, 0.40944, -0.2779, -0.3595, -0.10539, -0.2755, -0.3134, -0.1851, -0.2845, -0.385, -0.2682, -0.2959, -0.5048, -0.26812, -0.2984, -0.5502, -0.1871, -0.307, -0.633, -0.186, -0.2887, -0.4787, -0.1045, -0.29, -0.513, -0.049, -0.2628, -0.1791, -0.1773, -0.2863, -0.3395, -0.362, -0.2458, -0.1264, -0.09497, -0.2757, -0.1259, -0.24745, -0.1217, -0.2947, -0.3886, -0.4548, -0.3553, -0.40796, -0.6, -0.262, -0.0011, -0.665, -0.217, 0.026, -0.682, -0.301, -0.047, -0.298, 0.12, 0.45529, -0.188, 0.161, 0.485, -0.28, 0.123, 0.399, -0.2812, 0.162, 0.6992, -0.259, 0.196, 0.754, -0.372, 0.193, 0.672, -0.1987, -0.6086, 0.5542, -0.214, -0.604, 0.611, -0.303, -0.659, 0.527, -0.2141, -0.6414, 0.3117, -0.12, -0.669, 0.339, -0.32, -0.671, 0.336, -0.7594, -0.1236, 0.1177, -0.736, -0.047, 0.13, -0.886, -0.145, 0.1, -0.1772, -0.5861, 0.1378, -0.181, -0.604, 0.193, -0.196, -0.514, 0.139, 0.0396, 0.1253, 0.59056, -0.2757, -0.1259, 0.75255, -0.1217, -0.2947, 0.6114, -0.4548, -0.3553, 0.59204, -0.702, -0.12, 0.54471, -0.5452, -0.6447, 0.40796, 0.1987, -0.3914, 0.4458, 0.28596, 0.08164, 0.42551, -0.28596, -0.08164, -0.42551, -0.21667, -0.40913, -0.57416, -0.78333, -0.59087, -0.42584, -0.71404, 0.08164, 0.42551, 0.21667, -0.59087, 0.57416, -0.812, -0.161, 0.515, -0.72, -0.123, 0.601, 0.214, -0.396, 0.389, 0.303, -0.341, 0.473 };
	std::string res;

	ASSERT_NO_THROW({ res = FindMoleculesInCell(cell, symm, types, xyz).first; });
	EXPECT_TRUE(res.find(";") == std::string::npos);
}

TEST(FindMoleculesWithoutCellTest, C3H6O2) {
	p_distances = &testdistances;
	std::vector<int>  types {1,8,8,6,6,1,1,6,1,1,1};

	std::vector<float>  xyz { 0.438335628 ,     0.000000000  ,   -4.767885621 ,
	0.438335628 ,    -0.587240971  ,   -3.993208195 ,
	0.438335628 ,     1.434717404  ,   -3.001401644 ,
	0.438335628 ,     0.227304724  ,   -2.907173504 ,
	0.438335628 ,    -0.566448428  ,   -1.619029840 ,
	1.310639628 ,    -1.231202352  ,   -1.642431518 ,
   -0.433968372 ,    -1.231202352  ,   -1.642431518 ,
	0.438335628 ,     0.318369069  ,   -0.375358357 ,
	0.438335628 ,    -0.294739868  ,    0.530163989 ,
	1.319085628 ,     0.965221387  ,   -0.354947074 ,
   -0.442414372 ,     0.965221387  ,   -0.354947074 };
	auto res = FindMoleculesWithoutCell(types, xyz);
}
TEST(FindMoleculesWithoutCellTest, MonoO) {
	p_distances = &testdistances;
	std::vector<int>  types{ 6 };

	std::vector<float>  xyz{ 0.438335628 , -0.587240971, -3.993208195 };
	auto res = FindMoleculesWithoutCell(types, xyz);
	ASSERT_FALSE(res.empty());
}

TEST(FindDistanceTest, C3H6O2) {
	p_distances = &testdistances;
	std::vector<int>types { 1,8,8,6,6,1,1,6,1,1,1 };

	std::vector<float> xyz { 0.438335628 ,     0.000000000  ,   -4.767885621 ,
	0.438335628 ,    -0.587240971  ,   -3.993208195 ,
	0.438335628 ,     1.434717404  ,   -3.001401644 ,
	0.438335628 ,     0.227304724  ,   -2.907173504 ,
	0.438335628 ,    -0.566448428  ,   -1.619029840 ,
	1.310639628 ,    -1.231202352  ,   -1.642431518 ,
   -0.433968372 ,    -1.231202352  ,   -1.642431518 ,
	0.438335628 ,     0.318369069  ,   -0.375358357 ,
	0.438335628 ,    -0.294739868  ,    0.530163989 ,
	1.319085628 ,     0.965221387  ,   -0.354947074 ,
   -0.442414372 ,     0.965221387  ,   -0.354947074 };
	auto res = FindDistanceWC(types, xyz, { 8,8 }, { 1, 5 });
	ASSERT_FALSE(res.empty());
}

TEST(FindAngleTest, C3H6O2) {
	p_distances = &testdistances;
	std::vector<int>types{ 1,8,8,6,6,1,1,6,1,1,1 };

	std::vector<float> xyz{ 0.438335628 ,     0.000000000  ,   -4.767885621 ,
                        	0.438335628 ,    -0.587240971  ,   -3.993208195 ,
                        	0.438335628 ,     1.434717404  ,   -3.001401644 ,
                        	0.438335628 ,     0.227304724  ,   -2.907173504 ,
                        	0.438335628 ,    -0.566448428  ,   -1.619029840 ,
                        	1.310639628 ,    -1.231202352  ,   -1.642431518 ,
                           -0.433968372 ,    -1.231202352  ,   -1.642431518 ,
                        	0.438335628 ,     0.318369069  ,   -0.375358357 ,
                        	0.438335628 ,    -0.294739868  ,    0.530163989 ,
                        	1.319085628 ,     0.965221387  ,   -0.354947074 ,
                           -0.442414372 ,     0.965221387  ,   -0.354947074 };
	auto res = FindAngleWC(types, xyz, { 8,6,8 }, { std::make_pair(1.0f, 2.0f), std::make_pair(0.0f, 2.0f), }, std::make_pair(-180.f, 180.f));
	ASSERT_FALSE(res.empty());
}

TEST(FindTorsionTest, C3H6O2) {
	p_distances = &testdistances;
	std::vector<int>types{ 1,8,8,6,6,1,1,6,1,1,1 };

	std::vector<float> xyz{ 0.438335628 ,     0.000000000  ,   -4.767885621 ,
							0.438335628 ,    -0.587240971  ,   -3.993208195 ,
							0.438335628 ,     1.434717404  ,   -3.001401644 ,
							0.438335628 ,     0.227304724  ,   -2.907173504 ,
							0.438335628 ,    -0.566448428  ,   -1.619029840 ,
							1.310639628 ,    -1.231202352  ,   -1.642431518 ,
						   -0.433968372 ,    -1.231202352  ,   -1.642431518 ,
							0.438335628 ,     0.318369069  ,   -0.375358357 ,
							0.438335628 ,    -0.294739868  ,    0.530163989 ,
							1.319085628 ,     0.965221387  ,   -0.354947074 ,
						   -0.442414372 ,     0.965221387  ,   -0.354947074 };
	auto res = FindTorsionWC(types, xyz, { 8,6,6,6 }, 
							 { std::make_pair(1.0f, 2.0f), std::make_pair(0.0f, 2.0f), std::make_pair(0.0f, 2.0f) }, 
							 { std::make_pair(0.f, 180.f), std::make_pair(0.f, 180.f) },
							 std::make_pair(-180.f, 180.f));
	ASSERT_FALSE(res.empty());
}
TEST(FindDistanceTest, EKESIW) {
	p_distances = &testdistances;

	std::array<float, 6> cell{ 6.83, 11.339, 14.627, 94.58, 95.92, 102.54 };
	std::vector<const char*> symm{ "x,y,z", "-x,-y,-z" };
	std::vector<int>  types{ 29, 19, 19, 7, 6, 6, 7, 6, 1, 6, 1, 6, 6, 8, 8, 8, 8, 7, 6, 6, 7, 6, 1, 6, 1, 6, 6, 8, 8, 8, 8, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 8, 8, 8, 8, 8, 8, 19, 19, 19, 19, 19, 19, 1, 1, 1, 1 };
	std::vector<float> xyz{ -0.27059, -0.2371, 0.00097, -0.28596, -0.08164, 0.57449, -0.21667, -0.40913, 0.42584, -0.2366, -0.1124, 0.10899, -0.2214, -0.1596, 0.1874, -0.2208, -0.0927, 0.2724, -0.2343, 0.0234, 0.27433, -0.2488, 0.0702, 0.1945, -0.257, 0.151, 0.195, -0.2521, 0.0034, 0.111, -0.265, 0.039, 0.056, -0.2054, -0.2902, 0.1788, -0.2086, -0.1443, 0.3643, -0.2252, -0.3425, 0.0965, -0.1747, -0.3422, 0.24852, -0.3715, -0.1982, 0.38902, -0.0396, -0.1253, 0.40944, -0.2779, -0.3595, -0.10539, -0.2755, -0.3134, -0.1851, -0.2845, -0.385, -0.2682, -0.2959, -0.5048, -0.26812, -0.2984, -0.5502, -0.1871, -0.307, -0.633, -0.186, -0.2887, -0.4787, -0.1045, -0.29, -0.513, -0.049, -0.2628, -0.1791, -0.1773, -0.2863, -0.3395, -0.362, -0.2458, -0.1264, -0.09497, -0.2757, -0.1259, -0.24745, -0.1217, -0.2947, -0.3886, -0.4548, -0.3553, -0.40796, -0.6, -0.262, -0.0011, -0.665, -0.217, 0.026, -0.682, -0.301, -0.047, -0.298, 0.12, 0.45529, -0.188, 0.161, 0.485, -0.28, 0.123, 0.399, -0.2812, 0.162, 0.6992, -0.259, 0.196, 0.754, -0.372, 0.193, 0.672, -0.1987, -0.6086, 0.5542, -0.214, -0.604, 0.611, -0.303, -0.659, 0.527, -0.2141, -0.6414, 0.3117, -0.12, -0.669, 0.339, -0.32, -0.671, 0.336, -0.7594, -0.1236, 0.1177, -0.736, -0.047, 0.13, -0.886, -0.145, 0.1, -0.1772, -0.5861, 0.1378, -0.181, -0.604, 0.193, -0.196, -0.514, 0.139, 0.0396, 0.1253, 0.59056, -0.2757, -0.1259, 0.75255, -0.1217, -0.2947, 0.6114, -0.4548, -0.3553, 0.59204, -0.702, -0.12, 0.54471, -0.5452, -0.6447, 0.40796, 0.1987, -0.3914, 0.4458, 0.28596, 0.08164, 0.42551, -0.28596, -0.08164, -0.42551, -0.21667, -0.40913, -0.57416, -0.78333, -0.59087, -0.42584, -0.71404, 0.08164, 0.42551, 0.21667, -0.59087, 0.57416, -0.812, -0.161, 0.515, -0.72, -0.123, 0.601, 0.214, -0.396, 0.389, 0.303, -0.341, 0.473 };

	auto res = FindDistanceIC(cell,symm, types, xyz, { 8, 6 }, { 1, 2 });
	EXPECT_TRUE(true);
}

TEST(GenSymmTest, E28) {
	p_distances = &testdistances;
	std::vector<std::tuple<int,float,float,float>> arg = {{16, 0.55638, 0.01587, 0.36247}, {6, 0.3188, -0.0361, 0.9365}, {8, 0.482, 0.22765, 0.8132}, {7, 0.3478, 0.0758, 0.9261}, {8, 0.9234, 0.03692, 0.6186}, {6, 0.4816, 0.1263, 0.8315}, {7, 0.6105, 0.05617, 0.7551}, {6, 0.457, -0.1057, 0.8582}, {1, 0.4495, -0.1831, 0.8662}, {7, 0.1671, -0.074, 1.018}, {8, 0.9404, -0.1496, 0.2756}, {6, 0.8764, 0.0011, 0.448}, {1, 0.9622, 0.0507, 0.3884}, {6, 0.5993, -0.0573, 0.7721}, {1, 0.6937, -0.1024, 0.7223}, {6, 0.5781, 0.1417, 0.4813}, {1, 0.6446, 0.2021, 0.4287}, {1, 0.4209, 0.164, 0.4952}, {6, 0.7468, 0.1111, 0.6457}, {1, 0.8252, 0.1788, 0.6996}, {6, 0.9693, -0.1166, 0.4419}, {1, 1.1389, -0.1199, 0.4998}, {1, 0.8807, -0.1671, 0.4964}, {1, 0.092, -0.026, 1.064}, {1, 0.257, 0.119, 0.97}, {1, 0.773, -0.176, 0.22}, {1, 0.133, -0.153, 1.019}, {17, 0.02131, 0.18155, 0.12771}};
	std::vector<PointType> points;
	std::vector<AtomTypeData> atoms;
	for (size_t i = 0; i < arg.size(); i++)
	{
		atoms.emplace_back(std::get<0>(arg[i]));
		points.emplace_back(std::get<1>(arg[i]), std::get<2>(arg[i]), std::get<3>(arg[i]));
	}
	bool mvtc = false;

	std::vector<const char*> symms {"x+1,y,z","x-1,y,z", "x,y+1,z" };
	std::vector<SymmType> symm1; symm1.emplace_back(symms[0]);
	std::vector<SymmType> symm2; symm2.emplace_back(symms[1]);
	std::vector<SymmType> symm3; symm3.emplace_back(symms[2]);

	FAMStructType famstr(std::move(atoms), std::move(points));
	FAMCellType fcell(CellType(10, 10, 10, 90, 90, 90, true));
	fcell.GenerateSymm(famstr, symm1, mvtc);
	std::cout << "After x+1: " << famstr.sizePoints << std::endl;
	fcell.GenerateSymm(famstr, symm2, mvtc);
	std::cout << "After x-1: " << famstr.sizePoints << std::endl;
	fcell.GenerateSymm(famstr, symm3, mvtc);
	std::cout << "After y+1: " << famstr.sizePoints << std::endl;

	EXPECT_TRUE(true);
}

TEST(SearchMainTest, multytype1) {
	const char search[] {"1 7 7 6 0 6 0 6 0 6 0 6 0 6 0 -1 0 1 2 1 6 1 7 2 3 3 4 4 5 5 6 -1 6 7 8 0"};
	std::vector<const char*> dat ={ "1 7 7 6 0 6 0 6 0 6 0 6 0 6 0 6 0 1 2 1 6 1 7 2 3 3 4 4 5 5 6"
								 , "2 7 7 6 0 6 0 6 0 6 0 6 0 6 0 7 0 1 2 1 6 1 7 2 3 3 4 4 5 5 6"
								 , "3 7 7 6 0 6 0 6 0 6 0 6 0 6 0 5 0 1 2 1 6 1 7 2 3 3 4 4 5 5 6"
};

	std::vector<int> res;
	res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);
	EXPECT_EQ(res.size(), 2);
}
TEST(SearchMainTest, multytype2) {
	const char search[]{ "1 3 2 6 3 6 2 -1 0 1 2 2 3 -1 1 9 0" };
	std::vector<const char*> dat = { "1 2 1 6 3 6 3 1 2",
									 "2 3 2 6 3 6 2 1 0 1 2 2 3",
									 "3 3 2 9 0 6 2 6 3 1 2 2 3" };

	std::vector<int> res;
	res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);
	EXPECT_EQ(res.size(), 3);
}
TEST(SearchMainTest, multytype3) {
	const char search[] {"1 9 9 6 0 6 0 6 0 6 0 6 0 6 0 -1 0 6 0 6 3 1 2 1 6 1 7 2 3 3 4 4 5 5 6 7 8 8 9 -1 8 16 32 0"};
	std::vector<const char*> dat = {"1 9 9 6 0 6 0 6 0 6 0 6 0 6 0 8 0 6 0 6 3 1 2 1 6 1 7 2 3 3 4 4 5 5 6 7 8 8 9",
							        "2 9 9 6 0 6 0 6 0 6 0 6 0 6 0 16 0 6 0 6 3 1 2 1 6 1 7 2 3 3 4 4 5 5 6 7 8 8 9"};

	std::vector<int> res;
	res = SearchMain(cpplib::MoleculeGraph<cpplib::currents::AtomTypeRequest>::_ParseOldInputString(search).data(), std::move(dat), 1, false);
	EXPECT_EQ(res.size(), 2);
}


TEST(databaseSearch10k, d10k) {
	const char search1[]{ "1 5 4 6 2 6 2 6 0 6 2 6 2 1 2 2 3 3 4 4 5" };
	const char search2[]{ "1 5 4 6 2 6 2 7 0 6 2 6 2 1 2 2 3 3 4 4 5" };
	const char search3[]{ "1 5 4 6 2 6 2 8 0 6 2 6 2 1 2 2 3 3 4 4 5" };
	const char search4[]{ "1 5 4 6 2 6 2 -1 0 6 2 6 2 1 2 2 3 3 4 4 5 -1 6 7 8 0" };
	std::ifstream db("../../../../../cpplib/Tests/d10k.datt");
	if (!db.is_open())
		FAIL() << "Could not open the dataset";
	int s;
	db >> s;
	std::vector<std::string> datstr(s);
	std::vector<const char*> dat(s);
	std::vector<const char*> temp;
	std::getline(db, datstr[0]);
	for (size_t i = 0; i < s; i++)
	{
		std::getline(db, datstr[i]);
		dat[i] = datstr[i].c_str();
	}

	temp = dat;
	std::vector<int> res1 = SearchMain(search1, std::move(temp), 4, false);
	std::cout << "res1 = " << res1.size() << std::endl;
	temp = dat;
	std::vector<int> res2 = SearchMain(search2, std::move(temp), 4, false);
	std::cout << "res2 = " << res2.size() << std::endl;
	temp = dat;
	std::vector<int> res3 = SearchMain(search3, std::move(temp), 4, false);
	std::cout << "res3 = " << res3.size() << std::endl;
	temp = dat;
	std::vector<int> res4 = SearchMain(search4, std::move(temp), 4, false);
	std::cout << "res4 = " << res4.size() << std::endl;
	EXPECT_GE(res1.size() + res2.size() + res3.size(), res4.size());

	std::sort(res1.begin(), res1.end());
	std::sort(res2.begin(), res2.end());
	std::sort(res3.begin(), res3.end());
	std::sort(res4.begin(), res4.end());

	size_t i1 = 0;
	size_t i2 = 0;
	size_t i3 = 0;

	for (size_t i = 0; i < res4.size(); i++)
	{
		bool flag = false;
		for (; i1 < res1.size(); i1++)
		{
			if (res1[i1] == res4[i]) {
				flag = true;
				i1++;
				break;
			}
			if (res1[i1] > res4[i]) {
				break;
			}
		}
		for (; i2 < res2.size(); i2++)
		{
			if (res2[i2] == res4[i]) {
				flag = true;
				i2++;
				break;
			}
			if (res2[i2] > res4[i]) {
				break;
			}
		}
		for (; i3 < res3.size(); i3++)
		{
			if (res3[i3] == res4[i]) {
				flag = true;
				i3++;
				break;
			}
			if (res3[i3] > res4[i]) {
				break;
			}
		}
		EXPECT_TRUE(flag) << res4[i] << std::endl;
	}
	EXPECT_TRUE(res1.size() == i1);
	EXPECT_TRUE(res2.size() == i2);
	EXPECT_TRUE(res3.size() == i3);
}