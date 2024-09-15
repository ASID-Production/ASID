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
#include "../src/Classes/Engine.h"
//#include <vector>

using namespace cpplib;
using namespace cpplib::currents;

// XAtom section
TEST(EngineTest, XAtomType) {
	ASSERT_NO_THROW({ XAtom atom; });
	ASSERT_NO_THROW({
		XAtom::SimpleAtomType a(6);
		XAtom atom(a);
	});

	XAtom atom(-1);
	atom.AddType(1);
	atom.AddType(9);
	auto bits = atom.get_bitset();

	decltype(bits) bits2;
	bits2.set(1);
	bits2.set(9);

	EXPECT_EQ(bits, bits2);

	auto simple = atom.get_simple();
	auto simple2 = static_cast<char>(atom.get_simple());
	auto simple3 = static_cast<AtomTypeData>(atom);

	EXPECT_EQ(simple, -1);
	EXPECT_EQ(simple2, simple3);
	EXPECT_TRUE(atom.simple_eq(-1));
	EXPECT_TRUE(atom == XAtom::SimpleAtomType(9));
	EXPECT_TRUE(atom == XAtom::SimpleAtomType(1));
	EXPECT_FALSE(atom == XAtom::SimpleAtomType(6));
}

// Coord section
TEST(EngineTest, CoordType) {
	ASSERT_NO_THROW({ Coord c; });
	ASSERT_NO_THROW({ Coord c(0); });
	ASSERT_NO_THROW({ Coord c(1,2); Coord c2(0,14); });
	Coord mono(3);
	Coord duo1(1, 2);
	Coord duo2(2, 3);
	Coord duo3(2, 4);
	Coord duo4(3, 4);
	Coord duo5(4, 4);
	EXPECT_FALSE(duo1.right_in_left(mono));
	EXPECT_TRUE(duo2.right_in_left(mono));
	EXPECT_TRUE(duo3.right_in_left(mono));
	EXPECT_TRUE(duo4.right_in_left(mono));
	EXPECT_FALSE(duo5.right_in_left(mono));
}


// NeighboursType section
TEST(EngineTest, NeighboursType) {
	ASSERT_NO_THROW({ NeighboursType a; });
	ASSERT_NO_THROW({ NeighboursType::ShiftType shift = 1; });
	NeighboursType a;
	NeighboursType::ShiftType shift = 1;
	EXPECT_EQ(a.size(), 0);
	a.push_back(shift);
	EXPECT_EQ(a.size(), 1);
	a.push_back(-shift);
	EXPECT_EQ(a.size(), 2);
	EXPECT_EQ(a[0], shift);
	EXPECT_EQ(a[1], -shift);

	ASSERT_NO_THROW({ a.simpleSort(); });
	EXPECT_EQ(a.size(), 2);
	EXPECT_EQ(a[0], -shift);
	EXPECT_EQ(a[1], shift);
	ASSERT_NO_THROW({ a.erase(0); });
	EXPECT_EQ(a[0], shift);
	EXPECT_EQ(a.size(), 1);
	EXPECT_TRUE(a.exchange(shift, -shift));
	EXPECT_FALSE(a.exchange(0, shift));
	EXPECT_EQ(a[0], -shift);
	EXPECT_EQ(a.size(), 1);
	ASSERT_NO_THROW({ a.addShift(shift); });
	ASSERT_NO_THROW({ a.addShift(-shift); });
	EXPECT_EQ(a[0], -shift);
	EXPECT_EQ(a[1], shift);
	EXPECT_EQ(a[2], -shift);
	EXPECT_EQ(a.size(), 3);

	// NeighboursType::sort<A> was not tested
}


// Node Section
using NodeData = cpplib::Node<AtomTypeData>;
using NodeRequest = cpplib::Node<AtomTypeRequest>;
TEST(NodeRequestTest, Constructors) {

	EXPECT_NO_THROW({ NodeRequest a1(1, 1, 1); });
	NodeRequest a1(1, 1, 1);
	//NodeRequest::NeighboursType v1;
	//v1.push_back(&a1);
	EXPECT_NO_THROW({ NodeRequest a2(1, 1, 1); });
	EXPECT_NO_THROW({ NodeRequest a3; });
	//EXPECT_NO_THROW({NodeRequest a6(1,1,1,std::move(v1)); });
}
TEST(NodeDataTest, Constructors) {

	EXPECT_NO_THROW({ NodeData a1(1, 1, 1); });
	NodeData a1(1, 1, 1);
	//NodeData::NeighboursType v1;
	//v1.push_back(&a1);
	EXPECT_NO_THROW({ NodeData a2(1, 1, 1); });
	EXPECT_NO_THROW({ NodeData a3; });
	//EXPECT_NO_THROW({NodeData a6(1,1,1,std::move(v1)); });
}
TEST(NodeRequestTest, EqualWithoutNeighbours) {
	std::vector<NodeRequest> v1;
	v1.emplace_back(1, 0, 1);
	v1.emplace_back(2, 0, 2);
	v1.emplace_back(1, 1, 3);
	v1.emplace_back(1, 0, 4);

	EXPECT_FALSE(v1[0] == v1[1]);
	EXPECT_FALSE(v1[0] == v1[2]);
	EXPECT_TRUE(v1[0] == v1[3]);
	EXPECT_FALSE(v1[1] == v1[2]);
	EXPECT_FALSE(v1[1] == v1[3]);
	EXPECT_FALSE(v1[2] == v1[3]);
}
TEST(NodeDataTest, EqualWithoutNeighbours) {
	std::vector<NodeData> v1;
	v1.emplace_back(1, 0, 1);
	v1.emplace_back(2, 0, 2);
	v1.emplace_back(1, 1, 3);
	v1.emplace_back(1, 0, 4);

	EXPECT_FALSE(v1[0] == v1[1]);
	EXPECT_FALSE(v1[0] == v1[2]);
	EXPECT_TRUE(v1[0] == v1[3]);
	EXPECT_FALSE(v1[1] == v1[2]);
	EXPECT_FALSE(v1[1] == v1[3]);
	EXPECT_FALSE(v1[2] == v1[3]);
}
