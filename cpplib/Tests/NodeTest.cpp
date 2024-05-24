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
#include "../src/Functions/AllInOneAndCurrent.h"
#include <vector>

using namespace cpplib::currents;
using NodeData = cpplib::Node<AtomTypeData>;
using NodeRequest = cpplib::Node<AtomTypeRequest>;

TEST(NodeRequestTest, Constructors) {

	EXPECT_NO_THROW({NodeRequest a1(1, 1, 1); });
	NodeRequest a1(1, 1, 1);
	NodeRequest::NeighboursType v1;
	v1.push_back(&a1);
	EXPECT_NO_THROW({NodeRequest a2(1, 1, 1); });
	EXPECT_NO_THROW({NodeRequest a3; });
	EXPECT_NO_THROW({NodeRequest a4(a1); });
	EXPECT_NO_THROW({NodeRequest a5(std::move(a1)); });
	EXPECT_NO_THROW({NodeRequest a6(1,1,1,std::move(v1)); });
}
TEST(NodeDataTest, Constructors) {

	EXPECT_NO_THROW({NodeData a1(1, 1, 1); });
	NodeData a1(1, 1, 1);
	NodeData::NeighboursType v1;
	v1.push_back(&a1);
	EXPECT_NO_THROW({NodeData a2(1, 1, 1); });
	EXPECT_NO_THROW({NodeData a3; });
	EXPECT_NO_THROW({NodeData a4(a1); });
	EXPECT_NO_THROW({NodeData a5(std::move(a1)); });
	EXPECT_NO_THROW({NodeData a6(1,1,1,std::move(v1)); });
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
TEST(NodeDataTest, NeighboursSorting) {
	std::vector<NodeData> node;
	node.emplace_back(0, 0, 0);
	node.emplace_back(17, 0, 1);
	node.emplace_back(6, 0, 2);
	node.emplace_back(6, 1, 3);
	node.emplace_back(6, 0, 4);
	node.emplace_back(6, 1, 5);
	node.emplace_back(6, 0, 6);
	node.emplace_back(6, 1, 7);
	node.emplace_back(6, 1, 8);
	node.emplace_back(17, 0, 9);
	node.emplace_back(7, 0, 10);
	node.emplace_back(7, 0, 11);
	node.emplace_back(6, 0, 12);
	node.emplace_back(6, 0, 13);
	node.emplace_back(8, 0, 14);
	node.emplace_back(6, 3, 15);
	node.emplace_back(7, 0, 17);
	node.emplace_back(7, 0, 16);
	node.emplace_back(7, 0, 19);
	node.emplace_back(7, 1, 20);
	node.emplace_back(6, 1, 18);
	node.emplace_back(6, 0, 21);
	node.emplace_back(8, 0, 22);
	node.emplace_back(6, 3, 23);

	node[1].addNeighboursVector({&(node[2])});
	node[2].addNeighboursVector({&(node[1]), &(node[3]), &(node[4])});
	node[3].addNeighboursVector({&(node[5]), &(node[2])});
	node[4].addNeighboursVector({&(node[7]), &(node[6]), &(node[2])});
	node[5].addNeighboursVector({&(node[3]), &(node[8])});
	node[6].addNeighboursVector({&(node[9]), &(node[8]), &(node[4])});
	node[7].addNeighboursVector({&(node[10]), &(node[4])});
	node[8].addNeighboursVector({&(node[5]), &(node[6])});
	node[9].addNeighboursVector({&(node[6])});
	node[10].addNeighboursVector({&(node[11]), &(node[7])});
	node[11].addNeighboursVector({&(node[10]), &(node[12]), &(node[13])});
	node[12].addNeighboursVector({&(node[14]), &(node[11]), &(node[15])});
	node[13].addNeighboursVector({&(node[11]), &(node[17]), &(node[16])});
	node[14].addNeighboursVector({&(node[12])});
	node[15].addNeighboursVector({&(node[12])});
	node[16].addNeighboursVector({&(node[18]), &(node[13])});
	node[17].addNeighboursVector({&(node[19]), &(node[20]), &(node[13])});
	node[18].addNeighboursVector({&(node[16]), &(node[20])});
	node[19].addNeighboursVector({&(node[17]), &(node[21])});
	node[20].addNeighboursVector({&(node[17]), &(node[18])});
	node[21].addNeighboursVector({&(node[22]), &(node[19]), &(node[23])});
	node[22].addNeighboursVector({&(node[21])});
	node[23].addNeighboursVector({&(node[21])});

	std::vector<NodeData> node2;
	node2.emplace_back(0, 0, 0);
	node2.emplace_back(17, 0, 1);
	node2.emplace_back(6, 0, 2);
	node2.emplace_back(6, 1, 3);
	node2.emplace_back(6, 0, 4);
	node2.emplace_back(6, 1, 5);
	node2.emplace_back(6, 0, 6);
	node2.emplace_back(6, 1, 7);
	node2.emplace_back(6, 1, 8);
	node2.emplace_back(17, 0, 9);
	node2.emplace_back(7, 0, 10);
	node2.emplace_back(7, 0, 11);
	node2.emplace_back(6, 0, 12);
	node2.emplace_back(6, 0, 13);
	node2.emplace_back(8, 0, 14);
	node2.emplace_back(6, 3, 15);
	node2.emplace_back(7, 0, 17);
	node2.emplace_back(7, 0, 16);
	node2.emplace_back(7, 0, 19);
	node2.emplace_back(7, 1, 20);
	node2.emplace_back(6, 1, 18);
	node2.emplace_back(6, 0, 21);
	node2.emplace_back(8, 0, 22);
	node2.emplace_back(6, 3, 23);

	node2[1].addBondWithSort(node2[2]);
	node2[2].addBondWithSort(node2[3]);
	node2[2].addBondWithSort(node2[4]);
	node2[3].addBondWithSort(node2[5]);
	node2[4].addBondWithSort(node2[6]);
	node2[4].addBondWithSort(node2[7]);
	node2[5].addBondWithSort(node2[8]);
	node2[6].addBondWithSort(node2[8]);
	node2[6].addBondWithSort(node2[9]);
	node2[7].addBondWithSort(node2[10]);
	node2[10].addBondWithSort(node2[11]);
	node2[11].addBondWithSort(node2[12]);
	node2[11].addBondWithSort(node2[13]);
	node2[12].addBondWithSort(node2[14]);
	node2[12].addBondWithSort(node2[15]);
	node2[13].addBondWithSort(node2[16]);
	node2[13].addBondWithSort(node2[17]);
	node2[16].addBondWithSort(node2[18]);
	node2[17].addBondWithSort(node2[19]);
	node2[17].addBondWithSort(node2[20]);
	node2[18].addBondWithSort(node2[20]);
	node2[19].addBondWithSort(node2[21]);
	node2[21].addBondWithSort(node2[22]);
	node2[21].addBondWithSort(node2[23]);

	for (size_t i = 0; i < 23; i++) {
		ASSERT_EQ(node[i].neighboursSize(), node2[i].neighboursSize()) << "if i = " << i;
		for (AtomIndex j = 0; j < node[i].neighboursSize(); j++) {
			const auto& nei1 = (*(node[i].getNeighbour(j)));
			const auto& nei2 = (*(node2[i].getNeighbour(j)));
			if (!(nei1 == nei2))
				FAIL();
		}
	}
}

