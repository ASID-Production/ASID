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
#include "../src/BaseHeaders/Currents.h"
#include <vector>

using namespace cpplib;
using namespace cpplib::currents;

// SimpleAtom
TEST(SimpleAtomTest, BasicFunctionality) {
    SimpleAtom a1(1);
    SimpleAtom a2(2);

    EXPECT_EQ(a1.get_bitset().count(), 1);
    EXPECT_TRUE(a1.contains(1));
    EXPECT_FALSE(a1.contains(2));
    EXPECT_TRUE(a1.intersect(a1));
    EXPECT_FALSE(a1.intersect(a2));
    EXPECT_EQ(static_cast<SimpleAtom::AtomTypeBase>(a1), 1);
}
TEST(SimpleAtomTest, Comparison) {
    SimpleAtom a1(1), a2(2), a3(1);

    EXPECT_EQ(a1, a3);
    EXPECT_NE(a1, a2);
    EXPECT_LT(a1, a2);
    EXPECT_GT(a2, a1);
}

// CompositeAtom
TEST(CompositeAtomTest, BasicOperations) {
    CompositeAtom atom(1);
    atom.AddType(2);
    atom.AddType(3);

    EXPECT_TRUE(atom.contains(1));
    EXPECT_TRUE(atom.contains(2));
    EXPECT_TRUE(atom.contains(3));
    EXPECT_FALSE(atom.contains(4));
    EXPECT_EQ(atom.get_bitset().count(), 3);
}
TEST(CompositeAtomTest, Intersection) {
    CompositeAtom a1(1);
    a1.AddType(2);

    CompositeAtom a2(2);
    a2.AddType(3);

    CompositeAtom a3(3);
    a3.AddType(4);

    EXPECT_TRUE(a1.intersect(a2));
    EXPECT_FALSE(a1.intersect(a3));
}

// Coord
TEST(CoordTest, Intersection) {
    Coord c1(5);
    Coord c2(3, 7);
    Coord c3(8, 10);

    EXPECT_TRUE(c1.intersect(c2));
    EXPECT_FALSE(c1.intersect(c3));
    EXPECT_EQ(c1.getLow(), 5);
    EXPECT_EQ(c1.getHigh(), 5);
    EXPECT_EQ(c2.getLow(), 3);
    EXPECT_EQ(c2.getHigh(), 7);
}

// NeighboursType
TEST(NeighboursTypeTest, BasicOperations) {
    NeighboursType neighbors;

    neighbors.push_back(10);
    neighbors.push_back(20);
    neighbors.push_back(30);

    EXPECT_EQ(neighbors.size(), 3);
    EXPECT_EQ(neighbors[0], 10);
    EXPECT_EQ(neighbors[1], 20);
    EXPECT_EQ(neighbors[2], 30);

    neighbors.erase(1);
    EXPECT_EQ(neighbors.size(), 2);
    EXPECT_EQ(neighbors[0], 10);
    EXPECT_EQ(neighbors[1], 30);

    EXPECT_TRUE(neighbors.exchange(10, 15));
    EXPECT_EQ(neighbors[0], 15);
}
TEST(NeighboursTypeTest, SortingAndShifting) {
    NeighboursType neighbors;
    neighbors.push_back(30);
    neighbors.push_back(10);
    neighbors.push_back(20);

    neighbors.simpleSort();
    EXPECT_EQ(neighbors[0], 10);
    EXPECT_EQ(neighbors[1], 20);
    EXPECT_EQ(neighbors[2], 30);

    neighbors.addShift(5);
    EXPECT_EQ(neighbors[0], 15);
    EXPECT_EQ(neighbors[1], 25);
    EXPECT_EQ(neighbors[2], 35);
}
TEST(NeighboursTypeTest, SomeEdgeCases) {
    NeighboursType neighbors;

    // Up to max
    for (int i = 0; i < NeighboursType::maxNeighbours; i++) {
        neighbors.push_back(i);
    }

    // wrong exchange
    EXPECT_FALSE(neighbors.exchange(-1, 10));
}

// Node
class NodeTest : public ::testing::Test {
protected:
    void SetUp() override {
        nodes.resize(3);
        nodes[0] = Node<SimpleAtom>(SimpleAtom(1), 0, 0);
        nodes[1] = Node<SimpleAtom>(SimpleAtom(2), 0, 1);
        nodes[2] = Node<SimpleAtom>(SimpleAtom(3), 0, 2);
    }

    std::vector<Node<SimpleAtom>> nodes;
};
TEST_F(NodeTest, BondManagement) {
    nodes[0].addBondSimple(nodes[1]);
    nodes[0].addBondSimple(nodes[2]);

    EXPECT_EQ(nodes[0].neighboursSize(), 2);
    EXPECT_EQ(nodes[1].neighboursSize(), 1);
    EXPECT_EQ(nodes[2].neighboursSize(), 1);

    EXPECT_TRUE(nodes[0].isNeighbour(nodes[1]));
    EXPECT_FALSE(nodes[1].isNeighbour(nodes[2]));

    nodes[0].deleteBond(nodes[1]);
    EXPECT_EQ(nodes[0].neighboursSize(), 1);
    EXPECT_EQ(nodes[1].neighboursSize(), 0);
}
TEST_F(NodeTest, ComparisonOperators) {
    Node<SimpleAtom> n1(SimpleAtom(1), 2, 0);
    Node<SimpleAtom> n2(SimpleAtom(1), 1, 1);
    Node<SimpleAtom> n3(SimpleAtom(2), 1, 2);

    EXPECT_TRUE(n1.RawLess(n3));
    EXPECT_TRUE(n1.RawMore(n2));
    EXPECT_LT(n1, n3);
    EXPECT_GT(n1, n2);
}
TEST_F(NodeTest, SwapOperation) {
    nodes[0].addBondSimple(nodes[1]);
    auto original_neighbor = nodes[0].getNeighbour(0);

    std::swap(nodes[0], nodes[2]);

    EXPECT_EQ(nodes[2].getID(), 2);
    EXPECT_EQ(original_neighbor->getID(), 1);
    EXPECT_EQ(nodes[2].neighboursSize(), 1);
    EXPECT_EQ(nodes[0].neighboursSize(), 0);
}
TEST_F(NodeTest, AddBondWithSort) {
    nodes[0].addBondWithSort(nodes[1]);
    nodes[0].addBondWithSort(nodes[2]);

    // Sort Check
    const auto& neighbours = nodes[0].getNeighboursVector();
    EXPECT_TRUE(std::is_sorted(neighbours.begin(), neighbours.end()));
}
TEST_F(NodeTest, FindNeighbour) {
    nodes[0].addBondSimple(nodes[1]);
    auto* neighbour = nodes[0].getNeighbour(0);

    // Positive
    EXPECT_EQ(nodes[0].findNeighbour(neighbour), 0);

    // Negative
    EXPECT_EQ(nodes[0].findNeighbour(&nodes[2]), Node<SimpleAtom>::AtomIndex(-1));
}
TEST_F(NodeTest, HasNeighbours) {
    EXPECT_FALSE(nodes[0].hasNeighbours());
    nodes[0].addBondSimple(nodes[1]);
    EXPECT_TRUE(nodes[0].hasNeighbours());
}
TEST_F(NodeTest, CalculateCoordEdgeCases) {
    Node<SimpleAtom> node(SimpleAtom(1), 0, 0);

    // 0 neighbours
    node.calculateCoord();
    EXPECT_EQ(node.getCoord().getLow(), 0);

    // max value
    for (int i = 0; i < constants::maxNeighbours; i++) {
        node.addBondSimple(nodes[1]);
    }
    node.calculateCoord();
    EXPECT_EQ(node.getCoord().getHigh(), constants::maxNeighbours);
}

//  Node with CompositeAtom
TEST(CompositeNodeTest, BasicOperations) {
    CompositeAtom compAtom(1);
    compAtom.AddType(2);
    compAtom.AddType(3);

    Node<CompositeAtom> node(compAtom, 2, 100);

    EXPECT_TRUE(node.getType().contains(1));
    EXPECT_TRUE(node.getType().contains(2));
    EXPECT_FALSE(node.getType().contains(4));
    EXPECT_EQ(node.getHAtoms(), 2);
    EXPECT_EQ(node.getID(), 100);
}
TEST(CompositeNodeTest, NeighbourManagement) {
    std::vector<Node<CompositeAtom>> nodes(3);

    CompositeAtom types[3] = {
        CompositeAtom(1),
        CompositeAtom(2),
        CompositeAtom(3)
    };

    for (int i = 0; i < 3; ++i) {
        types[i].AddType(i + 10); // Добавляем дополнительный тип
        nodes[i] = Node<CompositeAtom>(types[i], i, i);
    }

    nodes[0].addBondSimple(nodes[1]);
    nodes[0].addBondSimple(nodes[2]);

    EXPECT_EQ(nodes[0].neighboursSize(), 2);
    EXPECT_EQ(nodes[1].neighboursSize(), 1);

    // Проверка связи между разными типами атомов
    EXPECT_TRUE(nodes[0].isNeighbour(nodes[1]));
    EXPECT_FALSE(nodes[1].isNeighbour(nodes[2]));
}
TEST(CompositeNodeTest, ComparisonOperations) {
    CompositeAtom type1(1);
    type1.AddType(10);

    CompositeAtom type2(2);
    type2.AddType(20);

    Node<CompositeAtom> node1(type1, 2, 1);
    Node<CompositeAtom> node2(type2, 3, 2);
    Node<CompositeAtom> node3(type1, 2, 3); // Такой же тип как node1

    // Проверка сравнения по базовому типу
    EXPECT_TRUE(node1.RawLess(node2));
    EXPECT_TRUE(node2.RawMore(node1));

    // Проверка операторов сравнения
    EXPECT_LT(node1, node2);
    EXPECT_GT(node2, node1);

    // Узлы с одинаковым базовым типом
    EXPECT_EQ(node1 == node3, true);
    EXPECT_NE(node1, node2);
}
TEST(CompositeNodeTest, SwapOperation) {
    std::vector<Node<CompositeAtom>> nodes(3);

    for (int i = 0; i < 3; ++i) {
        CompositeAtom atom(i + 1);
        atom.AddType(i + 10);
        nodes[i] = Node<CompositeAtom>(atom, i, i);
    }

    nodes[0].addBondSimple(nodes[1]);
    nodes[0].addBondSimple(nodes[2]);

    auto original_neighbor1 = nodes[0].getNeighbour(0);
    auto original_neighbor2 = nodes[0].getNeighbour(1);

    // Сохраняем оригинальные связи
    const auto neighbors_before = nodes[0].getNeighboursVector();

    // Выполняем swap
    std::swap(nodes[0], nodes[1]);

    // Проверяем обновление связей
    EXPECT_EQ(nodes[1].neighboursSize(), 2);
    EXPECT_EQ(nodes[0].neighboursSize(), 1);

    // Проверяем что соседи обновили свои ссылки
    for (size_t i = 0; i < nodes[1].neighboursSize(); ++i) {
        Node<CompositeAtom>* neighbor = nodes[1].getNeighbour(i);
        EXPECT_TRUE(neighbor->isNeighbour(nodes[1]));
    }
}
TEST(CompositeNodeTest, NotExactCompareFunctionality) {
    // Подготовка композитного атома
    CompositeAtom compType(1);
    compType.AddType(2);
    compType.AddType(3);

    Node<CompositeAtom> compNode(compType, 1, 100);
    compNode.calculateCoord(); // Координата (1+1=2, 2)

    // Подходящий простой атом (входит в композитный)
    Node<SimpleAtom> matchingSimpleNode(SimpleAtom(2), 1, 200);
    matchingSimpleNode.calculateCoord(); // Координата (1, 1)

    // Неподходящий простой атом (тип не входит)
    Node<SimpleAtom> nonMatchingTypeNode(SimpleAtom(4), 1, 300);

    // Неподходящий по количеству H-атомов
    Node<SimpleAtom> tooManyHAtomsNode(SimpleAtom(2), 2, 400);

    // Неподходящий по координатам
    Node<SimpleAtom> nonMatchingCoordNode(SimpleAtom(2), 1, 500);
    nonMatchingCoordNode.setCoord(Coord(10, 10));

    // Проверки
    EXPECT_TRUE(NotExactCompare(compNode, matchingSimpleNode));
    EXPECT_FALSE(NotExactCompare(compNode, nonMatchingTypeNode));
    EXPECT_FALSE(NotExactCompare(compNode, tooManyHAtomsNode));
    EXPECT_FALSE(NotExactCompare(compNode, nonMatchingCoordNode));
}
TEST(CompositeNodeTest, EdgeCases) {
    // Узел с пустым композитным атомом
    CompositeAtom emptyAtom;
    Node<CompositeAtom> emptyNode(emptyAtom, 0, 0);

    // Узел с одним типом
    CompositeAtom singleType(5);
    Node<CompositeAtom> singleTypeNode(singleType, 3, 1);

    // Проверка contains
    EXPECT_FALSE(emptyNode.getType().contains(1));
    EXPECT_TRUE(singleTypeNode.getType().contains(5));
    EXPECT_FALSE(singleTypeNode.getType().contains(1));

    // Проверка сравнения с пустым узлом
    Node<CompositeAtom> anotherEmptyNode(emptyAtom, 0, 2);
    EXPECT_EQ(emptyNode == anotherEmptyNode, true);

    // Проверка NotExactCompare с пустым узлом
    Node<SimpleAtom> simpleNode(SimpleAtom(1), 0, 3);
    EXPECT_FALSE(NotExactCompare(emptyNode, simpleNode));
}

// Bond
TEST(BondTest, BasicOperations) {
    Bond b1(1, 2);
    Bond b2(2, 1);
    Bond b3(3, 4);

    EXPECT_EQ(b1.first, 1);
    EXPECT_EQ(b1.second, 2);
    EXPECT_NE(b1, b2);
    b2.validate();
    EXPECT_EQ(b1, b2);
    EXPECT_NE(b1, b3);
    EXPECT_LT(b1, b3);

    EXPECT_EQ(b1.ToStr(), "(1, 2)");
}

// BondEx
TEST(BondExTest, ExtendedFunctionality) {
    BondEx b1(1, 2, 1.5f);
    BondEx b2(2, 1, 1.5f);
    BondEx b3(1, 2, 2.0f);

    EXPECT_EQ(b1, b2);
    EXPECT_EQ(b1.length, 1.5f);
    EXPECT_LT(b1, b3);
    EXPECT_EQ(b1.ToStr(), "(1, 2, {\"distance\": 1.500000})");
}

// NotExactCompare
TEST(NotExactCompareTest, Basic) {
    Node<CompositeAtom> compNode;
    CompositeAtom type(1);
    type.AddType(2);
    compNode.setType(type);
    compNode.setHAtoms(1);
    compNode.setCoord(Coord(1, 2));

    Node<SimpleAtom> simpleNode(SimpleAtom(2), 2, 0);

    EXPECT_TRUE(NotExactCompare(compNode, simpleNode));

    Node<SimpleAtom> simpleNode2(SimpleAtom(3), 1, 0);
    EXPECT_FALSE(NotExactCompare(compNode, simpleNode2));

    Node<SimpleAtom> simpleNode3(SimpleAtom(1), 3, 0);
    EXPECT_FALSE(NotExactCompare(compNode, simpleNode3));
}

// ExactCompare
TEST(ExactCompareTest, Basic) {
    CompositeAtom compType(1);
    compType.AddType(2);
    Node<CompositeAtom> compNode(compType, 2, 100);
    compNode.calculateCoord();

    // exactMatch
    Node<SimpleAtom> exactMatch(SimpleAtom(2), 2, 200);
    exactMatch.calculateCoord();
    EXPECT_TRUE(ExactCompare(compNode, exactMatch));

    // hMismatch
    Node<SimpleAtom> hMismatch(SimpleAtom(2), 1, 300);
    EXPECT_FALSE(ExactCompare(compNode, hMismatch));
}