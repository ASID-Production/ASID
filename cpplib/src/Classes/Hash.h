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
#include "MoleculeGraph.h"
#include "Engine.h"
#include <cassert>
#include <vector>
#include <cstdint>

template<class A, class H, class AI>
class Hash {
public:
	// Definitions
	using hash8 = uint8_t;
	using hash16 = uint16_t;
	using hash32 = uint32_t;
	using hash64 = uint64_t;

	// Current settings
	using hash_single = hash16;
	using hash_full = hash64;

	using NodeType = Node<A, H, AI>;
	using size_type = AI;

private:
	// Data
	std::vector<hash_full> hash;
	template <class MI>
	explicit constexpr Hash(const MoleculeGraph<A, H, AI, MI>& nodes) = delete;

public:
	// Constructors
	explicit constexpr Hash(const std::vector<NodeType>& nodes) {
		AI s = nodes.size();
		hash.resize(s, 0);
		const std::vector<hash_single> monohash = createMonohash(nodes);
		// Second-Fourth levels
		for (AI i = 0; i < s; i++) {
			const auto& t = i;
			hash[t] = hash_recursive(t, nodes[t], monohash, 3);
		}
		std::sort(hash.begin(), hash.end());
	}
	explicit constexpr Hash(const std::vector<AI>& ai, const std::vector<NodeType>& nodes) {
		AI s = ai.size();
		hash.resize(s, 0);
		const std::vector<hash_single> monohash = createMonohash(nodes);
		// Second-Fourth levels
		for (AI i = 0; i < s; i++) {
			hash.operator[](i) = hash_recursive(ai[i], nodes[ai[i]], monohash, 3);
		}
		std::sort(hash.begin(), hash.end());
	}
	bool operator==(const Hash& other) const {
		size_type s = hash.size();
		for (size_type i = 0; i < s; i++) {
			if (hash[i] != other.hash[i])
				return false;
		}
		return true;
	}
private:
	constexpr std::vector<hash_single> createMonohash(const std::vector<NodeType>& nodes) const {
		auto size = nodes.size();
		std::vector<hash_single> monohash(size);
		for (size_type i = 0; i < size; i++) {
			monohash[i] = static_cast<hash_single>(nodes[i].getType()) + static_cast<hash_single>(static_cast<hash_single>(nodes[i].getHAtoms()) << 8) + static_cast<hash_single>(nodes[i].neighboursSize() << 12);
		}
		return monohash;
	}
	constexpr std::vector<hash_single> createMonohash(const std::vector<AI>& ai, const std::vector<NodeType>& nodes) const {

		auto size = nodes.size(); 
		std::vector<hash_single> monohash(size);
		for (size_type i = 0; i < size; i++) {
			monohash[i] = static_cast<hash_single>(nodes[ai[i]].getType()) + static_cast<hash_single>(static_cast<hash_single>(nodes[ai[i]].getHAtoms()) << 8) + static_cast<hash_single>(nodes[ai[i]].neighboursSize() << 12);
		}
		return monohash;
	}

	hash_full hash_recursive(const AI cur, const NodeType& node, const std::vector<hash_single>& mono, const char depth) {
		if (depth == 0)
			return mono[cur];
		hash_full res = mono[cur];
		auto size = node.neighboursSize();

		for (decltype(size) i = 0; i < size; i++) {
			res += hash_recursive(node.getNeighbour(i)->getID(), *(node.getNeighbour(i)), mono, depth - 1) << 16;
		}
		return res;
	}
};
