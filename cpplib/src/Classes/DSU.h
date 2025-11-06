#pragma once
#include <vector>
#include "../BaseHeaders/BaseTypes.h"
namespace cpplib {
    // Disjoint Set Union (DSU) data structure
    // Used for finding connected components in a graph
    class DSU {
    public:
        using AtomIndex = basic_types::AtomIndex;
    private:
        std::vector<AtomIndex> parent;
        std::vector<AtomIndex> set_id; // O(1) time to find set_id
        std::vector<std::vector<AtomIndex>> set_elements;
        AtomIndex count_components;

    public:
        explicit DSU(AtomIndex n) : parent(n), set_id(n), set_elements(n), count_components(n) {
            for (AtomIndex i = 0; i < n; ++i) {
                parent[i] = i;
                set_id[i] = i;
                set_elements[i] = { i };
            }
        }

        AtomIndex find(AtomIndex x) {
            if (parent[x] != x) {
                parent[x] = find(parent[x]);
            }
            return parent[x];
        }

        AtomIndex unite(AtomIndex x, AtomIndex y) {
            AtomIndex rootX = find(x);
            AtomIndex rootY = find(y);

            if (rootX == rootY) return set_id[rootX];

            AtomIndex new_set;
            AtomIndex old_set;
            if (set_id[rootX] < set_id[rootY]) {
                new_set = set_id[rootX];
                old_set = set_id[rootY];
                parent[rootY] = rootX;
            }
            else {
                new_set = set_id[rootY];
                old_set = set_id[rootX];
                parent[rootX] = rootY;
            }

            for (AtomIndex elem : set_elements[old_set]) {
                set_id[elem] = new_set;
            }

            set_elements[new_set].insert(
                set_elements[new_set].end(),
                set_elements[old_set].begin(),
                set_elements[old_set].end()
            );
            set_elements[old_set].clear();

            count_components--;
            return new_set;
        }

        inline AtomIndex get_set_id(AtomIndex element) const {
            return set_id[element];
        }

        inline const std::vector<AtomIndex>& get_elements(AtomIndex setid) const {
            return set_elements[setid];
        }
        inline AtomIndex get_count_components() const noexcept {
            return count_components;
        }
        const std::vector<std::vector<AtomIndex>>& get_components_ref() const noexcept {
            return set_elements;
        }
    };
} // namespace cpplib