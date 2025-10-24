#pragma once

#include <concepts>
#include <utility>
#include "BaseTypes.h"

namespace cpplib {
    // Concepts
    // Universal functions
    namespace concept_support {
        template<typename T>
        constexpr auto getAtom1(const T& bond) {
            if constexpr (requires { bond.atom1; }) return bond.atom1;
            else if constexpr (requires { bond.first; }) return bond.first;
            else if constexpr (requires { bond.first(); }) return bond.first();
            else if constexpr (requires { bond.getFirstAtom(); }) return bond.getFirstAtom();
        }

        template<typename T>
        constexpr auto getAtom2(const T& bond) {
            if constexpr (requires { bond.atom2; }) return bond.atom2;
            else if constexpr (requires { bond.second; }) return bond.second;
            else if constexpr (requires { bond.second(); }) return bond.second();
            else if constexpr (requires { bond.getSecondAtom(); }) return bond.getSecondAtom();
        }
    }
    template<typename T>
    concept BondConcept = requires(const T & bond) {
        requires std::integral<decltype(concept_support::getAtom1(bond))>;
        requires std::integral<decltype(concept_support::getAtom2(bond))>;
    };

    template<typename T>
    concept AtomTypeConcept = requires(const T & a, const T & b, basic_types::AtomTypeBase base) {
        // Constructors
            requires std::is_default_constructible_v<T>&& std::is_copy_constructible_v<T>;
    T(base); // AtomTypeBase constructor

    // Operators
    { a == b } -> std::same_as<bool>;
    { a <=> b } -> std::convertible_to<std::strong_ordering>;

    // Methods
    { a.get_bitset() } noexcept -> std::same_as<basic_types::TypeBitset>;
    { a.contains(base) } -> std::same_as<bool>;
    { a.intersect(b) } -> std::same_as<bool>;
    { static_cast<typename basic_types::AtomTypeBase>(a) } -> std::same_as<basic_types::AtomTypeBase>;
    };
}