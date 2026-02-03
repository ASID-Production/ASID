#pragma once

#include <concepts>
#include <utility>
#include "BaseTypes.h"

namespace cpplib {
    // Concepts
    /**
     * Access the first atom from a bond-like object.
     * @tparam T Type of the bond.
     * @param bond Bond-like object that provides a `first` member.
     * @returns The value of `bond.first`.
     */
    /**
     * Access the second atom from a bond-like object.
     * @tparam T Type of the bond.
     * @param bond Bond-like object that provides a `second` member.
     * @returns The value of `bond.second`.
     */
    namespace concept_support {
        template<typename T>
        constexpr /**
         * Obtain the first atom from a bond-like object.
         *
         * @tparam T Type of the bond-like object providing a `first` member.
         * @param bond Bond-like object whose `first` member will be returned.
         * @returns The `first` member of `bond`.
         */
        auto getAtom1(const T& bond) {
            if constexpr (requires { bond.first; }) return bond.first;
        }

        template<typename T>
        constexpr /**
         * Retrieve the second atom from a bond-like object.
         *
         * @tparam T Type of the bond-like object.
         * @param bond Bond-like object that exposes a member named `second`.
         * @returns The `second` member of `bond`.
         */
        auto getAtom2(const T& bond) {
            if constexpr (requires { bond.second; }) return bond.second;
        }
    }
    template<typename T>
    concept BondConcept = requires(const T & bond) {
        requires std::integral<decltype(concept_support::getAtom1(bond))>;
        requires std::same_as<decltype(concept_support::getAtom1(bond)), decltype(concept_support::getAtom2(bond))>;
        requires std::constructible_from<T,
            decltype(concept_support::getAtom1(bond)),
            decltype(concept_support::getAtom2(bond))>;
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