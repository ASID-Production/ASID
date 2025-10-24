#pragma once

#include <cstdint>
#include <bitset>
#include <array>
#include <type_traits>

namespace cpplib::constants {
    static constexpr int maxNeighbours = 100;
    constexpr char mend_size = 119;
    constexpr std::array<const char*, mend_size> mend{ "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Me", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };
} // namespace cpplib::constants

namespace cpplib::basic_types {
    // Fundamental Types
    using AtomIndex = int_fast32_t;
    using MoleculeIndex = int_fast32_t;
    using HType = int_fast8_t;
    using AtomTypeBase = int_fast8_t;
    using FloatingPointType = float;
    using TypeBitset = ::std::bitset<constants::mend_size>;

    // Maximum size type
    using size_type = std::conditional_t<(sizeof(AtomIndex) >= sizeof(MoleculeIndex)),
        AtomIndex, MoleculeIndex>;

    // Assets
    static_assert(std::is_integral_v<AtomTypeBase>, "AtomTypeBase must be integral");
    static_assert(std::is_integral_v<HType>, "HType must be integral");
    static_assert(std::is_integral_v<AtomIndex>, "AtomIndex must be integral");
    static_assert(sizeof(AtomIndex) >= sizeof(int32_t), "AtomIndex must be at least 32-bit");
    static_assert(std::is_integral_v<MoleculeIndex>, "MoleculeIndex must be integral");
    static_assert(sizeof(MoleculeIndex) >= sizeof(int32_t), "MoleculeIndex must be at least 32-bit");
    static_assert(TypeBitset().size() == constants::mend_size, "Size of TypeBitset must be constants::mend_size");
    static_assert(std::is_floating_point_v<FloatingPointType>, "FloatingPointType must be floating point");

} // namespace cpplib::base_types
