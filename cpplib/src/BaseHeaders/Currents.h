#pragma once
#include <cstdint>
#include <type_traits>
#include <bitset> // for TypeBitset
#include <vector>
#include <concepts>
#include "Support.h"

namespace cpplib {

	// Definitions (currents)
	namespace currents {
		using AtomIndex = int_fast32_t;
		using MoleculeIndex = int_fast32_t;
		using HType = int_fast8_t;
		using AtomTypeBase = int_fast8_t;
		using FloatingPointType = float;
		using TypeBitset = ::std::bitset<mend_size>;
	}

	// Concepts
	template<typename T>
	concept AtomTypeConcept = requires(const T & a, const T & b, currents::AtomTypeBase base) {
		// Constructors
			requires std::is_default_constructible_v<T>&& std::is_copy_constructible_v<T>;
	T(base); // AtomTypeBase constructor

	// Operators
	{ a == b } -> std::same_as<bool>;
	{ a <=> b } -> std::convertible_to<std::strong_ordering>;

	// Methods
	{ a.get_bitset() } noexcept -> std::same_as<currents::TypeBitset>;
	{ a.contains(base) } -> std::same_as<bool>;
	{ a.intersect(b) } -> std::same_as<bool>;
	{ static_cast<currents::AtomTypeBase>(a) } -> std::same_as<currents::AtomTypeBase>;
	};


	// Forward declarations
	class SimpleAtom;
	class CompositeAtom;
	struct Bond;
	struct BondEx;
	class Coord;
	template<AtomTypeConcept A> class Node;
	template<AtomTypeConcept A> class MoleculeCore;
	class SearchGraph;
	class Distances;
	class FindMolecules;
	class FindGeometry;
	class FAM_Struct;
	class FAM_Cell;
	class SearchDataInterface;
	class ParseData;
	namespace geometry {
		template<class T> struct Point;
		template<class T> class Matrix;
		template<class T> struct Cell;
		template<class T> struct Symm;
	}


	// Type alliases (currents)
	namespace currents {
		using AtomTypeRequest = CompositeAtom;
		using AtomTypeData = SimpleAtom;
		using size_type = ::std::conditional_t<(sizeof(AtomIndex) > sizeof(MoleculeIndex)), AtomIndex, MoleculeIndex>;
		using PointType = geometry::Point<FloatingPointType>;
		using BondType = Bond;
		using BondExType = BondEx;
		using DistancesType = Distances;
		using ParseIndexType = ::std::vector<size_type>;
		using CellType = geometry::Cell<FloatingPointType>;
		using SearchGraphType = SearchGraph;
		using SearchDataInterfaceType = SearchDataInterface;
		using FindMoleculesType = FindMolecules;
		using FindGeometryType = FindGeometry;
		using FAMStructType = FAM_Struct;
		using FAMCellType = FAM_Cell;
		using ParseDataType = ParseData;
		using SymmType = geometry::Symm<FloatingPointType>;
	}

	// Asserts
	namespace currents {
		static_assert(::std::is_integral_v<AtomTypeBase>, "AtomTypeBase must be integral");
		static_assert(::std::is_integral_v<HType>, "HType must be integral");
		static_assert(::std::is_integral_v<AtomIndex> && (sizeof(AtomIndex) >= sizeof(int32_t)), "AtomIndex must be integral");
		static_assert(::std::is_integral_v<MoleculeIndex> && (sizeof(MoleculeIndex) >= sizeof(int32_t)), "MoleculeIndex must be integral");
		static_assert(currents::TypeBitset().size() == mend_size, "Size of TypeBitset must be mend_size");
		static_assert(::std::is_floating_point_v<FloatingPointType>, "FloatingPointType must be floating point");
	}
}