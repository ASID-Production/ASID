#pragma once
#include <cstdint>
#include <type_traits>

namespace cpplib {
	class XAtom;
	struct Bond;
	struct BondEx;
	class Coord;
	template<class A> class Node;
	template<class A> class MoleculeGraph;
	class SearchGraph;
	class Distances;
	class FindMolecules;
	class FindGeometry;
	class FAM_Struct;
	class FAM_Cell;
	class SearchDataInterface;
	class ParseData;

	namespace geometry {
		template<class T> class Point;
		template<class T> class Matrix;
		template<class T> struct Cell;
		template<class T> struct Symm;
	}

	namespace currents {
		using AtomIndex = int_fast32_t;
		using MoleculeIndex = int_fast32_t;
		using DistancesIndexType = int_fast8_t;
		using HType = int_fast8_t;
		using AtomTypeRequest = XAtom;
		using AtomTypeData = int_fast8_t;
		using FloatingPointType = float;

		using size_type = ::std::conditional<(sizeof(AtomIndex) > sizeof(MoleculeIndex)), AtomIndex, MoleculeIndex>::type;
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
}