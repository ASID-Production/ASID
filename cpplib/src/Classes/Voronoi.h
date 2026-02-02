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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/Bond.h"
#include "../Classes/Geometry.h"

namespace cpplib::voronoi {
	class Vertex;
	class Edge;
	class Face;
	class Cell;

	template<class T>
	using Container = ::std::unordered_set<T>;

	enum class State : char {
		DELETE = 0,
		VALID = 1,
		INVALID = 2,
		MODIFICATION = 3
	};

	class Object {
		size_t id;
		State state;
	public:
		inline size_t get_id() const noexcept {
			return id;
		}
		inline void set_id(size_t i) noexcept {
			id = i;
		}
		inline State get_state() const noexcept {
			return state;
		}
		inline void set_state(State s) noexcept {
			state = s;
		}
		constexpr explicit Object(size_t ID, State s = State::INVALID) noexcept
			: id(ID), state(s) {}
	};

	struct Vertex : public Object{
		// Types
		using PointType = geometry::Point<basic_types::FloatingPointType>;


		// Constants
		static constexpr typename PointType::value_type COMPARISON_EPSILON = 1.0 / (1 << 16);

		// Data
		PointType point;
		basic_types::FloatingPointType distance = basic_types::FloatingPointType(0.0);

		// Data (not owning)
		Container<Edge*> edges;
		Container<Face*> faces;

		// Constructors
		explicit Vertex(size_t ID) noexcept: Object(ID, State::VALID) {}
		Vertex(size_t ID, const PointType& p) noexcept : Object(ID, State::VALID), point(p) {}
		Vertex(size_t ID, PointType&& p) noexcept : Object(ID, State::VALID), point(std::move(p)) {}

		inline void add_edge(Edge* edge) {
			edges.insert(edge);
		}

		// Getters
		inline const PointType& get_point() const {
			return point;
		}
		inline const Container<Edge*>& get_edges() const {
			return edges;
		}

		// Special comparison
		inline friend bool operator==(const Vertex& a, const Vertex& b) {
			return std::abs(a.point[0] - b.point[0]) < COMPARISON_EPSILON &&
				std::abs(a.point[1] - b.point[1]) < COMPARISON_EPSILON &&
				std::abs(a.point[2] - b.point[2]) < COMPARISON_EPSILON;
		}


	};

	struct Edge : public Object {
		// Data (not owning)
		Container<Vertex*> vertices; // max 2
		Container<Face*> faces; // max 2

		// Types
		using PointType = geometry::Point<basic_types::FloatingPointType>;
		using PlaneType = geometry::Plane<basic_types::FloatingPointType>;

		Edge(size_t ID, Vertex* v1, Vertex* v2) : Object(ID), vertices({v1, v2}) {
			assert(v1 != NULL && v1 != nullptr);
			assert(v2 != NULL && v2 != nullptr);
			v1->edges.emplace(this);
			v2->edges.emplace(this);
		}


		inline State calculateState() noexcept {
			using enum State;
			if (vertices.size() == 2 && faces.size() == 2) {
				size_t counter = 0;
				for (const auto v : vertices)
				{
					auto curstate = v->get_state();
					if (curstate == VALID || curstate == MODIFICATION)
						counter++;
				}
				switch (counter) {
					case 2:
						set_state(VALID);
						break;
					case 1:
						set_state(MODIFICATION);
						break;
					case 0:
						set_state(DELETE);
						break;
					default:
						// Impossible
						break;
				}
			} else {
				set_state(INVALID);
			}
			return get_state();
		}

		PointType intersectSegmentPlane(PlaneType plane) {
			assert(calculateState() == State::MODIFICATION);
			auto it = vertices.begin();
			auto v1 = *it;
			it++;
			auto v2 = *it;

			PointType direction = v2->get_point() - v1->get_point();
			auto unnormalized_normal = PointType(plane.a[0], plane.a[1], plane.a[2]);
			auto denom = PointType::Scalar(unnormalized_normal, direction);

			// Check that Edge is not parallel to plane
			assert(std::abs(denom) >= 1e-6);


			auto t = -(plane.a[0] * v1->get_point()[0] +
					   plane.a[1] * v1->get_point()[1] +
					   plane.a[2] * v1->get_point()[2] +
					   plane.a[3]) / denom;

			assert(t > 0.0 && t < 1.0);

			return v1->get_point() + direction * t;
		}

	};

	struct Face : public Object {
		// Data (not owning)
		size_t owner_id;
		size_t other_id;
		char other_shiftcode;
		Container<Vertex*> vertices;
		Container<Edge*> edges;
		explicit Face(size_t ID) : Object(ID) {}
		inline State calculateState() {
			using enum State;

			if (vertices.size() == edges.size()) {
				size_t counter = 0;
				set_state(VALID);
				for (const auto e : edges)
				{
					auto estate = e->get_state();
					if (estate == VALID) {
						counter++;
					} else if (estate == MODIFICATION) {
						counter++;
						set_state(MODIFICATION);
					}
				}
				if (counter == 0) {
					set_state(DELETE);
				} else if (counter == 1) { // ERROR STATE!!!
					set_state(INVALID);
				}
			} else {
				set_state(INVALID);
			}
			return get_state();
		}

	};
	struct Cell {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;
		using ShiftType = geometry::Point<int8_t>;
		using BondWithShift = BondWithPoint<ShiftType>;
		using PlaneType = geometry::Plane<FloatingPointType>;

		static constexpr FloatingPointType limit = 2 * ::std::numeric_limits<FloatingPointType>::epsilon();

		// Data (Owning)
		std::vector<std::unique_ptr<Vertex>> vertices;
		std::vector<std::unique_ptr<Edge>> edges;
		std::vector<std::unique_ptr<Face>> faces;
		PointType center;
		int id = -1;


		static constexpr std::array<PointType, 8> base_vertices = {{
			PointType{-0.5, -0.5, -0.5}, // 0
			PointType{ 0.5, -0.5, -0.5}, // 1
			PointType{ 0.5,  0.5, -0.5}, // 2
			PointType{-0.5,  0.5, -0.5}, // 3
			PointType{-0.5, -0.5,  0.5}, // 4
			PointType{ 0.5, -0.5,  0.5}, // 5
			PointType{ 0.5,  0.5,  0.5}, // 6
			PointType{-0.5,  0.5,  0.5}  // 7
		}};

		// Indexes for each face of the cube (conter-clockwise from outside)
		static constexpr std::array<std::array<int, 4>, 6> face_indices = {{
			{4, 7, 6, 5}, // front face
			{0, 1, 2, 3}, // back face
			{0, 3, 7, 4}, // left face
			{1, 5, 6, 2}, // right face
			{0, 4, 5, 1}, // bottom face
			{3, 2, 6, 7}  // top face
		}};
		static constexpr std::array<char, 6> face_shiftcodes = {{
			22, // front face
			 4, // back face
			12, // left face
			14, // right face
			10, // bottom face
			16  // top face
		}};

		// Indexes for each edge of the cube
		static constexpr std::array<std::array<int, 2>, 12> edge_indices = {{
			{4, 7}, // edge 0: front-left
			{7, 6}, // edge 1: front-top
			{6, 5}, // edge 2: front-right
			{5, 4}, // edge 3: front-bottom
			{0, 1}, // edge 4: back-bottom
			{1, 2}, // edge 5: back-right
			{2, 3}, // edge 6: back-top
			{3, 0}, // edge 7: back-left
			{0, 4}, // edge 8: left-bottom
			{3, 7}, // edge 9: left-top
			{1, 5}, // edge 10: right-bottom
			{2, 6}  // edge 11: right-top
		}};
		static constexpr std::array<std::array<int, 4>, 6> face_edge_indices = {{
			{ 0,  1,  2,  3},  // front face
			{ 4,  5,  6,  7},  // back face
			{ 7,  9,  0,  8},  // left face
			{10,  2, 11,  5},  // right face
			{ 8,  3, 10,  4},  // bottom face
			{ 9,  1, 11,  6}   // top face
		}};
		

	public:
		Cell() = default;
		Cell(const PointType& c, int i) : center(c), id(i){
			vertices.reserve(128);
			edges.reserve(128);
            faces.reserve(64);
			for (size_t vert_id = 0; vert_id < 8; vert_id++) {
				vertices.emplace_back(std::make_unique<Vertex>(vert_id, base_vertices[vert_id] + center));
			}
			for (size_t edge_id = 0; edge_id < 12; edge_id++) {
				auto vertex1_ptr = vertices[edge_indices[edge_id][0]].get();
				auto vertex2_ptr = vertices[edge_indices[edge_id][1]].get();
				auto owner_ptr = std::make_unique<Edge>(edge_id, vertex1_ptr, vertex2_ptr); // smart pointer
				auto edge_ptr = owner_ptr.get(); // raw pointer
				edges.emplace_back(std::move(owner_ptr)); // smart pointer becomes invalid
				edge_ptr->set_state(State::VALID);
			}
			for (size_t face_id = 0; face_id < 6; face_id++) {
				auto face = std::make_unique<Face>(face_id);
				face->owner_id = id;
				face->other_id = id;
				face->other_shiftcode = face_shiftcodes[face_id];

				// Fill vertices
				for (int j = 0; j < 4; j++) {
					auto vertex_ptr = vertices[face_indices[face_id][j]].get();
					face->vertices.emplace(vertex_ptr);
					vertex_ptr->faces.emplace(face.get());
				}

				// Fill edges
				for (int j = 0; j < 4; j++) {
					auto edge_ptr = edges[face_edge_indices[face_id][j]].get();
					face->edges.emplace(edge_ptr);
					edge_ptr->faces.emplace(face.get());
				}

				face->set_state(State::VALID);
				faces.push_back(std::move(face));
			}
		}
		Vertex* add_vertex(const PointType& p) {
			auto temp = std::make_unique<Vertex>(vertices.size(), p);
			for (auto& v : vertices) {
				if (*v == *temp)
					return v.get();
			}
			auto simple_ptr = temp.get();
			vertices.emplace_back(std::move(temp));
			return simple_ptr;
		}
		void clipByPlaneAndAddNewFace(const PlaneType& clipping_plane, size_t id_of_another_cell, char another_shiftcode) {
			using enum State;

			// Check if the side is correct
			assert(clipping_plane.side(center) > 0);

			// Check cutting
			bool modified = false;
			// 1. Separate vertices to sides of plane
			for (auto& v : vertices)
			{
				if (v->get_state() == DELETE) {
					continue;
				}
				// calculate side:
				auto side = clipping_plane.side(v->get_point());
				if (side < -limit) {
					// Point cutted off
					v->set_state(DELETE);
					modified = true;
				} else if (side < limit) {
					// Point on the Face
					v->set_state(MODIFICATION);
				}
			}
			if (modified == false) {
				for (auto& v : vertices)
				{
					if (v->get_state() == MODIFICATION) {
						v->set_state(VALID);
					}
				}
				return;
			}

			// 2. Calculate States of edges and faces
			for (auto& e : edges)
			{
				if (e->get_state() == DELETE) {
					continue;
				}
				e->calculateState();
			}
			for (auto& f : faces)
			{
				if (f->get_state() == DELETE) {
					continue;
				}
				f->calculateState();
				assert(f->get_state() != INVALID);
			}

			// 3. Cut edges
			// So we need check all Edges, which ask about MODIFICATION
			for (auto& e : edges)
			{
				if (e->get_state() != MODIFICATION) {
					continue;
				}
				auto intersection = e->intersectSegmentPlane(clipping_plane);
				auto new_vertex_ptr = add_vertex(intersection);

				// delete vertex from set
				std::erase_if(e->vertices,
							  [](auto* ptr) {
								  if (ptr->get_state() == DELETE) {
									  return true;
								  }
								  return false;
							  });
				// and add new vertex to set
				e->vertices.insert(new_vertex_ptr);


				new_vertex_ptr->edges.insert(e.get()); // insert this edge to set of new vertex
				new_vertex_ptr->faces.insert(e->faces.begin(), e->faces.end()); // copy set of faces from the edge
				for (auto& f : e->faces) {
					f->vertices.emplace(new_vertex_ptr);
				}
				new_vertex_ptr->set_state(MODIFICATION); // Set, that Vertex is on the NEW FACE

				// Now, new vertex complete, so edge is VALID too:
				e->set_state(VALID);
			}

			// 4. Cut Faces, which need modification
			for (auto& f : faces) {
				if (f->get_state() != MODIFICATION) continue;

				// 4.1. Find and delete all unnesessary edges and vertices
				std::erase_if(f->edges, [](const auto* ptr) {
					return ptr->get_state() == DELETE;
					});
				std::erase_if(f->vertices, [](const auto* ptr) {
					return ptr->get_state() == DELETE;
					});

				// 4.2 Modify the face: add Edge

				auto iter_vertex = f->vertices.cbegin();
				Vertex* v1 = nullptr;
				Vertex* v2 = nullptr;

				// Find new vertices
				while (iter_vertex != f->vertices.cend()) {
					v1 = *iter_vertex;
					if (v1->get_state() == MODIFICATION) {
						iter_vertex++;
						break;
					}
					iter_vertex++;
				}
				while (iter_vertex != f->vertices.cend()) {
					v2 = *iter_vertex;
					if (v2->get_state() == MODIFICATION) {
						break;
					}
					iter_vertex++;
				}
				assert(iter_vertex != f->vertices.cend());

				const auto& new_edge = edges.emplace_back(std::make_unique<Edge>(edges.size(), v1, v2));
				new_edge->faces.emplace(f.get());
				f->edges.emplace(new_edge.get());
				new_edge->set_state(MODIFICATION);
				f->set_state(VALID);

			}

			// 5. Create new Face
			const auto& new_face = faces.emplace_back(std::make_unique<Face>(faces.size())); // smart pointer ref
			auto raw_face_ptr = new_face.get(); // raw pointer
			for (const auto& e : edges)
			{
				if (e->get_state() == MODIFICATION) {
					raw_face_ptr->edges.emplace(e.get());
					e->faces.emplace(raw_face_ptr);
					e->set_state(VALID);
					assert(e->calculateState() == VALID);
				}
			}
			for (const auto& v : vertices)
			{
				if (v->get_state() == MODIFICATION) {
					raw_face_ptr->vertices.emplace(v.get());
					v->faces.emplace(raw_face_ptr);
					v->set_state(VALID);
				}
			}
			raw_face_ptr->owner_id = id;
			raw_face_ptr->other_id = id_of_another_cell;
			raw_face_ptr->other_shiftcode = another_shiftcode;
			raw_face_ptr->set_state(VALID);
			assert(raw_face_ptr->calculateState() == VALID);

			// Cleanup after modifications is not needed right now. 
			// It may be done after last cut.

		}
		const Vertex* update_vertices_distances(const geometry::Matrix<FloatingPointType>& fractocart) {
			const Vertex* ret = nullptr;
			FloatingPointType m = FloatingPointType(0.0);
			auto s = vertices.size();

			for (size_t i = 0; i < s; i++)
			{
				if (vertices[i]->get_state() == State::DELETE)
					continue;
				auto& curdist = vertices[i]->distance;
				if (curdist == FloatingPointType(0.0)) {
					curdist = (fractocart * (vertices[i]->point - center)).rSq();
				}
				if (curdist > m) {
					m = curdist;
					ret = vertices[i].get();
				}
			}
			assert(ret != nullptr);
			return ret;
		}
	};

	class VoronoiDiagram {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;
		using VoronCell = voronoi::Cell;
		using PointVector = ::std::vector<PointType>;
		using CellVector = ::std::vector<VoronCell>;
		using BoolVector = ::std::vector<bool>;
		using SpatialGrid = geometry::SpatialGrid<FloatingPointType>;
		using Matrix = geometry::Matrix<FloatingPointType>;
		using PointsSorted = std::vector<std::tuple<int, char, FloatingPointType>>;
		using PlaneType = geometry::Plane<FloatingPointType>;

	private:
		//Data
		BoolVector flags_;
		CellVector cells_;
	public:
		explicit VoronoiDiagram(const PointVector& points_in_unit01,
								const std::vector<SpatialGrid::BondWithShift>& bonds,
								const Matrix& FtoC,
								const BoolVector& flags = BoolVector()) : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points_in_unit01.size(), true);
			} else if (points_in_unit01.size() != flags_.size()) {
				flags_.resize(points_in_unit01.size(), false);
			}
			add_points(points_in_unit01, flags_);
			auto vec = find_interactions(bonds);
			calculate_and_sort(vec, points_in_unit01, FtoC);
			for (size_t i = 0; i < cells_.size(); i++)
			{
				if (flags_[i] == false)
					continue;
				manager(cells_[i], vec[i], points_in_unit01, FtoC);
			}
		}

		CellVector extractCells() noexcept {
			return std::move(cells_);
		}
	private:
		void add_points(const PointVector& points, const BoolVector& flags) noexcept {
			cells_.reserve(points.size());
			for (size_t i = 0; i < points.size(); i++)
			{
				if (flags[i]) {
					cells_.emplace_back(points[i], i);
				} else {
					cells_.emplace_back();
				}
			}
		}

		std::vector<PointsSorted> find_interactions(const std::vector<SpatialGrid::BondWithShift>& bonds) noexcept {
			std::vector<PointsSorted> ret(flags_.size());
			for (auto& bond : bonds) {
				// Skip incorrect bonds
				if (bond.first == bond.second)
					continue;


				if (flags_[bond.first] == true) {
					ret[bond.first].emplace_back(bond.second, bond.shiftcode, FloatingPointType(0.0));
				}
				if (flags_[bond.second] == true) {
					ret[bond.second].emplace_back(bond.first, SpatialGrid::inverse_code(bond.shiftcode), FloatingPointType(0.0));
				}
			}
			return ret;
		}
		void calculate_and_sort(std::vector<PointsSorted>& vec, const PointVector& points_in_unit01, const Matrix& FtoC) {
			for (size_t i = 0; i < vec.size(); i++)
			{
				// 1. Calculate distances
				if (flags_[i] == false)
					continue;
				for (auto& [second, code, lengthsq] : vec[i])
				{
					lengthsq = (FtoC * (points_in_unit01[i] -
									  points_in_unit01[second] -
									  SpatialGrid::decompress_shift(code))).rSq();
				}

				// 2. Sort
				std::sort(vec[i].begin(), vec[i].end(),
						  [](const typename PointsSorted::value_type& a,
							 const typename PointsSorted::value_type& b) {
								 return std::get<2>(a) < std::get<2>(b);
						  });
			}
		}
		// NOTE: this function modifies only one cell  
		void manager(VoronCell& cell, const PointsSorted& vec, const PointVector& points_in_unit01, const Matrix& FtoC) const {
			const Vertex* maxVert = cell.update_vertices_distances(FtoC);
			auto maxVertDoubleDistanceSq = maxVert->distance * 4; // Squared double distance
			for (auto& [second, code, lengthsq] : vec) {
				if (maxVert->get_state() == State::DELETE) {
					maxVert = cell.update_vertices_distances(FtoC);
					maxVertDoubleDistanceSq = maxVert->distance * 4;
				}

				if (lengthsq > maxVertDoubleDistanceSq) {
					// Early exit archived
					return;
				}

				// calculate plane
				auto sumsecond = points_in_unit01[second] + SpatialGrid::decompress_shift(code);
				auto inter = (cell.center + sumsecond) * FloatingPointType(0.5);
				auto normal = cell.center - sumsecond;
				normal /= normal.r();

				PlaneType plane(inter, normal);

				cell.clipByPlaneAndAddNewFace(plane, second, code);
			}
		}
	};

	class VoronoiFused {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;

		struct PolygonIn {
			::std::vector<::std::size_t>   vert_ids;
			::std::vector<::std::size_t>   edge_ids;
			::std::array<::std::size_t, 2> atom_ids;

		};
		struct EdgeIn {
			::std::vector<::std::size_t>   atom_ids;
			::std::array<::std::size_t, 2> vert_ids;
		};


		static constexpr FloatingPointType EPSILON = 0.0001;
	public:
		//Data
		::std::vector<PointType> vertices;
		::std::vector<EdgeIn> edges;
		::std::vector<PolygonIn> polygons;
	public:
		explicit VoronoiFused(const ::std::vector<voronoi::Cell>& cells) {

			size_t count_vertices = 0;
			size_t count_edges = 0;
			size_t count_pol = 0;

			for (const auto& cell : cells) {
				count_vertices += cell.vertices.size();
				count_edges += cell.edges.size();
				count_pol += cell.faces.size();
			}
			// reserve memory
			vertices.reserve(count_vertices);
			edges.reserve(count_edges);
			polygons.reserve(count_pol);

			struct SortEntry {
				FloatingPointType key;
				uint32_t vIdx;   // glogal index
				uint32_t lIdx;   // local index
				const Vertex* ptr;
			};
			
			std::vector<SortEntry> sortentries;
			sortentries.reserve(count_vertices);

			// Merge living vertices
			for (const auto& cell : cells) {
				for (const auto& vert : cell.vertices) {
					if (vert->get_state() == State::DELETE)
						continue;
					sortentries.emplace_back(vert->point[0] + vert->point[1] + vert->point[2],
											 sortentries.size(),
											 vert->get_id(),
											 vert.get());
					
				}
			}
			// Sort by key
			std::sort(sortentries.begin(), sortentries.end(), [](auto& a, auto& b) {
				return a.key < b.key;
	        });
			vertices.reserve(sortentries.size());
			
			// Two-eyes comparator
			size_t cur_size = 0;
			size_t right = 0;
			for (size_t left = 0; left < sortentries.size(); left++) {
				for (; right < sortentries.size(); right++) {
					if (sortentries[right].key - sortentries[left].key >= EPSILON)
						break;
				}
				for (size_t iter = left + 1; iter < right; iter++) {
					if (PointType::distanceSq(sortentries[left].ptr->point, sortentries[iter].ptr->point) < EPSILON*EPSILON) {

					}
				}
			}

		}
	};
}
