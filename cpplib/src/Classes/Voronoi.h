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

/// @brief Voronoi diagram construction and analysis namespace
///
/// This namespace contains classes and utilities for constructing and manipulating
/// Voronoi diagrams in 3D space. The Voronoi diagram partitions space based on
/// proximity to a set of input points (typically atom positions).
namespace cpplib::voronoi {
	// Forward declarations
	class Vertex;
	class Edge;
	class Face;
	class Cell;

	/// @brief Container type for non-owning pointers in Voronoi structures
	/// @tparam T Pointer type to store
	template<class T>
	using Container = ::std::unordered_set<T>;

	/// @brief State of a Voronoi geometric object during construction
	///
	/// Objects transition through states during plane clipping operations:
	/// - DELETE: Object is outside the valid region
	/// - VALID: Object is fully valid and unchanged
	/// - INVALID: Object is in an inconsistent state (error condition)
	/// - MODIFICATION: Object is being modified by a clipping operation
	enum class State : char {
		DELETE = 0,       ///< Object should be deleted
		VALID = 1,        ///< Object is valid
		INVALID = 2,      ///< Object is in invalid state
		MODIFICATION = 3  ///< Object is being modified
	};

	/// @brief Base class for Voronoi geometric objects with state and ID management
	///
	/// Provides common functionality for tracking the state and unique identifier
	/// of Voronoi vertices, edges, and faces during diagram construction.
	class Object {
		State state;     ///< Current state of the object
		uint32_t id;     ///< Unique identifier
	public:
		/// @brief Get the unique identifier
		/// @return Object's ID
		inline uint32_t get_id() const noexcept {
			return id;
		}

		/// @brief Set the unique identifier
		/// @param i New ID value
		inline void set_id(uint32_t i) noexcept {
			id = i;
		}

		/// @brief Get the current state
		/// @return Current State
		inline State get_state() const noexcept {
			return state;
		}

		/// @brief Set the current state
		/// @param s New State value
		inline void set_state(State s) noexcept {
			state = s;
		}

		/// @brief Construct an Object with ID and optional state
		/// @param ID Unique identifier
		/// @param s Initial state (default: INVALID)
		constexpr explicit Object(uint32_t ID, State s = State::INVALID) noexcept
			: state(s), id(ID) {}
	};

	/// @brief Epsilon for comparing vertex positions
	///
	/// Vertices within this distance are considered equal to handle
	/// floating-point precision issues.
	constexpr basic_types::FloatingPointType EPSILON = (2<<16) * std::numeric_limits<basic_types::FloatingPointType>::epsilon();

	/// @brief Represents a vertex in a Voronoi diagram
	///
	/// A Vertex stores a 3D point position, distance metric, and maintains
	/// non-owning pointers to connected edges and faces. Vertices are compared
	/// using a spatial epsilon for numerical stability.
	struct Vertex : public Object {
		/// @brief 3D point type for vertex position
		using PointType = geometry::Point<basic_types::FloatingPointType>;


		/// @brief 3D position of this vertex
		PointType point;

		/// @brief Distance metric for this vertex (typically squared distance to cell center)
		basic_types::FloatingPointType distance = basic_types::FloatingPointType(0.0);

		/// @brief Edges connected to this vertex (non-owning pointers)
		Container<Edge*> edges;

		/// @brief Faces connected to this vertex (non-owning pointers)
		Container<Face*> faces;

		/// @brief Default constructor with ID 0
		Vertex() noexcept : Object(0, State::VALID) {}

		/// @brief Construct vertex with given ID
		/// @param ID Unique identifier for this vertex
		explicit Vertex(uint32_t ID) noexcept : Object(ID, State::VALID) {}

		/// @brief Construct vertex with ID and position
		/// @param ID Unique identifier
		/// @param p Position point
		Vertex(uint32_t ID, const PointType& p) noexcept : Object(ID, State::VALID), point(p) {}

		/// @brief Construct vertex with ID and position (move)
		/// @param ID Unique identifier
		/// @param p Position point (moved)
		Vertex(uint32_t ID, PointType&& p) noexcept : Object(ID, State::VALID), point(std::move(p)) {}

		/// @brief Compare vertices for equality using spatial epsilon
		/// @param a First vertex
		/// @param b Second vertex
		/// @return True if vertices are spatially equivalent
		///
		/// Uses EPSILON to handle floating-point precision.
		inline friend bool operator==(const Vertex& a, const Vertex& b) {
			return std::abs(a.point[0] - b.point[0]) < EPSILON &&
				std::abs(a.point[1] - b.point[1]) < EPSILON &&
				std::abs(a.point[2] - b.point[2]) < EPSILON;
		}
	};

	/// @brief Represents an edge in a Voronoi diagram
	///
	/// An Edge connects two vertices and is shared by (at most) two faces.
	/// Edges maintain non-owning pointers to their vertices and adjacent faces.
	struct Edge : public Object {
		/// @brief Vertices at the endpoints of this edge (max 2, non-owning)
		Container<Vertex*> vertices;

		/// @brief Faces adjacent to this edge (max 2, non-owning)
		Container<Face*> faces;

		/// @brief 3D point type
		using PointType = geometry::Point<basic_types::FloatingPointType>;
		/// @brief Plane type for intersection calculations
		using PlaneType = geometry::Plane<basic_types::FloatingPointType>;

		/// @brief Construct an edge connecting two vertices
		/// @param ID Unique identifier
		/// @param v1 First vertex (must not be null)
		/// @param v2 Second vertex (must not be null)
		///
		/// Automatically registers this edge with both vertices.
		Edge(uint32_t ID, Vertex* v1, Vertex* v2) : Object(ID), vertices({v1, v2}) {
			assert(v1 != nullptr);
			assert(v2 != nullptr);
			v1->edges.emplace(this);
			v2->edges.emplace(this);
		}

		/// @brief Calculate and update the state of this edge
		/// @return Updated State
		///
		/// An edge is VALID if it has 2 vertices and 2 faces, with both vertices valid.
		/// An edge is MODIFICATION if it's properly formed but one vertex is being modified.
		/// An edge is DELETE if both vertices are deleted.
		/// Otherwise, it's INVALID.
		inline State calculateState() noexcept {
			using enum State;
			if (vertices.size() == 2 && faces.size() == 2) {
				uint32_t counter = 0;
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

		/// @brief Calculate intersection point of this edge with a plane
		/// @param plane Clipping plane
		/// @return Intersection point
		///
		/// Requires that this edge is in MODIFICATION state (one vertex on each side of plane).
		/// The edge must not be parallel to the plane.
		PointType intersectSegmentPlane(PlaneType plane) {
			assert(calculateState() == State::MODIFICATION);
			auto it = vertices.begin();
			auto v1 = *it;
			it++;
			auto v2 = *it;

			PointType direction = v2->point - v1->point;
			auto unnormalized_normal = PointType(plane.a[0], plane.a[1], plane.a[2]);
			auto denom = PointType::Scalar(unnormalized_normal, direction);

			// Check that Edge is not parallel to plane
			assert(std::abs(denom) >= 1e-6);

			auto t = -(plane.a[0] * v1->point[0] +
					   plane.a[1] * v1->point[1] +
					   plane.a[2] * v1->point[2] +
					   plane.a[3]) / denom;

			assert(t > 0.0 && t < 1.0);

			return v1->point + direction * t;
		}

		/// @brief Get the other vertex of this edge
		/// @param v One vertex of the edge
		/// @return The other vertex
		///
		/// Requires that v is one of the two vertices of this edge.
		Vertex* get_second_vertex(const Vertex* v) const {
			for (auto& i : vertices) {
				if (i != v)
					return i;
			}
			assert(false); // Code should not reach here
			return nullptr;
		}
	};

	/// @brief Represents a face in a Voronoi diagram
	///
	/// A Face is a polygon formed by vertices and edges, separating two Voronoi cells.
	/// Each face has an owner cell and an "other" cell (possibly shifted by periodic boundaries).
	struct Face : public Object {
		/// @brief ID of the cell that owns this face
		uint32_t owner_id;

		/// @brief ID of the neighboring cell on the other side of this face
		uint32_t other_id;

		/// @brief Shift code for periodic boundary conditions
		///
		/// Indicates which periodic image the "other" cell is in.
		geometry::ShiftCode other_shiftcode;

		/// @brief Vertices forming this face (non-owning pointers)
		Container<Vertex*> vertices;

		/// @brief Edges forming this face (non-owning pointers)
		Container<Edge*> edges;

		/// @brief Construct a face with given ID
		/// @param ID Unique identifier
		explicit Face(uint32_t ID) : Object(ID) {}

		/// @brief Calculate and update the state of this face
		/// @return Updated State
		///
		/// A face is VALID if it has equal numbers of vertices and edges, and all edges are valid.
		/// A face is MODIFICATION if any edge is being modified.
		/// A face is DELETE if all edges are deleted.
		/// A face with exactly one valid edge is INVALID (error condition).
		inline State calculateState() {
			using enum State;

			if (vertices.size() == edges.size()) {
				uint32_t counter = 0;
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

	/// @brief Represents a complete Voronoi cell in 3D space
	///
	/// A Cell contains vertices, edges, and faces that define a Voronoi polyhedron.
	/// Supports progressive clipping by planes to construct the final cell geometry.
	/// Cells are initialized as cubes centered at an atomic position and then
	/// iteratively clipped by planes perpendicular to neighboring atoms.
	struct Cell {
	public:
		/// @brief Floating-point type for coordinates
		using FloatingPointType = basic_types::FloatingPointType;
		/// @brief 3D point type
		using PointType = geometry::Point<FloatingPointType>;
		/// @brief Integer shift type for periodic boundaries
		using ShiftType = geometry::Point<int8_t>;
		/// @brief Bond with periodic boundary shift
		using BondWithShift = BondWithPoint<ShiftType>;
		/// @brief Plane type for clipping operations
		using PlaneType = geometry::Plane<FloatingPointType>;

		/// @brief Epsilon for geometric comparisons
		///
		/// Used to determine if a vertex lies exactly on a clipping plane.
		static constexpr FloatingPointType limit = 2 * ::std::numeric_limits<FloatingPointType>::epsilon();

		/// @brief Owned vertices in this cell (smart pointers)
		std::vector<std::unique_ptr<Vertex>> vertices;

		/// @brief Owned edges in this cell (smart pointers)
		std::vector<std::unique_ptr<Edge>> edges;

		/// @brief Owned faces in this cell (smart pointers)
		std::vector<std::unique_ptr<Face>> faces;

		/// @brief Center point of this cell (typically an atom position)
		PointType center;

		/// @brief Cell identifier (typically atom index)
		int id = -1;

		/// @brief Base cube vertices (8 corners, relative to center)
		///
		/// Initial geometry for a Voronoi cell before clipping.
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

		/// @brief Vertex indexes for each face of the cube (counter-clockwise from outside)
		static constexpr std::array<std::array<int, 4>, 6> face_indexes = {{
			{4, 7, 6, 5}, // front face  (+Z)
			{0, 1, 2, 3}, // back face   (-Z)
			{0, 3, 7, 4}, // left face   (-X)
			{1, 5, 6, 2}, // right face  (+X)
			{0, 4, 5, 1}, // bottom face (-Y)
			{3, 2, 6, 7}  // top face    (+Y)
		}};

		/// @brief Shift codes for each face of the base cube
		///
		/// Indicates which periodic image each cube face points toward.
		static constexpr std::array<int8_t, 6> face_shiftcodes = {{
			22, // front face
			 4, // back face
			12, // left face
			14, // right face
			10, // bottom face
			16  // top face
		}};

		/// @brief Vertex indexes for each edge of the cube
		static constexpr std::array<std::array<int, 2>, 12> edge_indexes = {{
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

		/// @brief Edge indexes for each face of the cube
		static constexpr std::array<std::array<int, 4>, 6> face_edge_indexes = {{
			{ 0,  1,  2,  3},  // front face
			{ 4,  5,  6,  7},  // back face
			{ 7,  9,  0,  8},  // left face
			{10,  2, 11,  5},  // right face
			{ 8,  3, 10,  4},  // bottom face
			{ 9,  1, 11,  6}   // top face
		}};

	public:
		/// @brief Default constructor
		Cell() = default;

		/// @brief Construct a Voronoi cell as a cube centered at given point
		/// @param c Center point (typically atom position)
		/// @param i Cell identifier (typically atom index)
		///
		/// Initializes the cell as a cube with 8 vertices, 12 edges, and 6 faces.
		Cell(const PointType& c, int i) : center(c), id(i) {
			vertices.reserve(128);
			edges.reserve(128);
			faces.reserve(64);

			// Create base vertices
			for (uint32_t vert_id = 0; vert_id < 8; vert_id++) {
				vertices.emplace_back(std::make_unique<Vertex>(vert_id, base_vertices[vert_id] + center));
			}

			// Create edges connecting vertices
			for (uint32_t edge_id = 0; edge_id < 12; edge_id++) {
				auto vertex1_ptr = vertices[edge_indexes[edge_id][0]].get();
				auto vertex2_ptr = vertices[edge_indexes[edge_id][1]].get();
				auto owner_ptr = std::make_unique<Edge>(edge_id, vertex1_ptr, vertex2_ptr);
				auto edge_ptr = owner_ptr.get();
				edges.emplace_back(std::move(owner_ptr));
				edge_ptr->set_state(State::VALID);
			}

			// Create faces
			for (uint32_t face_id = 0; face_id < 6; face_id++) {
				auto face = std::make_unique<Face>(face_id);
				face->owner_id = id;
				face->other_id = id;
				face->other_shiftcode = face_shiftcodes[face_id];

				// Fill vertices for this face
				for (int j = 0; j < 4; j++) {
					auto vertex_ptr = vertices[face_indexes[face_id][j]].get();
					face->vertices.emplace(vertex_ptr);
					vertex_ptr->faces.emplace(face.get());
				}

				// Fill edges for this face
				for (int j = 0; j < 4; j++) {
					auto edge_ptr = edges[face_edge_indexes[face_id][j]].get();
					face->edges.emplace(edge_ptr);
					edge_ptr->faces.emplace(face.get());
				}

				face->set_state(State::VALID);
				faces.push_back(std::move(face));
			}
		}

		/// @brief Add a new vertex to the cell (with deduplication)
		/// @param p Vertex position
		/// @return Pointer to the vertex (new or existing)
		///
		/// If a vertex at this position already exists (within epsilon), returns
		/// the existing vertex. Otherwise, creates a new vertex.
		Vertex* add_vertex(const PointType& p) {
			auto temp = std::make_unique<Vertex>(static_cast<uint32_t>(vertices.size()), p);
			for (auto& v : vertices) {
				if (*v == *temp)
					return v.get();
			}
			auto simple_ptr = temp.get();
			vertices.emplace_back(std::move(temp));
			return simple_ptr;
		}

		/// @brief Clip cell by a plane and add the resulting face
		/// @param clipping_plane Plane to clip by
		/// @param id_of_another_cell ID of the neighboring cell
		/// @param another_shiftcode Shift code for periodic boundaries
		///
		/// Progressively clips the Voronoi cell by a plane perpendicular to a neighboring atom.
		/// This algorithm:
		/// 
		/// 1. Classifies vertices as deleted, modified, or valid based on their side of the plane
		/// 
		/// 2. Updates edge and face states accordingly
		/// 
		/// 3. Creates new vertices at edge-plane intersections
		/// 
		/// 4. Adds new edges and a new face where the plane cuts the cell
		void clipByPlaneAndAddNewFace(const PlaneType& clipping_plane, uint32_t id_of_another_cell, geometry::ShiftCode another_shiftcode) {
			using enum State;

			// Check if the side is correct (center must be on positive side)
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
				auto side = clipping_plane.side(v->point);
				if (side < -limit) {
					// Point cutted off
					v->set_state(DELETE);
					modified = true;
				} else if (side < limit) {
					// Point on the Face
					v->set_state(MODIFICATION);
				}
			}

			// Early exit if no vertices were cut
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
			// Check all Edges which need MODIFICATION
			for (auto& e : edges)
			{
				if (e->get_state() != MODIFICATION) {
					continue;
				}
				auto intersection = e->intersectSegmentPlane(clipping_plane);
				auto new_vertex_ptr = add_vertex(intersection);

				// Erase deleted vertex from set
				std::erase_if(e->vertices,
							  [](auto* ptr) {
								  if (ptr->get_state() == DELETE) {
									  return true;
								  }
								  return false;
							  });
				// Add new vertex to set
				e->vertices.insert(new_vertex_ptr);

				new_vertex_ptr->edges.insert(e.get());
				new_vertex_ptr->faces.insert(e->faces.begin(), e->faces.end());
				for (auto& f : e->faces) {
					f->vertices.emplace(new_vertex_ptr);
				}
				new_vertex_ptr->set_state(MODIFICATION); // Mark as on the NEW FACE

				// Now, new vertex complete, so edge is VALID too:
				e->set_state(VALID);
			}

			// 4. Cut Faces which need modification
			for (auto& f : faces) {
				if (f->get_state() != MODIFICATION) continue;

				// 4.1. Find and delete all unnecessary edges and vertices
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

				// Find new vertices (those on the clipping plane)
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
				if(iter_vertex == f->vertices.cend())
					assert(iter_vertex != f->vertices.cend());

				// Create edge connecting the two new vertices
				const auto& new_edge = edges.emplace_back(std::make_unique<Edge>(static_cast<uint32_t>(edges.size()), v1, v2));
				new_edge->faces.emplace(f.get());
				f->edges.emplace(new_edge.get());
				new_edge->set_state(MODIFICATION);
				f->set_state(VALID);
			}

			// 5. Create new Face (the clipping plane becomes a face)
			const auto& new_face = faces.emplace_back(std::make_unique<Face>(static_cast<uint32_t>(faces.size())));
			auto raw_face_ptr = new_face.get();

			// Add all MODIFICATION edges to the new face
			for (const auto& e : edges)
			{
				if (e->get_state() == MODIFICATION) {
					raw_face_ptr->edges.emplace(e.get());
					e->faces.emplace(raw_face_ptr);
					e->set_state(VALID);
					assert(e->calculateState() == VALID);
				}
			}

			// Add all MODIFICATION vertices to the new face
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
			raw_face_ptr->other_shiftcode = another_shiftcode.get_code();
			raw_face_ptr->set_state(VALID);
			assert(raw_face_ptr->calculateState() == VALID);

			// Cleanup after modifications is not needed right now. 
			// It may be done after last cut.
		}

		/// @brief Update vertex distances from the cell center
		/// @param fractocart Transformation matrix (fractional to Cartesian)
		/// @return Pointer to the vertex with maximum distance
		///
		/// Calculates squared distances for vertices that haven't been computed yet,
		/// and returns the vertex farthest from the center. Used for early exit
		/// optimization during cell construction.
		const Vertex* update_vertices_distances(const geometry::Matrix<FloatingPointType>& fractocart) {
			const Vertex* ret = nullptr;
			FloatingPointType m = FloatingPointType(0.0);
			auto s = static_cast<uint32_t>(vertices.size());

			for (uint32_t i = 0; i < s; i++)
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

	/// @brief Main class for constructing Voronoi diagrams from atomic positions
	///
	/// VoronoiDiagram constructs Voronoi cells for a set of input points (atoms) by:
	/// 
	/// 1. Initializing each cell as a cube centered at the point
	/// 
	/// 2. Finding neighboring points via a spatial grid and bonds
	/// 
	/// 3. Iteratively clipping each cell by planes perpendicular to its neighbors
	/// 
	/// 4. Extracting the final cell geometry
	///
	/// Note: Only cells marked with flags=true are fully computed. Others are left empty.
	class VoronoiDiagram {
	public:
		/// @brief Floating-point type
		using FloatingPointType = basic_types::FloatingPointType;
		/// @brief 3D point type
		using PointType = geometry::Point<FloatingPointType>;
		/// @brief Voronoi cell type
		using VoronCell = voronoi::Cell;
		/// @brief Vector of points
		using PointVector = ::std::vector<PointType>;
		/// @brief Vector of cells
		using CellVector = ::std::vector<VoronCell>;
		/// @brief Vector of flags
		using BoolVector = ::std::vector<bool>;
		/// @brief Spatial grid for neighbor finding
		using SpatialGrid = geometry::SpatialGrid<FloatingPointType>;
		/// @brief Transformation matrix type
		using Matrix = geometry::Matrix<FloatingPointType>;
		/// @brief Sorted list of neighbor interactions (index, shiftcode, distance?)
		using PointsSorted = std::vector<std::tuple<int, geometry::ShiftCode, FloatingPointType>>;
		/// @brief Plane type
		using PlaneType = geometry::Plane<FloatingPointType>;

	private:
		//Data
		/// @brief Flags indicating which cells to compute
		BoolVector flags_;
		/// @brief Computed Voronoi cells
		CellVector cells_;

	public:
		/// @brief Construct a Voronoi diagram from points and bonds
		/// @param points_in_unit01 Atomic positions in fractional coordinates [0,1]
		/// @param bonds List of bonds with periodic shifts
		/// @param FtoC Transformation matrix (fractional to Cartesian coordinates)
		/// @param flags Optional flags indicating which cells to compute (default: all true)
		///
		/// The constructor performs the complete Voronoi construction:
		/// - Initializes cells for all points
		/// - Finds interactions via bonds
		/// - Sorts neighbors by distance
		/// - Clips each cell by planes to neighboring atoms
		///
		/// Only cells with flags[i]=true are fully computed. This allows computing only
		/// the asymmetric unit cells while using symmetry-expanded points for boundaries.
		explicit VoronoiDiagram(const PointVector& points_in_unit01,
								const std::vector<SpatialGrid::BondWithShift>& bonds,
								const geometry::Cell<FloatingPointType>& unitcell,
								const BoolVector& flags = BoolVector()) : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points_in_unit01.size(), true);
			} else if (points_in_unit01.size() != flags_.size()) {
				flags_.resize(points_in_unit01.size(), false);
			}
			add_points(points_in_unit01, flags_);
			auto vec = find_interactions(bonds);
			calculate_and_sort(vec, points_in_unit01, unitcell.fracToCart());
			for (uint32_t i = 0; i < static_cast<uint32_t>(cells_.size()); i++)
			{
				if (flags_[i] == false)
					continue;
				manager(cells_[i], vec[i], points_in_unit01, unitcell);
			}
		}

		/// @brief Extract the computed Voronoi cells
		/// @return Vector of cells (moved out)
		///
		/// After calling this, the diagram is left in an empty state.
		CellVector extractCells() noexcept {
			return std::move(cells_);
		}

	private:
		/// @brief Initialize cells for all points
		/// @param points Atomic positions
		/// @param flags Flags indicating which cells to compute
		void add_points(const PointVector& points, const BoolVector& flags) noexcept {
			cells_.reserve(points.size());
			for (uint32_t i = 0; i < static_cast<uint32_t>(points.size()); i++)
			{
				if (flags[i]) {
					cells_.emplace_back(points[i], i);
				} else {
					cells_.emplace_back();
				}
			}
		}

		/// @brief Find neighbor interactions from bonds
		/// @param bonds List of bonds with shift codes
		/// @return Vector of neighbor lists for each point
		///
		/// For each cell, creates a list of (neighbor_id, shiftcode, distance?) tuples.
		std::vector<PointsSorted> find_interactions(const std::vector<SpatialGrid::BondWithShift>& bonds) noexcept {
			std::vector<PointsSorted> ret(flags_.size());
			for (auto& bond : bonds) {
				// Skip incorrect bonds (self-bonds)
				if (bond.first == bond.second)
					continue;

				// Add neighbor to first atom's list
				if (flags_[bond.first] == true) {
					ret[bond.first].emplace_back(bond.second, bond.shiftcode, FloatingPointType(0.0));
				}
				// Add neighbor to second atom's list (with inverse shift)
				if (flags_[bond.second] == true) {
					ret[bond.second].emplace_back(bond.first, geometry::ShiftCode::inverse(bond.shiftcode), FloatingPointType(0.0));
				}
			}
			return ret;
		}

		/// @brief Calculate distances and sort neighbors for each cell
		/// @param vec Neighbor lists (modified in-place)
		/// @param points_in_unit01 Atomic positions
		/// @param FtoC Fractional to Cartesian transformation
		///
		/// Computes squared distances to neighbors in Cartesian space and sorts
		/// neighbors by increasing distance for each cell.
		void calculate_and_sort(std::vector<PointsSorted>& vec, const PointVector& points_in_unit01, const Matrix& FtoC) {
			auto vec_s = static_cast<uint32_t>(vec.size());
			for (uint32_t i = 0; i < vec_s; i++)
			{
				// 1. Calculate distances
				if (flags_[i] == false)
					continue;
				for (auto& [second, code, lengthsq] : vec[i])
				{
					lengthsq = (FtoC * (points_in_unit01[i] -
										points_in_unit01[second] -
										code.get_shift())).rSq();
				}

				// 2. Sort by distance
				std::sort(vec[i].begin(), vec[i].end(),
						  [](const typename PointsSorted::value_type& a,
							 const typename PointsSorted::value_type& b) {
								 return std::get<2>(a) < std::get<2>(b);
						  });
			}
		}

		/// @brief Construct a single Voronoi cell by clipping
		/// @param cell Cell to modify
		/// @param vec Sorted neighbor list
		/// @param points_in_unit01 Atomic positions
		/// @param FtoC Fractional to Cartesian transformation
		///
		/// Iteratively clips the cell by planes perpendicular to neighbors (sorted by distance).
		/// Uses early exit optimization: stops when the nearest unprocessed neighbor is farther
		/// than twice the distance to the farthest vertex of the current cell.
		void manager(VoronCell& cell, const PointsSorted& vec, const PointVector& points_in_unit01, const geometry::Cell<FloatingPointType>& unitcell) const {
			const Vertex* maxVert = cell.update_vertices_distances(unitcell.fracToCart());
			auto maxVertDoubleDistanceSq = maxVert->distance * 4; // Squared double distance

			for (auto& [second, code, lengthsq] : vec) {
				// Update max vertex if it was deleted
				if (maxVert->get_state() == State::DELETE) {
					maxVert = cell.update_vertices_distances(unitcell.fracToCart());
					maxVertDoubleDistanceSq = maxVert->distance * 4;
				}

				// Early exit: neighbor is too far to affect the cell
				if (lengthsq > maxVertDoubleDistanceSq) {
					return;
				}

				// Calculate clipping plane
				auto sumsecond = points_in_unit01[second] + code.get_shift();

				auto cartA = unitcell.fracToCart() * cell.center;
				auto cartB = unitcell.fracToCart() * sumsecond;

				auto normal_cart = cartA - cartB;

				auto inter_frac = (cell.center + sumsecond) * FloatingPointType(0.5); // Midpoint
				auto hkl = unitcell.fracToCart().TransposeMultiply(normal_cart);

				PlaneType plane(inter_frac, hkl);


				auto inter_cart = (cartA + cartB) * FloatingPointType(0.5);

				PlaneType plane_cart(inter_cart, normal_cart);

				PointType b = inter_cart;
				PointType c = inter_cart;

				auto nx = plane_cart.a[0];
				auto ny = plane_cart.a[1];
				auto nz = plane_cart.a[2];

				if (std::abs(nx) > EPSILON) {
					b[1] += 1.0;
					b[0] -= ny / nx;

					c[2] += 1.0;
					c[0] -= nz / nx;
				} else if (std::abs(ny) > 1e-9) {
					b[0] += 1.0;

					c[2] += 1.0; 
					c[1] -= nz / ny; 
				} else {
					b[0] += 1.0;
					c[1] += 1.0;
				}

				PlaneType plane_other(inter_frac, unitcell.cartToFrac() * b, unitcell.cartToFrac() * c);
				if (plane_other.side(cell.center) < 0) {
					plane_other.a[0] = -plane_other.a[0];
					plane_other.a[1] = -plane_other.a[1];
					plane_other.a[2] = -plane_other.a[2];
					plane_other.a[3] = -plane_other.a[3];
				}
				
				cell.clipByPlaneAndAddNewFace(plane_other, second, code);
			}
		}
	};

	/// @brief Fused Voronoi structure for merged/unified cells
	///
	/// VoronoiFused takes a collection of Voronoi cells and merges coincident vertices,
	/// edges, and faces to create a unified polyhedral structure. This is useful for
	/// visualization and further analysis of the Voronoi tessellation.
	class VoronoiFused {
	public:
		/// @brief Floating-point type
		using FloatingPointType = basic_types::FloatingPointType;
		/// @brief 3D point type
		using PointType = geometry::Point<FloatingPointType>;

		/// @brief Polygon (face) in the fused structure
		struct PolygonIn {
			geometry::ShiftCode::ShiftPoint second_shift; ///< Shift code of the second atom
			FloatingPointType area;                       ///< Area of Polygon
			FloatingPointType solidangle;                 ///< Solid angle of Polygon
			::std::vector<uint32_t>   vert_ids;           ///< Vertex indexes forming the polygon
			::std::vector<uint32_t>   edge_ids;           ///< Edge indexes forming the polygon
			::std::array<uint32_t, 2> atom_ids;           ///< IDs of the two atoms separated by this face
		};

		/// @brief Edge in the fused structure
		struct EdgeIn {
			::std::array<uint32_t, 2> vert_ids;   ///< Vertex indexes at endpoints
		};

		/// @brief Polyhedron (cell) in the fused structure
		struct Polyhedron {
			FloatingPointType volume;
			::std::vector<uint32_t>   vert_ids;   ///< Vertex indexes in this polyhedron
			::std::vector<uint32_t>   edge_ids;   ///< Edge indexes in this polyhedron
			::std::vector<uint32_t>   poly_ids;   ///< Polygon indexes in this polyhedron
			PointType center;
		};

		/// @brief Helper structure for sorting and merging vertices
		struct SortEntry {
			bool is_merged;                  ///< Whether this vertex was merged with others
			std::vector<uint32_t> cIdx;      ///< Cell indexes this vertex belongs to
			FloatingPointType key;           ///< Sort key (sum of coordinates)
			Vertex* ptr;                     ///< Pointer to original vertex
		};

	public:
		//Data
		/// @brief Unified vertex positions
		::std::vector<PointType> vertices;
		/// @brief Unified edges
		::std::vector<EdgeIn> edges;
		/// @brief Unified polygons (faces)
		::std::vector<PolygonIn> polygons;
		/// @brief Polyhedra (cells)
		::std::vector<Polyhedron> polyhedra;

	public:
		/// @brief Construct a fused Voronoi structure from cells
		/// @param cells Vector of Voronoi cells to merge
		///
		/// Performs the following operations:
		/// 
		/// 1. Merges coincident vertices across cells
		/// 
		/// 2. Merges duplicate edges
		/// 
		/// 3. Merges duplicate faces (polygons)
		/// 
		/// 4. Builds unified polyhedra (cells) referencing the merged geometry
		explicit VoronoiFused(::std::vector<voronoi::Cell>& cells, const geometry::Matrix<FloatingPointType>& FtoC) {
			uint32_t count_vertices = 0;
			polyhedra.resize(cells.size());
			for (uint32_t i = 0; i < cells.size(); i++) {
				polyhedra[i].center = cells[i].center;
			}

			for (const auto& cell : cells) {
				count_vertices += cell.vertices.size();
			}


			std::vector<SortEntry> sortentries;
			sortentries.reserve(count_vertices);

			// Collect all living vertices
			for (const auto& cell : cells) {
				for (const auto& vert : cell.vertices) {
					if (vert->get_state() == State::DELETE)
						continue;
					sortentries.emplace_back(false,
											 std::vector<uint32_t>(1, cell.id),
											 vert->point[0] + vert->point[1] + vert->point[2],
											 vert.get());
					sortentries.back().cIdx.reserve(cells.size());
				}
			}

			// Sort by coordinate sum (for efficient proximity search)
			std::sort(sortentries.begin(), sortentries.end(), [](auto& a, auto& b) {
				return a.key < b.key;
			});

			// Merge coincident vertices using two-pointer technique
			count_vertices = static_cast<uint32_t>(sortentries.size());
			uint32_t right = 0;
			for (uint32_t left = 0; left < count_vertices; left++) {
				if (sortentries[left].ptr == nullptr)
					continue;
				// Advance right pointer to include all vertices within EPSILON of left
				for (; right < count_vertices; right++) {
					if (sortentries[right].key - sortentries[left].key >= EPSILON)
						break;
				}
				// Check for exact matches within the window
				for (uint32_t iter = left + 1; iter < right; iter++) {
					if (sortentries[iter].ptr == nullptr)
						continue;
					if (PointType::isSame(sortentries[left].ptr->point, sortentries[iter].ptr->point, EPSILON)) {
						// Merge iter into left
						unite_vertices(sortentries[left].ptr, sortentries[iter].ptr);
						sortentries[iter].ptr = nullptr;
						sortentries[left].is_merged = true;
						sortentries[left].cIdx.insert(sortentries[left].cIdx.end(),
													  sortentries[iter].cIdx.begin(),
													  sortentries[iter].cIdx.end());
						sortentries[iter].cIdx.clear();
					}
				}
			}

			// Remove deleted entries and finalize vertices
			std::erase_if(sortentries, [](const auto& entry) {
				return entry.ptr == nullptr;
			});
			count_vertices = static_cast<uint32_t>(sortentries.size());
			vertices.reserve(sortentries.size());
			for (uint32_t i = 0; i < count_vertices; ++i) {
				sortentries[i].ptr->set_id(i);
				vertices.emplace_back(sortentries[i].ptr->point);
				std::sort(sortentries[i].cIdx.begin(), sortentries[i].cIdx.end());
				for (auto& c : sortentries[i].cIdx) {
					polyhedra[c].vert_ids.push_back(i);
				}
			}

			// Merge edges (find and unite duplicates)
			uint32_t count_edges = 0;
			for (uint32_t i = 0; i < count_vertices; i++)
			{
				if (sortentries[i].is_merged == true) {
					count_edges += find_dublicate_and_count_edges(sortentries[i].ptr);
				} else {
					count_edges += sortentries[i].ptr->edges.size();
				}
			}
			count_edges >>= 1; // Divide by 2 (each edge counted twice)
			edges.reserve(count_edges);

			// Finalize edges
			for (const auto& cell : cells) {
				for (const auto& edge : cell.edges) {
					if (edge->get_state() == State::DELETE)
						continue;
					edge->set_id(edges.size());
					auto v1 = *(edge->vertices.begin());
					auto v2 = edge->get_second_vertex(v1);
					add_edge_to_polyhedra(v1->get_id(), v2->get_id(), edges.size(), sortentries);
					edges.emplace_back(EdgeIn{{v1->get_id(), v2->get_id()}});
				}
			}

			// Merge polygons (faces)
			for (const auto& cell : cells) {
				for (const auto& face : cell.faces) {
					if (face->get_state() == State::DELETE)
						continue;
					// Skip duplicate internal faces (keep only one copy)
					if (cells[face->other_id].vertices.empty() == false &&
					   face->other_shiftcode.get_code() == 13 &&
					   face->owner_id > face->other_id) {
						continue;
					}

					// Create polygon entry
					polyhedra[face->owner_id].poly_ids.push_back(polygons.size());
					if (cells[face->other_id].vertices.empty() == false &&
					    face->other_shiftcode.get_code() == 13) {
						polyhedra[face->other_id].poly_ids.push_back(polygons.size());
					}
					polygons.emplace_back();
					auto& cur_poly = polygons.back();
					cur_poly.vert_ids.resize(face->vertices.size());
					cur_poly.edge_ids.reserve(face->edges.size());
					cur_poly.atom_ids = {face->owner_id, face->other_id};
					for (const auto& edge : face->edges) {
						cur_poly.edge_ids.push_back(edge->get_id());
					}
					reorder_vertices_and_edges_in_polygon(cur_poly);
					cur_poly.second_shift = face->other_shiftcode.get_shift();

					// Calculate area and solid angle
					cur_poly.area = calculate_area(cur_poly, FtoC);
					cur_poly.solidangle = calculate_solid_angle(cur_poly, FtoC);
				}
			}
			// Calculate volumes
			for (auto& p : polyhedra) {
				p.volume = calculate_volume(p, FtoC);
			}
		}

		/// @brief Merge two vertices into one
		/// @param a Target vertex (will contain merged data)
		/// @param b Source vertex (will be invalidated)
		///
		/// Updates all edges and faces referencing b to reference a instead,
		/// then copies all connectivity information from b to a.
		void unite_vertices(Vertex* a, Vertex* b) const {
			// Check container type dependency
			static_assert(std::is_same_v<cpplib::voronoi::Container<Vertex*>, std::unordered_set<Vertex*>>,
						  "Method was written for case when Container == unordered_set. Rewrite method elsewhere.");
			if (a == b) {
				assert(false); // Should be unreachable
				return;
			}

			// Update refs in edges
			for (auto& edge : b->edges) {
				edge->vertices.erase(b);
				edge->vertices.insert(a);
			}
			// Update refs in faces
			for (auto& face : b->faces) {
				face->vertices.erase(b);
				face->vertices.insert(a);
			}
			// Copy connectivity from b to a
			a->edges.insert(b->edges.cbegin(), b->edges.cend());
			a->faces.insert(b->faces.cbegin(), b->faces.cend());
		}

		/// @brief Merge two edges into one
		/// @param a Target edge (will contain merged data)
		/// @param b Source edge (will be marked for deletion)
		///
		/// Updates all faces referencing b to reference a instead,
		/// then copies face connectivity from b to a. Vertex merging is not needed.
		void unite_edges(Edge* a, Edge* b) const {
			// Check container type dependency
			static_assert(std::is_same_v<cpplib::voronoi::Container<Vertex*>, std::unordered_set<Vertex*>>,
						  "Method was written for case when Container == unordered_set. Rewrite method elsewhere.");
			if (a == b) {
				assert(false); // Should be unreachable
				return;
			}
			// Vertex copying is not necessary (vertices already merged)

			// Update refs in faces
			for (auto& face : b->faces) {
				face->edges.erase(b);
				face->edges.insert(a);
			}
			a->faces.insert(b->faces.cbegin(), b->faces.cend());
			b->set_state(State::DELETE);
		}

		/// @brief Find and merge duplicate edges connected to a vertex
		/// @param v Vertex to check
		/// @return Count of unique edges after merging
		///
		/// Searches for edges connected to v that have the same endpoints
		/// (i.e., duplicate edges) and merges them.
		uint32_t find_dublicate_and_count_edges(const Vertex* v) const {
			uint32_t count = v->edges.size();
			auto left = v->edges.cbegin();
			auto end = v->edges.cend();
			for (; left != end; left++) {
				auto cur_left = *left;
				if (cur_left->get_state() == State::DELETE) {
					continue;
				}
				auto second1 = cur_left->get_second_vertex(v)->get_id();
				// Check all remaining edges
				auto right = left;
				right++;
				for (; right != end; right++) {
					auto cur_right = *right;
					if (cur_right->get_state() == State::DELETE) {
						continue;
					}
					auto second2 = cur_right->get_second_vertex(v)->get_id();
					if (second1 == second2) {
						// Found duplicate - merge them
						unite_edges(*left, *right);
						count--;
					}
				}
			}
			return count;
		}


		void add_edge_to_polyhedra(uint32_t a, uint32_t b, uint32_t edge, const std::vector<SortEntry>& sort_entries) {
			uint32_t i1 = 0;
			uint32_t i2 = 0;
			uint32_t s1 = sort_entries[a].cIdx.size();
			uint32_t s2 = sort_entries[b].cIdx.size();

			while (i1 < s1 && i2 < s2) {
				auto v1 = sort_entries[a].cIdx[i1];
				auto v2 = sort_entries[b].cIdx[i2];
				if (v1 == v2) {
					polyhedra[v1].edge_ids.push_back(edge);
					i1++;
					i2++;
				}
				else if (v1 < v2) {
					i1++;
				}
				else {
					i2++;
				}
			}
		}

	private:
		void reorder_vertices_and_edges_in_polygon(PolygonIn& p) {
			assert(p.edge_ids.size() >= 3);
			assert(p.edge_ids.size() == p.vert_ids.size());

			auto current_edge = (edges[p.edge_ids[0]].vert_ids);
			const auto base_vertex = current_edge[0];
			auto next_vertex = current_edge[1];
			const auto s = static_cast<uint32_t>(p.edge_ids.size());
			const auto s1 = s - 1;

			p.vert_ids[0] = base_vertex;
			for (uint32_t i = 1; i < s1; i++)
			{
				assert(next_vertex != base_vertex);

				p.vert_ids[i] = next_vertex;

				for (uint32_t j = i; j < s; j++) {
					auto v1 = edges[p.edge_ids[j]].vert_ids[0];
					auto v2 = edges[p.edge_ids[j]].vert_ids[1];
					if (v1 != next_vertex && v2 != next_vertex) {
						continue;
					}
					if (v1 == next_vertex) {
						next_vertex = v2;
					} else {
						next_vertex = v1;
					}
					std::swap(p.edge_ids[i], p.edge_ids[j]);
					break;
				}
			}
			p.vert_ids[s1] = next_vertex;
		}

		// Requires correct order of vertices in polygon
		FloatingPointType calculate_area(const PolygonIn& p, const geometry::Matrix<FloatingPointType>& FtoC) const {
			PointType center(0, 0, 0);

			for (auto v : p.vert_ids) {
				center += vertices[v];
			}
			center /= p.vert_ids.size();

			std::vector<PointType> cart_verts;
			cart_verts.reserve(p.vert_ids.size());

			for (auto v : p.vert_ids) {
				cart_verts.emplace_back(FtoC * (vertices[v] - center));
			}

			FloatingPointType area = PointType::Vector(cart_verts.front(), cart_verts.back()).r() ;
			for (uint32_t i = 1; i < p.vert_ids.size(); i++) {
				area += PointType::Vector(cart_verts[i], cart_verts[i - 1]).r();
			}

			return area * FloatingPointType(0.5);
		}

		// Requires correct order of vertices in polygon
		FloatingPointType calculate_solid_angle(const PolygonIn& p, const geometry::Matrix<FloatingPointType>& FtoC) const {

			const auto& realO = polyhedra[p.atom_ids[0]].center;

			std::vector<PointType> cart_verts;
			std::vector<FloatingPointType> cart_verts_r;
			const auto cart_size = p.vert_ids.size();
			cart_verts.reserve(cart_size);
			cart_verts_r.reserve(cart_size);

			for (uint32_t i = 0; i < cart_size; i++) {
				cart_verts.emplace_back(FtoC * (vertices[p.vert_ids[i]] - realO));
				cart_verts_r.emplace_back(cart_verts.back().r());
			}

			FloatingPointType angle = 0;
			for (uint32_t i = 2; i < cart_size; i++) {
				angle += atan2(abs(PointType::Scalar(cart_verts[0], PointType::Vector(cart_verts[i - 1], cart_verts[i]))),
							   cart_verts_r[0]* cart_verts_r[i-1] * cart_verts_r[i] +
							   PointType::Scalar(cart_verts[0], cart_verts[i - 1]) * cart_verts_r[i] +
							   PointType::Scalar(cart_verts[0], cart_verts[i]) * cart_verts_r[i - 1] +
							   PointType::Scalar(cart_verts[i - 1], cart_verts[i]) * cart_verts_r[0]);
			}

			return angle*2;
		}

		// Requires area to be calculated in polygons
		FloatingPointType calculate_volume(const Polyhedron& p, const geometry::Matrix<FloatingPointType>& FtoC) const {
			FloatingPointType volume = 0;

			for (auto poly : p.poly_ids) {
				auto i1 = polygons[poly].vert_ids[0];
				auto i2 = polygons[poly].vert_ids[1];
				auto i3 = polygons[poly].vert_ids[2];
				geometry::Plane plane(FtoC * vertices[i1], FtoC * vertices[i2], FtoC * vertices[i3]);
				volume += plane.distance(FtoC * p.center) * polygons[poly].area * FloatingPointType(1./3);
			}

			return volume;
		}
	};
	class Net {
    public:
		using ShiftType = geometry::Point<uint8_t>;
		using PointType = VoronoiFused::PointType;
		using PolyhedronPointer = typename VoronoiFused::Polyhedron*;

	private:
		std::vector<VoronoiFused::Polyhedron*> root_polyhedron;
		std::vector<ShiftType> shift;
		std::vector<PointType> real_point;
		std::vector<size_t> all_neighbours;
		std::vector<size_t> offsets;

	public:
		Net(size_t size) {
			root_polyhedron.reserve(size);
			shift.reserve(size);
			real_point.reserve(size);
			all_neighbours.reserve(size << 4);
			offsets.reserve(size);
		}
		void add_polyhedron(PolyhedronPointer p, ShiftType s, const PointType& r) {
			root_polyhedron.push_back(p);
			shift.push_back(s);
			real_point.push_back(r);
		}

	private:
		size_t pack_key(ShiftType shift, uint32_t atom_id) const noexcept {
			return (static_cast<size_t>(shift[0]))       |
				   (static_cast<size_t>(shift[1]) << 8)  |
				   (static_cast<size_t>(shift[2]) << 16) |
				   (static_cast<size_t>(atom_id) << 24);
		}
		std::pair<ShiftType, uint32_t> unpack_key(size_t key) const noexcept {
			return std::make_pair(ShiftType(static_cast<uint8_t>(key & 0xFF),
											static_cast<uint8_t>((key >> 8) & 0xFF),
											static_cast<uint8_t>((key >> 16) & 0xFF)),
								  static_cast<uint32_t>((key >> 24) & 0xFFFFFFFF));
		}



	};

}