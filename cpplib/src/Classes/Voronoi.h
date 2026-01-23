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
#include <cmath>
#include <cstdint>
#include <memory>
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

	struct Vertex {
		// Types
		using PointType = geometry::Point<basic_types::FloatingPointType>;


		// Constants
		static constexpr typename PointType::value_type COMPARISON_EPSILON = 1.0 / (1 << 16);

		// Data
		State state = State::VALID;
		PointType point;

		// Data (not owning)
		Container<Edge*> edges;
		Container<Face*> faces;

		// Constructors
		Vertex() = default;
		explicit Vertex(const PointType& p) noexcept : point(p) {}
		explicit Vertex(PointType&& p) noexcept : point(std::move(p)) {}

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
		inline State get_state() const noexcept {
			return state;
		}
		inline void set_state(State s) noexcept {
			state = s;
		}

	};

	struct Edge {
		// Data (not owning)
		State state = State::INVALID;
		Container<Vertex*> vertexes; // max 2
		Container<Face*> faces; // max 2

		// Types
		using PointType = geometry::Point<basic_types::FloatingPointType>;
		using PlaneType = geometry::Plane<basic_types::FloatingPointType>;

		inline State calculateState() noexcept {
			const auto vs = vertexes.size();
			const auto fs = faces.size();
			if (vs == 2 && fs == 2) {
				size_t counter = 0;
				for (const auto v : vertexes)
				{
					auto curstate = v->get_state();
					if (curstate == State::VALID || curstate == State::MODIFICATION)
						counter++;
				}
				switch (counter) {
					case 2:
						state = State::VALID;
						break;
					case 1:
						state = State::MODIFICATION;
						break;
					case 0:
						state = State::DELETE;
						break;
					default:
						// Impossible
						break;
				}
			}
			else {
				state = State::INVALID;
			}
			return state;
		}
		inline State get_state() const noexcept {
			return state;
		}
		inline void set_state(State s) noexcept {
			state = s;
		}

		PointType intersectSegmentPlane(PlaneType plane) {
			assert(calculateState() == State::MODIFICATION);
			auto it = vertexes.begin();
			auto v1 = *it;
			it++;
			auto v2 = *it;


			PointType direction = v2->get_point() - v1->get_point();
			auto normal = plane.normal();
			auto denom = PointType::Scalar(normal, direction);

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

	struct Face {
		// Data (not owning)
		State state = State::INVALID;
		size_t id_owner;
		size_t id_other;
		Container<Vertex*> vertexes; 
		Container<Edge*> edges; 

		inline void calculateState() {
			const auto vs = vertexes.size();
			const auto es = edges.size();
			if (vs ==  es) {
				size_t counter = 0;
				for (const auto e : edges)
				{
					auto estate = e->get_state();
					if (estate == State::MODIFICATION || estate == State::VALID)
						counter++;
				}
				if (counter == es) {  // All edges are VALID
					state = State::VALID;
				}
				else if(counter >= 2) { // At least 2 edges are VALID or ask for MODIFICATION
					state = State::MODIFICATION;
				}
				else if (counter == 1) { // ERROR STATE!!!
					state = State::INVALID;
				} 
				else {
					state = State::DELETE;
				}
			} else {
				state = State::INVALID;
			}
		}
		inline State get_state() const noexcept {
			return state;
		}
		inline void set_state(State s) noexcept {
			state = s;
		}
	};
	class Cell {
	public:
		using FloatingPointType = basic_types::FloatingPointType;
		using PointType = geometry::Point<FloatingPointType>;
		using ShiftType = geometry::Point<int8_t>;
		using BondWithShift = BondWithPoint<ShiftType>;
		using PlaneType = geometry::Plane<FloatingPointType>;

		static constexpr FloatingPointType limit = 2*::std::numeric_limits<FloatingPointType>::epsilon();
		
	private:
		// Data (Owning)
		Container<std::unique_ptr<Vertex>> vertexes;
		Container<std::unique_ptr<Edge>> edges;
		Container<std::unique_ptr<Face>> faces;

		PointType center;
		size_t id;

	public:
		Vertex* add_vertex(const PointType& p) {
			auto temp = std::make_unique<Vertex>(p);
			for (auto& v : vertexes) {
				if (*v == *temp)
					return v.get();
			}
			auto simple_ptr = temp.get();
			vertexes.emplace(std::move(temp));
			return simple_ptr;
		}
		void clipByPlaneAndAddNewFace(const PlaneType& clipping_plane, size_t id_of_another_cell) {

			// Check if the side is correct
			assert(clipping_plane.side(center) > 0);

			// 1. Separate vertexes to sides of plane
			for (auto& v : vertexes)
			{
				if (clipping_plane.side(v->get_point()) < -limit) {
					// Point cutted off
					v->set_state(State::DELETE);
				} 
				else if (clipping_plane.side(v->get_point()) < limit) {
					// Point on the Face
				}
			}

			// 2. Calculate States of edges and faces
			for (auto& e : edges)
			{
				e->calculateState();
			}
			for (auto& f : faces)
			{
				f->calculateState();
				assert(f->get_state() != State::INVALID);
			}

			// 3. Cut edges
			// So we need check all Edges, which ask about MODIFICATION
			for (auto& e : edges)
			{
				if (e->get_state() != State::MODIFICATION) {
					continue;
				}
				auto intersection = e->intersectSegmentPlane(clipping_plane);
				auto new_vertex_ptr = add_vertex(intersection);

				// delete vertex from set
				std::erase_if(e->vertexes, 
							  [](auto* ptr) {
								  if (ptr->get_state() == State::DELETE) {
								  	  return true;
								  }
								  return false;
							  });
				// and add new vertex to set
				e->vertexes.insert(new_vertex_ptr);


				new_vertex_ptr->edges.insert(e.get()); // insert this edge to set of new vertex
				new_vertex_ptr->faces.insert(e->faces.begin(), e->faces.end()); // copy set of faces from the edge
				new_vertex_ptr->set_state(State::MODIFICATION); // Set, that Vertex is on the NEW FACE

				// Now, new vertex complete
			}

			// 4. Cut Faces, which need modification
			// TODO: from here

			// 5. Create new Face

			// 6. Cleanup after modifications


		}
		
	};


	template<class T>
	class VoronoiCell {
	public:
		using PointType = geometry::Point<T>;
		using ShiftType = geometry::Point<int8_t>;
		using BondWithShift = BondWithPoint<ShiftType>;

		struct Face {
			using PolygonType = geometry::Polygon<T>;
			using PlaneType = typename PolygonType::PlaneType;

			// Data
			PolygonType poly;
			BondWithShift bond;

			// Constructors
			constexpr Face() noexcept = default;
			constexpr explicit Face(PolygonType&& p,
									BondWithShift&& b = BondWithShift())
				: poly(std::move(p)), bond(std::move(b)) {}

		};
		using PlaneType = typename Face::PlaneType;
		using FaceVector = ::std::vector<Face>;

	private:
		FaceVector faces_;
		PointType seed_ = PointType(0, 0, 0);

	public:
		// Default constructor creates cube around [0,0,0]
		constexpr VoronoiCell() noexcept {
			initiate_cube_faces_on_seed();
		}
		// Creates cube around seed
		constexpr explicit VoronoiCell(const PointType& seed, bool init = true) noexcept : seed_(seed) {
			if (init) {
				initiate_cube_faces_on_seed();
			}
		}

		/// <summary>
		/// Cut both VoronoiCells by each other
		/// </summary>
		/// <returns> 0 - if correct even if any of the cells are empty, 
		///           1 - if Cells are too close</returns>
		static int interact(VoronoiCell& a, VoronoiCell& b) {
			// Define if faces are not empty
			bool empty_a = a.faces_.empty();
			bool empty_b = b.faces_.empty();

			if (empty_a && empty_b)
				return 0; // Both cells are empty

			PointType d = b.seed_ - a.seed_;

			// Move interval "d" to [-0.5; 0.5]
			for (int i = 0; i < 3; ++i) {
				if (d[i] > 0.5) d[i] -= 1.0;
				else if (d[i] < -0.5) d[i] += 1.0;
			}

			// If points too close - stop
			T d_length = d.r();
			if (d_length < 1e-6)
				return 1;

			PointType normal = d / d_length;
			PointType half_d = d * 0.5;

			// Create plane and clip a by it
			if (!empty_a) {
				PointType midpoint_a = a.seed_ + half_d;
				PlaneType plane(midpoint_a, -normal);
				a.clipByPlaneAndAddNewFace(plane);
			}

			// Using inverted plane for b
			if (!empty_b) {
				PointType midpoint_b = b.seed_ - half_d;
				PlaneType inverted_plane(midpoint_b, normal);
				b.clipByPlaneAndAddNewFace(inverted_plane);
			}

			return 0;
		}
		constexpr const PointType& getSeed() const noexcept {
			return seed_;
		}
		constexpr const FaceVector& getFaces() const noexcept {
			return faces_;
		}

	private:
		void clipByPlaneAndAddNewFace(const PlaneType& clipping_plane) {
			std::vector<PointType> intersection_points;
			std::vector<Face> new_faces;

			// Collect all intersection points and create new faces
			for (auto& face : faces_) {
				// Save original vertices for intersection detection
				std::vector<PointType> original_vertices;
				for (size_t i = 0; i < face.poly.size(); ++i) {
					original_vertices.push_back(face.poly[i]);
				}

				// Clip the face
				face.poly.clipByPlane(clipping_plane);

				// If face remains valid, add it
				if (face.poly.size() >= 3) {
					new_faces.push_back(face);

					// Collect intersection points with this face
					collectIntersectionPoints(original_vertices, clipping_plane, intersection_points);
				}
				// If face becomes invalid (less than 3 vertices), it is not added - thus removed
			}

			// Replace old faces with new ones
			faces_ = std::move(new_faces);

			// Create new face from intersection points if there are enough
			if (intersection_points.size() >= 3) {
				createNewFaceFromIntersections(intersection_points, clipping_plane);
			}
		}

		void collectIntersectionPoints(const std::vector<PointType>& vertices,
									   const PlaneType& plane,
									   std::vector<PointType>& intersection_points) const {
			const size_t n = vertices.size();
			if (n < 3) return;

			PointType prev_vertex = vertices.back();
			T prev_dist = plane.side(prev_vertex);

			for (size_t i = 0; i < n; i++) {
				const PointType& current_vertex = vertices[i];
				const T current_dist = plane.side(current_vertex);

				// If edge intersects the plane, find intersection point
				if (prev_dist * current_dist < 0) {
					const T t = prev_dist / (prev_dist - current_dist);
					const PointType intersection = prev_vertex + (current_vertex - prev_vertex) * t;
					intersection_points.push_back(intersection);
				}

				prev_vertex = current_vertex;
				prev_dist = current_dist;
			}
		}

		void createNewFaceFromIntersections(std::vector<PointType>& points,
											const PlaneType& plane) {
			if (points.size() < 3) return;

			// Order points in correct order (counter-clockwise relative to normal)
			orderPointsOnPlane(points, plane);
			// Remove duplicates
			removeDuplicatePoints(points);

			// Create new face
			typename Face::PolygonType new_polygon(points, plane);
			if (new_polygon.isConvex()) {
				faces_.emplace_back(std::move(new_polygon));
			} else {
				// For test purposes
				// TODO: delete after tests
				return;
			}
		}

		void removeDuplicatePoints(std::vector<PointType>& points) const {
			const auto s = points.size();

			// Remove consecutive duplicates
			for (size_t i = s - 1; i > 0; --i)
			{
				if (PointType::distance(points[i], points[i - 1]) < static_cast<T>(1e-6)) {
					points.erase(points.begin() + i);
				}
			}

			// Check first and last element
			if (PointType::distance(points.front(), points.back()) < static_cast<T>(1e-6)) {
				points.pop_back();
			}
		}

		void orderPointsOnPlane(std::vector<PointType>& points, const PlaneType& plane) {
			if (points.size() < 3) return;

			// Find center of mass of points
			PointType center(0, 0, 0);
			for (const auto& p : points) {
				center += p;
			}
			center = center / static_cast<T>(points.size());

			// Get plane normal
			PointType normal = plane.normal();

			// Choose arbitrary vector in plane (perpendicular to normal)
			PointType reference_vector;
			if (std::abs(normal[0]) > std::abs(normal[1])) {
				reference_vector = PointType(-normal[2], 0, normal[0]);
			} else {
				reference_vector = PointType(0, normal[2], -normal[1]);
			}
			reference_vector = reference_vector / reference_vector.r();

			// Sort points by angle relative to center
			std::ranges::sort(points,
							  [&](const PointType& a, const PointType& b)
							  {
								  PointType vecA = a - center;
								  PointType vecB = b - center;

								  // Project onto plane
								  PointType projA = vecA - normal * PointType::Scalar(vecA, normal);
								  PointType projB = vecB - normal * PointType::Scalar(vecB, normal);

								  if (projA.r() < 1e-10 || projB.r() < 1e-10) {
									  return false; // Points too close to center
								  }

								  // Normalize
								  projA = projA / projA.r();
								  projB = projB / projB.r();

								  // Calculate angles using scalar and vector products
								  T cosA = PointType::Scalar(reference_vector, projA);
								  T sinA = PointType::Scalar(PointType::Vector(reference_vector, projA), normal);
								  T angleA = std::atan2(sinA, cosA);

								  T cosB = PointType::Scalar(reference_vector, projB);
								  T sinB = PointType::Scalar(PointType::Vector(reference_vector, projB), normal);
								  T angleB = std::atan2(sinB, cosB);

								  return angleA < angleB;
							  });
		}


		inline void clipByPlane(const PlaneType& clipping_plane) {
			for (auto& face : faces_) {
				face.poly.clipByPlane(clipping_plane);
			}
		}
		constexpr void initiate_cube_faces_on_seed() {
			std::array<PointType, 8> cube;
			for (int i = 0; i < 8; ++i) {
				cube[i] = base_vertices[i] + seed_;
			}
			faces_.reserve(6);
			for (int i = 0; i < 6; ++i) {
				faces_.emplace_back(Face::PolygonType({cube[face_indices[i][0]],
													  cube[face_indices[i][1]],
													  cube[face_indices[i][2]],
													  cube[face_indices[i][3]]}));
			}
		}

		// Base array of vertices for a cube
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
	};

	template<class T, class AI> class HashedSpace;

	template<class T>
	class VoronoiDiagram {
	public:
		using PointType = geometry::Point<T>;
		using VoronCell = VoronoiCell<T>;

		template <class AI>
		using BondList = ::std::vector<::std::pair<AI, AI>>;
		using PointVector = ::std::vector<PointType>;
		using CellVector = ::std::vector<VoronCell>;
		using BoolVector = ::std::vector<bool>;

		enum class State : unsigned char {
			Uninitialized = 0,
			Cubic_cells = 1,
			Correct_cells = 2
		};
	private:
		//Data
		BoolVector flags_;
		CellVector cells_;
		State state = State::Uninitialized;
	public:
		constexpr VoronoiDiagram() noexcept = default;
		constexpr explicit VoronoiDiagram(const PointVector& points, const BoolVector& flags = BoolVector()) noexcept : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points.size(), true);
			}
			addPoints(points, flags_);
		}

		template<class AI>
		constexpr VoronoiDiagram(const PointVector& points, const BondList<AI>& bonds, const BoolVector& flags = BoolVector()) noexcept : flags_(flags) {
			if (flags.empty()) {
				flags_.resize(points.size(), true);
			}
			addPoints(points, flags_);
			calculateFaces<AI>(bonds);
		}
		void addPoints(const PointVector& points, const BoolVector& flags) noexcept {
			cells_.reserve(points.size());
			for (size_t i = 0; i < points.size(); i++)
			{
				cells_.emplace_back(points[i], flags[i]);
			}
			state = State::Cubic_cells;
		}

		template <class AI>
		int calculateFaces(const BondList<AI>& bonds) noexcept {
			if (state == State::Uninitialized)
				return 1; // Error: VoronoiDiagram not initialized
			for (auto& bond : bonds) {
				auto interaction_result = VoronCell::interact(cells_[bond.first], cells_[bond.second]);
				if (interaction_result != 0)
					return 2; // Error: Cells too close
			}
			state = State::Correct_cells;
			return 0;
		}

		constexpr T calculateLongestDiagonal(const geometry::Matrix<T>& mat) const noexcept {
			T ret = 0;
			for (auto& vcell : cells_) {
				const PointType seed = vcell.getSeed();
				for (const auto& face : vcell.getFaces()) {
					for (size_t i = 0; i < face.poly.size(); i++)
					{
						T val = (mat * (face.poly[i] - seed)).r();
						if (val > ret) ret = val;
					}
				}
			}
			return ret * 2;
		}

		CellVector extractCells() noexcept {
			if (state == State::Uninitialized) {
				return {};
			}
			state = State::Uninitialized;
			return std::move(cells_);
		}
	};

	template<class T>
	class VoronoiFused {
	public:
		using PointType = geometry::Point<T>;

		class Polygon {
		public:
			bool is_inner = false;
			::std::vector<::std::size_t> vert_ids;

			void rotateToCanonical() {
				const size_t n = vert_ids.size();
				if (n <= 1) return;

				// 1. Find position of minimum, O(n)
				size_t min_idx = 0;
				for (size_t i = 1; i < n; ++i) {
					if (vert_ids[i] < vert_ids[min_idx]) {
						min_idx = i;
					}
				}

				// 2. Find right diraction of ring
				const size_t prev_idx = (min_idx == 0)?n - 1:min_idx - 1;
				const size_t next_idx = (min_idx == n - 1)?0:min_idx + 1;
				const bool need_reverse = (vert_ids[prev_idx] < vert_ids[next_idx]);

				// 3. Final rotation on possible reversion
				if (need_reverse) {
					std::reverse(vert_ids.begin(), vert_ids.end());
					// change minimum position after reverse
					const size_t new_min_idx = n - 1 - min_idx;
					if (new_min_idx != 0) {
						std::rotate(vert_ids.begin(), vert_ids.begin() + new_min_idx, vert_ids.end());
					}
				} else {
					if (min_idx != 0) {
						std::rotate(vert_ids.begin(), vert_ids.begin() + min_idx, vert_ids.end());
					}
				}
			}
		};

		using Polyhedra = ::std::vector<size_t>; // Polygon indexes, equal center index

		static constexpr T EPSILON = 0.0001;
	public:
		//Data
		::std::vector<PointType> centers;
		::std::vector<PointType> vertexes;
		::std::vector<Polygon> polygons;
		::std::vector<Polyhedra> polyhedra;
	public:
		void AddCells(const ::std::vector<VoronoiCell<T>>& cells) {

			auto cells_s = cells.size();

			vertexes.clear();
			polygons.clear();
			vertexes.reserve(120 * cells_s);
			polygons.reserve(30 * cells_s);

			polyhedra.clear();
			polyhedra.resize(cells_s);
			centers.clear();
			centers.resize(cells_s);

			// fill vetexes with coppies
			for (size_t i = 0; i < cells_s; i++)
			{
				centers[i] = cells[i].getSeed();
				auto& faces = cells[i].getFaces();
				auto faces_s = faces.size();
				for (size_t j = 0; j < faces_s; j++)
				{
					Polygon p;

					auto& vert = faces[j].poly.getVertixes();
					auto vert_s = vert.size();
					for (size_t k = 0; k < vert_s; k++)
					{
						auto iter = add_to_vertex_union(vertexes, vert[k]);
						p.vert_ids.push_back(iter);
					}
					p.rotateToCanonical();
					auto pgon_it = add_to_polygon_union(polygons, p);
					polyhedra[i].push_back(pgon_it);
				}
			}

		}

	private:
		inline size_t add_to_vertex_union(::std::vector<PointType>& v, const PointType& x) const {
			auto f_It = ::std::ranges::find_if(v,
											   [&x](const PointType& p) {
												   return PointType::distance(x, p) <= EPSILON;
											   });
			if (f_It != v.end())
				return ::std::distance(v.begin(), f_It);
			else {
				v.push_back(x);
				return v.size() - 1;
			}
		}
		inline size_t add_to_polygon_union(::std::vector<Polygon>& v, const Polygon& x) const {
			auto f_It = ::std::ranges::find_if(v,
											   [&x](const Polygon& p) {
												   if (p.is_inner) return false;
												   if (p.vert_ids.size() != x.vert_ids.size())
													   return false;
												   for (size_t i = 0; i < p.vert_ids.size(); ++i) {
													   if (p.vert_ids[i] != x.vert_ids[i])
														   return false;
												   }
												   return true;
											   });
			if (f_It != v.end()) {
				f_It->is_inner = true;
				return ::std::distance(v.begin(), f_It);
			} else {
				v.push_back(x);
				return v.size() - 1;
			}
		}
	};
}
