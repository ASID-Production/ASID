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

// ============================================================================
//  LEVEL 0: Pure Data (POD Structures)
//  Inert, dense, cache-aligned. No methods. Zero polymorphism.
// ============================================================================

namespace cpplib::voronoi {

	using IDtype = uint16_t;
	using LocalIdx = uint16_t;
	using FloatingPointType = basic_types::FloatingPointType;
	using PointType = geometry::Point<FloatingPointType>;
	using MatrixType = geometry::Matrix<FloatingPointType>;

	constexpr IDtype INVALID_ID = std::numeric_limits<IDtype>::max();

	// --- Vertex: 32 bytes ---
	struct alignas(32) Vertex {
		PointType  pos;      // Cartesian position
		FloatingPointType distance; // Distance to atom (service)
	};

	// --- HalfEdge: 8 bytes ---
	struct alignas(8) HalfEdge {
		IDtype origin_vertex_id;
		IDtype polygon_id;
		IDtype next_edge_id;
		IDtype twin_edge_id;
	};

	// --- Polygon: 8 bytes ---
	struct alignas(8) Polygon {
		IDtype first_edge_id;  // Stable global ID of first HalfEdge
		IDtype other_cell;     // ID of another Polyhedron
		ShiftCode other_shift; // Shift of another Polyhedron
	};

	// --- Polyhedron: 32 bytes ---
	struct alignas(32) Polyhedron {
		PointType  center;         // Polyhedron center (atom position)
		FloatingPointType volume;  // Computed volume
	};
	static_assert(sizeof(Polyhedron) <= 64, "Polyhedron should fit in one cache-line");

} // namespace cpplib::voronoi

// ============================================================================
//  LEVEL 1: Technical Storage (SoA Container)
//  Memory management only: add, remove (swap-and-pop O(1)), point access.
//  Knows nothing about geometry or physics.
// ============================================================================

namespace cpplib::voronoi {

	template <typename T, size_t MAX_POOL = 256>
	class Storage {
	public:
		explicit Storage(IDtype max_id = MAX_POOL - 1)
			: mask_(static_cast<size_t>(max_id), INVALID_ID) {
		}

		// Insert object: returns local dense index
		LocalIdx add(const T& obj) {
			assert(next_ < mask_.size());
			assert(mask_[next_] == INVALID_ID);
			LocalIdx idx = static_cast<LocalIdx>(data_.size());
			data_.push_back(obj);
			ids_.push_back(next_);
			mask_[next_] = idx;
			++next_;
			return idx;
		}

		// O(1) removal: swap-and-pop + single mask cell update for the moved element
		void remove_by_id(IDtype global_id) noexcept {
			if (global_id >= mask_.size()) return;
			LocalIdx idx = mask_[global_id];
			if (idx == INVALID_ID) return;

			LocalIdx last = static_cast<LocalIdx>(data_.size() - 1);
			if (idx != last) {
				data_[idx] = std::move(data_[last]);
				ids_[idx] = ids_[last];
				mask_[ids_[idx]] = idx;
			}
			data_.pop_back();
			ids_.pop_back();
			mask_[global_id] = INVALID_ID;
		}

		// Branchless point access by global_id: 1 L1 cycle
		const T& get_by_id(IDtype global_id) const noexcept {
			assert(global_id < mask_.size() && mask_[global_id] != INVALID_ID);

			return data_[mask_[global_id]];
		}
		T& get_by_id(IDtype global_id) noexcept {
			assert(global_id < mask_.size() && mask_[global_id] != INVALID_ID);

			return data_[mask_[global_id]];
		}

		// Dense range for iteration
		const T* data() const noexcept {
			return data_.data();
		}
		T* data() noexcept {
			return data_.data();
		}
		size_t size() const noexcept {
			return data_.size();
		}
		bool empty() const noexcept {
			return data_.empty();
		}

		// Iteration over all active elements
		template <typename Fn>
		void for_each(Fn&& fn) {
			for (size_t i = 0; i < data_.size(); ++i)
				fn(data_[i]);
		}

		IDtype get_id(LocalIdx index) const noexcept {
			return ids_[index];
		}

		void clear() noexcept {
			data_.clear();
			std::fill(mask_.begin(), mask_.end(), INVALID_ID);
		}

		void reserve(size_t n) {
			data_.reserve(n);
		}

	private:
		std::vector<T>        data_;
		std::vector<IDtype>   ids_;
		std::vector<LocalIdx> mask_; // global_id -> local dense index (L1-friendly)
		IDtype next_ = 0; // next empty ID
	};


// --- United Storage ---
	struct PolyhedronData {
		Storage<Vertex>   vert;
		Storage<HalfEdge> edge;
		Storage<Polygon>  face;

		void clear() noexcept {
			vert.clear(); edge.clear(); face.clear();
		}
		void reserve(size_t nv, size_t ne, size_t nf) {
			vert.reserve(nv); edge.reserve(ne); face.reserve(nf);
		}

		IDtype add_vertex(const PointType& p) {
			Vertex v{}; v.pos = p; v.distance = FloatingPointType(0);
			const LocalIdx idx = vert.add(v);
			return static_cast<IDtype>(idx);
		}
		IDtype add_halfedge(IDtype origin, IDtype poly, IDtype next, IDtype twin) {
			HalfEdge h{}; h.origin_vertex_id = origin;
			h.polygon_id = poly; h.next_edge_id = next; h.twin_edge_id = twin;
			const LocalIdx idx = edge.add(h);
			return static_cast<IDtype>(idx);
		}
		IDtype add_face(IDtype first_he, IDtype other_cell, ShiftCode shift) {
			Polygon p{}; p.first_edge_id = first_he;
			p.other_cell = other_cell; p.other_shift = shift;
			const LocalIdx idx = face.add(p);
			return static_cast<IDtype>(idx);
		}
	};
} // namespace cpplib::voronoi

// ============================================================================
//  LEVEL 2: Math Engine (Low-Level, Stateless)
//  Pure static inline functions. Zero fields. Ideal for -O3 / AVX2 / FMA.
//  All functions accept data by const-reference from outside.
// ============================================================================

namespace cpplib::voronoi {

	struct MathEngine {
		using PlaneType = geometry::Plane<FloatingPointType>;

		static constexpr FloatingPointType EPS = static_cast<FloatingPointType>(1.0e-10);

		// --- Triangle area from 3 points (cross-product magnitude / 2) ---
		static inline FloatingPointType triangle_area(const PointType& vector_AB,
													  const PointType& vector_AC) noexcept {
			return PointType::Vector(vector_AB, vector_AC).r()* FloatingPointType(0.5);
		}

		// --- Polygon normal via Newell's method (robust for concave polygons) ---
		static inline PointType polygon_normal(const PolyhedronData& d,
											   IDtype first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			PointType v0 = d.vert.get_by_id(e0.origin_vertex_id).pos;

			IDtype next_he = e0.next_edge_id;
			PointType v1 = d.vert.get_by_id(d.edge.get_by_id(next_he).origin_vertex_id).pos;

			FloatingPointType nx = 0, ny = 0, nz = 0;

			while (true) {
				nx += (v0[1] - v1[1]) * (v0[2] + v1[2]);
				ny += (v0[2] - v1[2]) * (v0[0] + v1[0]);
				nz += (v0[0] - v1[0]) * (v0[1] + v1[1]);

				if (next_he == first_he) break;

				v0 = v1;
				const HalfEdge& e = d.edge.get_by_id(next_he);
				next_he = e.next_edge_id;
				v1 = d.vert.get_by_id(d.edge.get_by_id(next_he).origin_vertex_id).pos;
			}

			const FloatingPointType len = std::sqrt(nx * nx + ny * ny + nz * nz);
			if (len > EPS) {
				nx /= len; ny /= len; nz /= len;
			}
			return PointType(nx, ny, nz);
		}

		// --- Polygon area (fan triangulation) ---
		static inline FloatingPointType polygon_area(const PolyhedronData& d,
													 IDtype first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const PointType p0 = d.vert.get_by_id(e0.origin_vertex_id).pos;

			const HalfEdge& e1 = d.edge.get_by_id(e0.next_edge_id);

			PointType v1 = d.vert.get_by_id(e1.origin_vertex_id).pos - p0;

			IDtype e2 = e1.next_edge_id;
			FloatingPointType area = FloatingPointType(0);

			while (e2 != first_he) {
				const HalfEdge& e2_he = d.edge.get_by_id(e2);
				const PointType p2 = d.vert.get_by_id(e2_he.origin_vertex_id).pos;
				PointType v2 = p2 - p0;
				area += triangle_area(v1, v2);
				v1 = v2;
				e2 = e2_he.next_edge_id;
			}
			return area;
		}

		static inline PlaneType polygon_plane(const PolyhedronData& d,
															 IDtype first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const IDtype e1_id = e0.next_edge_id;
			const HalfEdge& e1 = d.edge.get_by_id(e1_id);
			const IDtype e2_id = e1.next_edge_id;
			return PlaneType(
				d.vert.get_by_id(e0.origin_vertex_id).pos,
				d.vert.get_by_id(e1.origin_vertex_id).pos,
				d.vert.get_by_id(d.edge.get_by_id(e2_id).origin_vertex_id).pos);
		}

		// --- Solid angle of a polygon face at point 'origin' (van Oosterom-Strackee) ---
		static inline FloatingPointType polygon_solid_angle(const PolyhedronData& d,
															IDtype first_he,
															const PointType& origin) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const PointType v0 = d.vert.get_by_id(e0.origin_vertex_id).pos - origin;
			const FloatingPointType r0 = v0.r();

			IDtype next_he = e0.next_edge_id;
			const PointType v1 = d.vert.get_by_id(
				d.edge.get_by_id(next_he).origin_vertex_id).pos - origin;

			PointType         v_prev = v1;
			FloatingPointType r_prev = v1.r();
			FloatingPointType dot0_prev = PointType::Scalar(v0, v1);

			next_he = d.edge.get_by_id(next_he).next_edge_id;

			FloatingPointType angle = FloatingPointType(0);

			while (next_he != first_he) {
				const HalfEdge& he = d.edge.get_by_id(next_he);
				const PointType v_curr = d.vert.get_by_id(he.origin_vertex_id).pos - origin;

				const FloatingPointType r_curr = v_curr.r();
				const FloatingPointType dot0_curr = PointType::Scalar(v0, v_curr);

				const FloatingPointType cross = std::abs(
					PointType::Scalar(v0, PointType::Vector(v_prev, v_curr)));

				const FloatingPointType denom =
					r0 * r_prev * r_curr
					+ dot0_prev * r_curr
					+ dot0_curr * r_prev
					+ PointType::Scalar(v_prev, v_curr) * r0;

				if (cross > EPS) angle += std::atan2(cross, denom);

				v_prev = v_curr;
				r_prev = r_curr;
				dot0_prev = dot0_curr;

				next_he = he.next_edge_id;
			}

			return angle * FloatingPointType(2);
		}

		// --- Volume: sum of pyramid contributions per face ---
		//   volume = sum over faces: dist(center, face_plane) * face_area / 3
		static inline FloatingPointType cell_volume(const PolyhedronData& d,
													const PointType& center) noexcept {
			const Polygon* polys = d.face.data();
			const size_t   nf = d.face.size();

			FloatingPointType volume = FloatingPointType(0);
			for (size_t i = 0; i < nf; ++i) {
				const IDtype fh = polys[i].first_edge_id;
				const FloatingPointType area = polygon_area(d, fh);
				const auto              pl = polygon_plane(d, fh);
				volume += area * std::abs(pl.distance(center));
			}
			return volume * FloatingPointType(1.0 / 3.0);
		}

		// --- Edge length ---
		static inline FloatingPointType edge_length(const PointType& a,
													const PointType& b) noexcept {
			return PointType::distance(a, b);
		}
	};

} // namespace cpplib::voronoi

// ============================================================================
//  LEVEL 3: Topological Context (Mid-Level Manager)
//  OWNS all Level-1 storages. Knows PBC, supercell, spatial hashing.
//  Coordinates topology updates and defragmentation.
// ============================================================================

namespace cpplib::voronoi {

	struct VoronoiCell {
		Polyhedron     meta;   // {center, volume}
		PolyhedronData topo;   // {vert, edge, face}

		void clear() noexcept {
			meta = Polyhedron{};
			topo.clear();
		}
		void reserve(size_t nv, size_t ne, size_t nf) {
			topo.reserve(nv, ne, nf);
		}
	};

	struct VoronoiContext {
		std::vector<VoronoiCell> cells;

		MatrixType FtoC, CtoF;
		int sup_x = 1, sup_y = 1, sup_z = 1;
		uint32_t          max_neighbors = 12;
		FloatingPointType cutoff = FloatingPointType(0);

		struct SpatialHash {
			FloatingPointType cell_size = FloatingPointType(1);
			std::unordered_map<uint64_t, std::vector<LocalIdx>> grid;

			static inline uint64_t key(int64_t ix, int64_t iy, int64_t iz) noexcept {
				return (uint64_t(ix) * 73856093ull)
					^ (uint64_t(iy) * 19349663ull)
					^ (uint64_t(iz) * 83492791ull);
			}

			void rebuild(const PointType* pts, size_t count, FloatingPointType cs) {
				cell_size = cs;
				grid.clear();
				grid.reserve(count * 2);
				for (size_t i = 0; i < count; ++i) {
					const int64_t ix = static_cast<int64_t>(std::floor(pts[i][0] / cs));
					const int64_t iy = static_cast<int64_t>(std::floor(pts[i][1] / cs));
					const int64_t iz = static_cast<int64_t>(std::floor(pts[i][2] / cs));
					grid[key(ix, iy, iz)].push_back(static_cast<LocalIdx>(i));
				}
			}

			template <typename Fn>
			void query(const PointType& c, FloatingPointType radius, Fn&& fn) const {
				const int32_t r = static_cast<int32_t>(std::ceil(radius / cell_size));
				const int64_t cx = static_cast<int64_t>(std::floor(c[0] / cell_size));
				const int64_t cy = static_cast<int64_t>(std::floor(c[1] / cell_size));
				const int64_t cz = static_cast<int64_t>(std::floor(c[2] / cell_size));
				for (int32_t dx = -r; dx <= r; ++dx)
					for (int32_t dy = -r; dy <= r; ++dy)
						for (int32_t dz = -r; dz <= r; ++dz) {
							auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
							if (it != grid.end())
								for (auto idx : it->second) fn(idx);
						}
			}
		} hash;

		// ------------------------------------------------------------------
		void clear() noexcept {
			cells.clear();
			hash.grid.clear();
		}

		VoronoiCell& cell(size_t i)       noexcept {
			return cells[i];
		}
		const VoronoiCell& cell(size_t i) const noexcept {
			return cells[i];
		}

		static void recompute_cell_metrics(VoronoiCell& c) noexcept {
			c.meta.volume = MathEngine::cell_volume(c.topo, c.meta.center);
		}

		struct Stats {
			size_t n_cells = 0;
			size_t n_vertices = 0;
			size_t n_edges = 0;
			size_t n_faces = 0;
			FloatingPointType total_volume = FloatingPointType(0);
		};

		Stats get_stats() const noexcept {
			Stats s;
			s.n_cells = cells.size();
			for (const auto& c : cells) {
				s.n_vertices += c.topo.vert.size();
				s.n_edges += c.topo.edge.size();
				s.n_faces += c.topo.face.size();
				s.total_volume += c.meta.volume;
			}
			return s;
		}
	};

} // namespace cpplib::voronoi


// ============================================================================
//  LEVEL 4: Orchestrator / Top-Level Pipeline
//  "Pure reason": drives iteration loops, coordinates L3 <-> L2.
//  Declarative: prepare data -> compute -> update state.
// ============================================================================

namespace cpplib::voronoi {

	struct VoronoiPipeline {
		struct Config {
			uint32_t          max_neighbors = 12;
			FloatingPointType cutoff_scale = 1.0;
			bool              use_supercell = true;
			int               sup_x = 1, sup_y = 1, sup_z = 1;
			bool              compute_volumes = true;
			bool              compute_areas = true;
			bool              compute_solid_angles = true;
		};

		Config           config;
		VoronoiContext   context;

		// --- Main entry: build Voronoi diagram for a set of atoms ---
		void build(const PointType* atom_pos, size_t n_atoms,
				   const MatrixType& FtoC) {
			context.FtoC = FtoC;
			context.max_neighbors = config.max_neighbors;
			context.vertex_storage.reserve(n_atoms * 4);
			context.cell_storage.reserve(n_atoms);

			// Step 1: Spatial hashing + neighbor list (Level 3)
			FloatingPointType avg_dist = estimate_avg_distance(atom_pos, n_atoms);
			FloatingPointType cutoff = avg_dist * config.cutoff_scale;
			context.build_neighbor_list(atom_pos, n_atoms, cutoff);

			// Step 2: Build one cell per atom (Level 3 orchestrates, Level 2 computes)
			for (size_t i = 0; i < n_atoms; ++i)
				build_single_cell(static_cast<IDtype>(i), atom_pos, n_atoms, FtoC);

			// Step 3: Extract edges and vertices from face topology
			build_edges_and_vertices();

			// Step 4: Final metric recomputation (Level 2)
			if (config.compute_volumes)      recompute_volumes();
			if (config.compute_areas)        recompute_areas();
			if (config.compute_solid_angles) recompute_solid_angles();
		}

		// --- Incremental update after MD step (atom positions changed) ---
		void update(const PointType* new_pos, size_t n_atoms) {
			for (size_t i = 0; i < n_atoms; ++i) {
				IDtype cid = static_cast<IDtype>(i);
				if (Cell* cell = context.cell_storage.get_by_id(cid)) {
					cell->center = new_pos[i];
					context.update_cell_topology(cid);
				}
			}
			if (config.compute_volumes)      recompute_volumes();
			if (config.compute_areas)        recompute_areas();
			if (config.compute_solid_angles) recompute_solid_angles();
		}

		// --- Public accessors for external consumers ---
		const Cell* get_cell(IDtype id)   const noexcept {
			return context.cell_storage.get_by_id(id);
		}
		const Face* get_face(IDtype id)   const noexcept {
			return context.face_storage.get_by_id(id);
		}
		const Vertex* get_vertex(IDtype id) const noexcept {
			return context.vertex_storage.get_by_id(id);
		}

		size_t              n_cells() const noexcept {
			return context.cell_storage.size();
		}
		VoronoiContext::Stats stats() const noexcept {
			return context.get_stats();
		}

	private:
		// --- Estimate average nearest-neighbor distance (for cutoff heuristic) ---
		static inline FloatingPointType estimate_avg_distance(
			const PointType* pts, size_t n) noexcept {
			if (n < 2) return FloatingPointType(1.0);
			FloatingPointType sum = 0;
			size_t sample = std::min(n, size_t{100});
			for (size_t i = 0; i < sample; ++i) {
				FloatingPointType min_d = std::numeric_limits<FloatingPointType>::max();
				for (size_t j = 0; j < n; ++j) {
					if (i == j) continue;
					min_d = std::min(min_d, MathEngine::distance(pts[i], pts[j]));
				}
				sum += min_d;
			}
			return sum / static_cast<FloatingPointType>(sample);
		}

		// --- Build Voronoi cell for atom i ---
		void build_single_cell(IDtype cell_id, const PointType* atoms,
							   size_t n_atoms, const MatrixType& FtoC) {
			(void)FtoC;
			(void)n_atoms;

			Cell cell{};
			cell.id = cell_id;
			cell.atom_id = cell_id;
			cell.center = atoms[cell_id];
			cell.volume = 0;
			cell.n_faces = 0;
			cell.n_vertices = 0;
			cell.n_edges = 0;
			context.cell_storage.add(cell_id, cell);

			// Register atom center as a vertex in storage
			Vertex vcenter{};
			vcenter.pos = atoms[cell_id];
			vcenter.id = cell_id * 1024;  // local vertex ID namespace
			vcenter.cell_id = cell_id;
			context.vertex_storage.add(vcenter.id, vcenter);

			// Phase A: Collect neighbor atoms (Level 3 spatial hash)
			// Phase B: Construct bisecting planes (Level 2)
			// Phase C: Intersect planes -> face polygons (Level 2)
			// Phase D: Write Face records into context (Level 3)
			// [Full geometric construction elided — depends on specific
			//  Voronoi algorithm variant (Fortune, incremental, etc.)]
		}

		// --- Derive edges and unique vertices from face topology ---
		void build_edges_and_vertices() {
			context.face_storage.for_each([&](const Face& f) {
				const TopologyBuffer& fv = context.face_to_vertices[f.id];
				for (LocalIdx k = 0; k + 1 < fv.count; ++k) {
					IDtype v0 = fv.data[k];
					IDtype v1 = fv.data[k + 1];

					Vertex* p0 = context.vertex_storage.get_by_id(v0);
					Vertex* p1 = context.vertex_storage.get_by_id(v1);
					if (!p0 || !p1) continue;

					Edge e{};
					e.v0 = v0;
					e.v1 = v1;
					e.face_id = f.id;
					e.length = MathEngine::edge_length(p0->pos, p1->pos);

					IDtype e_id = v0 ^ (v1 << 1);
					if (!context.edge_storage.get_by_id(e_id))
						context.edge_storage.add(e_id, e);
				}
										  });
		}

		// --- Recompute all cell volumes (delegates to Level 2) ---
		void recompute_volumes() {
			context.cell_storage.for_each([&](const Cell& c) {
				const Face* faces = context.face_storage.data();
				FloatingPointType vol = MathEngine::cell_volume(
					c.center,
					faces,
					c.n_faces,
					context.vertex_storage.data(),
					context.face_to_vertices.data()
				);
				if (Cell* cell = context.cell_storage.get_by_id(c.id))
					cell->volume = vol;
										  });
		}

		// --- Recompute all face areas ---
		void recompute_areas() {
			context.face_storage.for_each([&](const Face& f) {
				const TopologyBuffer& fv = context.face_to_vertices[f.id];
				if (fv.count < 3) return;
				const Vertex* vdata = context.vertex_storage.data();
				std::array<PointType, 16> cart_buf;
				uint32_t n = std::min<uint32_t>(fv.count, 16);
				for (uint32_t i = 0; i < n; ++i)
					cart_buf[i] = vdata[fv.data[i]].pos;

				FloatingPointType area = MathEngine::polygon_area(cart_buf.data(), n);
				if (Face* face = context.face_storage.get_by_id(f.id))
					face->area = area;
										  });
		}

		// --- Recompute solid angles at cell centers ---
		void recompute_solid_angles() {
			context.face_storage.for_each([&](const Face& f) {
				const TopologyBuffer& fv = context.face_to_vertices[f.id];
				if (fv.count < 3) return;
				Cell* cell = context.cell_storage.get_by_id(f.cell_id);
				if (!cell) return;

				const Vertex* vdata = context.vertex_storage.data();
				std::array<PointType, 16> cart_buf;
				uint32_t n = std::min<uint32_t>(fv.count, 16);
				for (uint32_t i = 0; i < n; ++i)
					cart_buf[i] = vdata[fv.data[i]].pos;

				FloatingPointType sa = MathEngine::solid_angle(
					cart_buf.data(), n, cell->center);
				if (Face* face = context.face_storage.get_by_id(f.id))
					face->solid_angle = sa;
										  });
		}
	};

}

namespace cpplib::voronoi_old {
	// Forward declarations
	struct Vertex;
	struct Edge;
	struct Face;
	struct Cell;

	/// @brief Container type for non-owning pointers in Voronoi structures
	/// @tparam T Pointer type to store
	template<class T>
	using Container = ::std::unordered_set<T>;

	/// @brief State of a Voronoi geometric object during construction
	///
	/// Objects transition through states during plane clipping operations:
	/// 
	/// - DELETE: Object is outside the valid region
	/// 
	/// - VALID: Object is fully valid and unchanged
	/// 
	/// - INVALID: Object is in an inconsistent state (error condition)
	/// 
	/// - MODIFICATION: Object is being modified by a clipping operation
	/// 
	/// - ONPLANE: Object is laying directely on a clipping plane
	enum class State : char {
		DELETE = 0,       ///< Object should be deleted
		VALID = 1,        ///< Object is valid
		INVALID = 2,      ///< Object is in invalid state
		MODIFICATION = 3, ///< Object is being modified
		ONPLANE = 4       ///< Object is on a clipping plane
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
			: state(s), id(ID) {
		}
	};

	/// @brief Epsilon for comparing vertex positions
	///
	/// Vertices within this distance are considered equal to handle
	/// floating-point precision issues.
	constexpr basic_types::FloatingPointType EPSILON = 1E-6;

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
		Vertex() noexcept : Object(0, State::VALID) {
		}

/// @brief Construct vertex with given ID
/// @param ID Unique identifier for this vertex
		explicit Vertex(uint32_t ID) noexcept : Object(ID, State::VALID) {
		}

/// @brief Construct vertex with ID and position
/// @param ID Unique identifier
/// @param p Position point
		Vertex(uint32_t ID, const PointType& p) noexcept : Object(ID, State::VALID), point(p) {
		}

/// @brief Construct vertex with ID and position (move)
/// @param ID Unique identifier
/// @param p Position point (moved)
		Vertex(uint32_t ID, PointType&& p) noexcept : Object(ID, State::VALID), point(std::move(p)) {
		}

/// @brief Compare vertices for equality using spatial epsilon
/// @param a First vertex
/// @param b Second vertex
/// @return True if vertices are spatially equivalent
///
/// Uses EPSILON to handle floating-point precision.
		inline friend bool operator==(const Vertex& a, const Vertex& b) {
			static constexpr basic_types::FloatingPointType EPSILONSQ = EPSILON * EPSILON;
			return geometry::Point<basic_types::FloatingPointType>::distanceSq(a.point, b.point) < EPSILONSQ;
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
		/// An edge is MODIFICATION if it's vertices are separated by clipping plane.
		/// An edge is ONPLANE if both it's vertices are laying on clipping plane.
		/// An edge is DELETE if both vertices are deleted.
		/// Otherwise, it's INVALID.
		inline State calculateState() noexcept {
			using enum State;
			if (vertices.size() != 2 || faces.size() != 2) [[unlikely]] {
				set_state(INVALID);
				return get_state();
			}
			auto it = vertices.begin();
			State s1 = (*it)->get_state();
			++it;
			State s2 = (*it)->get_state();

			if (s1 == DELETE && s2 == DELETE) {
				set_state(DELETE);
			} else if ((s1 == VALID && s2 == DELETE) || (s1 == DELETE && s2 == VALID)) {
				set_state(MODIFICATION); // Intersection
			} else if ((s1 == VALID && s2 == ONPLANE) || (s1 == ONPLANE && s2 == VALID)) {
				set_state(VALID);
			} else if (s1 == ONPLANE && s2 == ONPLANE) {
				set_state(ONPLANE);
			} else if ((s1 == DELETE && s2 == ONPLANE) || (s1 == ONPLANE && s2 == DELETE)) {
				set_state(DELETE);
			} else {
				set_state(VALID);
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
			if (std::abs(denom) < EPSILON)
				assert(std::abs(denom) >= EPSILON);

			auto t = -(plane.a[0] * v1->point[0] +
					   plane.a[1] * v1->point[1] +
					   plane.a[2] * v1->point[2] +
					   plane.a[3]) / denom;

			if (t <= 0.0 || t >= 1.0)
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
		explicit Face(uint32_t ID) : Object(ID) {
		}

/// @brief Calculate and update the state of this face
/// @return Updated State
///
/// A face is VALID if it has equal numbers of vertices and edges, and all edges are valid.
/// A face is MODIFICATION if any edge is being modified.
/// A face is DELETE if all edges are deleted.
/// A face with exactly one valid edge is INVALID (error condition).
		inline State calculateState() {
			using enum State;

			uint32_t v_size = vertices.size();
			uint32_t v_valid = 0;
			uint32_t v_pln = 0;
			uint32_t v_inv = 0;
			for (auto& ver : vertices)
			{
				switch (ver->get_state()) {
				case VALID:
					v_valid++;
					break;
				case ONPLANE:
					v_pln++;
					break;
				}
			}

			if (v_valid == 0) {
				set_state(DELETE);
				for (auto& e : edges) {
					if (e->get_state() != DELETE) {
						e->faces.erase(this);
					}
				}
				return get_state();
			} else if (v_valid + v_pln < v_size) {
				set_state(MODIFICATION);
				return get_state();
			}

			set_state(VALID);
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
		static constexpr FloatingPointType limit = EPSILON;

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
					v->set_state(ONPLANE);
				}
			}

			// Early exit if no vertices were cut
			if (modified == false) {
				for (auto& v : vertices)
				{
					if (v->get_state() == ONPLANE) {
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
				if (f->get_state() == INVALID) {
					assert(f->get_state() != INVALID);
				}
			}

			// 3. Cut edges
			// Check all Edges which need MODIFICATION
			for (auto& e : edges)
			{
				if (e->get_state() != MODIFICATION) {
					continue;
				}

				auto v1 = *e->vertices.begin();
				auto v2 = e->get_second_vertex(v1);
				if (!((v1->get_state() == State::VALID && v2->get_state() == State::DELETE) ||
					  (v1->get_state() == State::DELETE && v2->get_state() == State::VALID))) {
					e->set_state(State::VALID);
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
				new_vertex_ptr->set_state(ONPLANE); // Mark as on the NEW FACE

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
					if (v1->get_state() == ONPLANE) {
						iter_vertex++;
						break;
					}
					iter_vertex++;
				}
				while (iter_vertex != f->vertices.cend()) {
					v2 = *iter_vertex;
					if (v2->get_state() == ONPLANE) {
						break;
					}
					iter_vertex++;
				}
				if (iter_vertex == f->vertices.cend())
					assert(iter_vertex != f->vertices.cend());

				// Create edge connecting the two new vertices
				const auto& new_edge = edges.emplace_back(std::make_unique<Edge>(static_cast<uint32_t>(edges.size()), v1, v2));
				new_edge->faces.emplace(f.get());
				f->edges.emplace(new_edge.get());
				new_edge->set_state(ONPLANE);
				f->set_state(VALID);
			}

			// 5. Create new Face (the clipping plane becomes a face)
			const auto& new_face = faces.emplace_back(std::make_unique<Face>(static_cast<uint32_t>(faces.size())));
			auto raw_face_ptr = new_face.get();

			// Add all ONPLANE edges to the new face
			for (const auto& e : edges)
			{
				if (e->get_state() == ONPLANE) {
					raw_face_ptr->edges.emplace(e.get());
					e->faces.emplace(raw_face_ptr);
					e->set_state(VALID);
				}
			}

			// Add all ONPLANE vertices to the new face
			for (const auto& v : vertices)
			{
				if (v->get_state() == ONPLANE) {
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
				std::cout << i << std::endl;
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
			geometry::ShiftCode::ShiftPoint second_shift; ///< Shift point of the second atom
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
				} else if (v1 < v2) {
					i1++;
				} else {
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

			FloatingPointType area = PointType::Vector(cart_verts.front(), cart_verts.back()).r();
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
							   cart_verts_r[0] * cart_verts_r[i - 1] * cart_verts_r[i] +
							   PointType::Scalar(cart_verts[0], cart_verts[i - 1]) * cart_verts_r[i] +
							   PointType::Scalar(cart_verts[0], cart_verts[i]) * cart_verts_r[i - 1] +
							   PointType::Scalar(cart_verts[i - 1], cart_verts[i]) * cart_verts_r[0]);
			}

			return angle * 2;
		}

		// Requires area to be calculated in polygons
		FloatingPointType calculate_volume(const Polyhedron& p, const geometry::Matrix<FloatingPointType>& FtoC) const {
			FloatingPointType volume = 0;

			for (auto poly : p.poly_ids) {
				auto i1 = polygons[poly].vert_ids[0];
				auto i2 = polygons[poly].vert_ids[1];
				auto i3 = polygons[poly].vert_ids[2];
				geometry::Plane plane(FtoC * vertices[i1], FtoC * vertices[i2], FtoC * vertices[i3]);
				volume += plane.distance(FtoC * p.center) * polygons[poly].area * FloatingPointType(1. / 3);
			}

			return volume;
		}
	};
}
