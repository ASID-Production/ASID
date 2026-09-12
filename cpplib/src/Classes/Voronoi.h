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
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
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

	// ------------------------------------------------------------------------
	//  ID types
	//
	//  LocalId  — ID of an element inside one Voronoi cell. Bounded by the
	//             per-cell element count (V ≈ 260, E ≈ 390, F ≈ 130 for a
	//             dense 100-neighbour shell), so uint16_t is ample.
	//
	//  LocalIdx — dense index inside a Storage<>. Same size as LocalId but
	//             semantically distinct: it never contains gaps.
	//
	//  AtomId   — global atom index in VoronoiContext::cells. Grows with the
	//             system size, up to 10^7 atoms, so uint32_t.
	// ------------------------------------------------------------------------
	using LocalId = uint16_t;
	using LocalIdx = uint16_t;
	using AtomId = uint32_t;


	using ExternalFloatType = basic_types::FloatingPointType;
	using ExternalPointType = geometry::Point<ExternalFloatType>;
	using ExternalUnitCell = geometry::Cell<ExternalFloatType>;

    // Internal working precision.
	using InternalFloatType = double;

	using PointType = geometry::Point<InternalFloatType>;
	using MatrixType = geometry::Matrix<InternalFloatType>;
	using ShiftCode = geometry::ShiftCode;
	using UnitCell = geometry::Cell<InternalFloatType>;

	constexpr LocalId INVALID_LOCAL_ID = std::numeric_limits<LocalId>::max();
	constexpr AtomId  INVALID_ATOM_ID = std::numeric_limits<AtomId>::max();

	// --- Vertex: 32 bytes ---
	//
	//  Layout (packed to 32 bytes so two vertices fit one 64-byte cache
	//  line without straddling):
	//
	//    pos       24 B — Cartesian position
	//    dist_sq    4 B — cached squared distance to the cell center. Used
	//                     by max_vertex_dist_sq, which runs several times
	//                     per cell during clipping; caching it turns a
	//                     10-op scan over 24-byte Points into a 1-op scan
	//                     over 4-byte floats. float precision (~1e-5 A^2
	//                     at 10 A) is far beyond what the early-exit
	//                     bound needs.
	//    edge_hint  2 B — LocalId of one incident half-edge. Reserved for
	//                     a future O(1) vertex -> incident-edge lookup;
	//                     not used yet.
	//    _flags     2 B — bit-packed state. Reserved.
	struct alignas(32) Vertex {
		PointType pos;
		float     dist_sq = 0.0f;
		LocalId   edge_hint = INVALID_LOCAL_ID;
		uint16_t  _flags = 0;
	};
	static_assert(sizeof(Vertex) == 32, "Vertex must stay 32 bytes");

	// --- HalfEdge: 8 bytes ---
	//  All four fields are LocalId: they refer to elements of the same cell.
	struct alignas(8) HalfEdge {
		LocalId origin_vertex_id;
		LocalId polygon_id;
		LocalId next_edge_id;
		LocalId twin_edge_id;
	};
	static_assert(sizeof(HalfEdge) == 8, "HalfEdge must fit in 8 bytes");

	// --- Polygon: 8 bytes ---
	//  Field order chosen so the struct packs into exactly 8 bytes with no
	//  internal padding: AtomId(4) + LocalId(2) + ShiftCode(1) + pad(1).
	struct alignas(8) Polygon {
		AtomId    other_cell;     // Global atom ID of the neighbouring cell
		LocalId   first_edge_id;  // Local half-edge ID, entry point of the cycle
		ShiftCode other_shift;    // Shift of the neighbour
		uint8_t   _pad = 0;       // Explicit padding
	};
	static_assert(sizeof(Polygon) == 8, "Polygon must fit in 8 bytes");

	// --- Polyhedron: 32 bytes ---
	struct alignas(32) Polyhedron {
		PointType  center;         // Polyhedron center (atom position)
		InternalFloatType volume;  // Computed volume
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
		explicit Storage(LocalId max_id = MAX_POOL - 1)
			: mask_(static_cast<size_t>(max_id), INVALID_LOCAL_ID) {
		}

		// Insert object: returns local dense index
		LocalIdx add(const T& obj) {
			assert(next_ < mask_.size());
			assert(mask_[next_] == INVALID_LOCAL_ID);
			LocalIdx idx = static_cast<LocalIdx>(data_.size());
			data_.push_back(obj);
			ids_.push_back(next_);
			mask_[next_] = idx;
			++next_;
			return idx;
		}

		// O(1) removal: swap-and-pop + single mask cell update for the moved element
		void remove_by_id(LocalId global_id) noexcept {
			if (global_id >= mask_.size()) return;
			LocalIdx idx = mask_[global_id];
			if (idx == INVALID_LOCAL_ID) return;

			LocalIdx last = static_cast<LocalIdx>(data_.size() - 1);
			if (idx != last) {
				data_[idx] = std::move(data_[last]);
				ids_[idx] = ids_[last];
				mask_[ids_[idx]] = idx;
			}
			data_.pop_back();
			ids_.pop_back();
			mask_[global_id] = INVALID_LOCAL_ID;
		}

		// Branchless point access by global_id: 1 L1 cycle
		const T& get_by_id(LocalId global_id) const noexcept {
			assert(global_id < mask_.size() && mask_[global_id] != INVALID_LOCAL_ID);
			return data_[mask_[global_id]];
		}
		T& get_by_id(LocalId global_id) noexcept {
			assert(global_id < mask_.size() && mask_[global_id] != INVALID_LOCAL_ID);
			return data_[mask_[global_id]];
		}
		LocalIdx index_of(LocalId global_id) const noexcept {
			assert(global_id < mask_.size() && mask_[global_id] != INVALID_LOCAL_ID);
			return mask_[global_id];
		}

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

		template <typename Fn>
		void for_each(Fn&& fn) {
			for (size_t i = 0; i < data_.size(); ++i)
				fn(data_[i]);
		}

		LocalId get_id(LocalIdx index) const noexcept {
			return ids_[index];
		}

		void clear() noexcept {
			data_.clear();
			ids_.clear();
			std::fill(mask_.begin(), mask_.end(), INVALID_LOCAL_ID);
			next_ = 0;
		}

		void reserve(size_t n) {
			data_.reserve(n);
			ids_.reserve(n);
		}

	private:
		std::vector<T>        data_;
		std::vector<LocalId>  ids_;
		std::vector<LocalIdx> mask_; // global_id -> local dense index
		LocalId next_ = 0;           // next empty ID
	};


	// --- United Storage ---
	//  Pool sizes follow Euler's formula for a convex polyhedron:
	//    V - E + F = 2,  E ≈ 3V/2.
	//  For F = 128 faces (a dense 100-neighbour shell): V ≈ 258, E ≈ 384.
	//  The pools below give a ×2 safety margin.
	struct PolyhedronData {
		Storage<Vertex, 1024> vert;
		Storage<HalfEdge, 2048> edge;
		Storage<Polygon, 1024> face;

		void clear() noexcept {
			vert.clear(); edge.clear(); face.clear();
		}
		void reserve(size_t nv, size_t ne, size_t nf) {
			vert.reserve(nv); edge.reserve(ne); face.reserve(nf);
		}

		// Add a fresh vertex at position p with the cell center given, so
		// the cached squared distance can be computed once here instead of
		// on every call to max_vertex_dist_sq.
		LocalId add_vertex(const PointType& p, const PointType& center) {
			Vertex v{};
			v.pos = p;
			v.dist_sq = static_cast<float>((p - center).rSq());
			v.edge_hint = INVALID_LOCAL_ID;
			v._flags = 0;
			const LocalIdx idx = vert.add(v);
			return static_cast<LocalId>(idx);
		}

		// Copy an existing vertex verbatim (including dist_sq).
		LocalId copy_vertex(const Vertex& src) {
			const LocalIdx idx = vert.add(src);
			return static_cast<LocalId>(idx);
		}

		LocalId add_halfedge(LocalId origin, LocalId poly, LocalId next, LocalId twin) {
			HalfEdge h{}; h.origin_vertex_id = origin;
			h.polygon_id = poly; h.next_edge_id = next; h.twin_edge_id = twin;
			const LocalIdx idx = edge.add(h);
			return static_cast<LocalId>(idx);
		}
		LocalId add_face(LocalId first_he, AtomId other_cell, ShiftCode shift) {
			Polygon p{}; p.first_edge_id = first_he;
			p.other_cell = other_cell; p.other_shift = shift;
			const LocalIdx idx = face.add(p);
			return static_cast<LocalId>(idx);
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
		using PlaneType = geometry::Plane<InternalFloatType>;

		static constexpr InternalFloatType EPS = static_cast<InternalFloatType>(1.0e-10);

		// --- Triangle area from two vectors sharing a common origin ---
		static inline InternalFloatType triangle_area(const PointType& v1,
													  const PointType& v2) noexcept {
			return PointType::Vector(v1, v2).r() * InternalFloatType(0.5);
		}

		// --- Polygon normal via Newell's method (robust for concave polygons) ---
		static inline PointType polygon_normal(const PolyhedronData& d,
											   LocalId first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			PointType v0 = d.vert.get_by_id(e0.origin_vertex_id).pos;

			LocalId next_he = e0.next_edge_id;
			PointType v1 = d.vert.get_by_id(d.edge.get_by_id(next_he).origin_vertex_id).pos;

			InternalFloatType nx = 0, ny = 0, nz = 0;

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

			const InternalFloatType len = std::sqrt(nx * nx + ny * ny + nz * nz);
			if (len > EPS) {
				nx /= len; ny /= len; nz /= len;
			}
			return PointType(nx, ny, nz);
		}

		// --- Polygon area (fan triangulation over a half-edge cycle) ---
		static inline InternalFloatType polygon_area(const PolyhedronData& d,
													 LocalId first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const PointType p0 = d.vert.get_by_id(e0.origin_vertex_id).pos;

			const HalfEdge& e1 = d.edge.get_by_id(e0.next_edge_id);
			PointType v1 = d.vert.get_by_id(e1.origin_vertex_id).pos - p0;

			LocalId e2 = e1.next_edge_id;
			InternalFloatType area = InternalFloatType(0);

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
											  LocalId first_he) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const LocalId e1_id = e0.next_edge_id;
			const HalfEdge& e1 = d.edge.get_by_id(e1_id);
			const LocalId e2_id = e1.next_edge_id;
			return PlaneType(
				d.vert.get_by_id(e0.origin_vertex_id).pos,
				d.vert.get_by_id(e1.origin_vertex_id).pos,
				d.vert.get_by_id(d.edge.get_by_id(e2_id).origin_vertex_id).pos);
		}

		// --- Solid angle of a polygon face at point `origin` (van Oosterom-Strackee) ---
		static inline InternalFloatType polygon_solid_angle(const PolyhedronData& d,
															LocalId first_he,
															const PointType& origin) noexcept {
			const HalfEdge& e0 = d.edge.get_by_id(first_he);
			const PointType v0 = d.vert.get_by_id(e0.origin_vertex_id).pos - origin;
			const InternalFloatType r0 = v0.r();

			LocalId next_he = e0.next_edge_id;
			const PointType v1 = d.vert.get_by_id(
				d.edge.get_by_id(next_he).origin_vertex_id).pos - origin;

			PointType         v_prev = v1;
			InternalFloatType r_prev = v1.r();
			InternalFloatType dot0_prev = PointType::Scalar(v0, v1);

			next_he = d.edge.get_by_id(next_he).next_edge_id;

			InternalFloatType angle = InternalFloatType(0);

			while (next_he != first_he) {
				const HalfEdge& he = d.edge.get_by_id(next_he);
				const PointType v_curr = d.vert.get_by_id(he.origin_vertex_id).pos - origin;

				const InternalFloatType r_curr = v_curr.r();
				const InternalFloatType dot0_curr = PointType::Scalar(v0, v_curr);

				const InternalFloatType cross = std::abs(
					PointType::Scalar(v0, PointType::Vector(v_prev, v_curr)));

				const InternalFloatType denom =
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

			return angle * InternalFloatType(2);
		}

		// --- Volume: sum of pyramid contributions per face ---
		//   volume = sum over faces: dist(center, face_plane) * face_area / 3
		static inline InternalFloatType cell_volume(const PolyhedronData& d,
													const PointType& center) noexcept {
			const Polygon* polys = d.face.data();
			const size_t   nf = d.face.size();

			InternalFloatType volume = InternalFloatType(0);
			for (size_t i = 0; i < nf; ++i) {
				const LocalId fh = polys[i].first_edge_id;
				const InternalFloatType area = polygon_area(d, fh);
				const auto              pl = polygon_plane(d, fh);
				volume += area * std::abs(pl.distance(center));
			}
			return volume * InternalFloatType(1.0 / 3.0);
		}
	};

} // namespace cpplib::voronoi

// ============================================================================
//  LEVEL 3: Topological Context (Mid-Level Manager)
// ============================================================================

namespace cpplib::voronoi {

// ============================================================================
//  VoronoiCell
//
//  One Voronoi polyhedron = metadata + private topology.
//  Fully self-contained: no references to other cells, except
//  Polygon::other_cell / other_shift — those are AtomId keys into
//  VoronoiContext::cells, not pointers.
// ============================================================================
	struct VoronoiCell {
		Polyhedron     meta;   // {center, volume}
		PolyhedronData topo;   // {vert, edge, face} of this cell

		void clear() noexcept {
			meta = Polyhedron{};
			topo.clear();
		}
		void reserve(size_t nv, size_t ne, size_t nf) {
			topo.reserve(nv, ne, nf);
		}
	};


	// ============================================================================
	//  VoronoiContext
	//
	//  Dense vector of independent cells (index = atom_id) + shared read-only
	//  spatial index SG<T>.
	// ============================================================================
	struct VoronoiContext {
		using GridType = geometry::SG<InternalFloatType>;

		// --- Cells: index = atom_id, dense, 1:1 with atoms ---
		std::vector<VoronoiCell> cells;

		// --- Shared spatial index. Built once; during the parallel phase
		//     it is only read (all query methods are const). ---
		GridType grid;

		// --- Configuration copied from the pipeline ---
		uint32_t max_neighbors = 12;             // reserved for future tuning
		bool     use_pbc = true;
		int      sup_x = 1, sup_y = 1, sup_z = 1; // reserved for supercell support

		void clear() noexcept {
			cells.clear();
			grid = GridType{};
		}

		VoronoiCell& cell(size_t i)       noexcept {
			return cells[i];
		}
		const VoronoiCell& cell(size_t i) const noexcept {
			return cells[i];
		}

		size_t n_cells() const noexcept {
			return cells.size();
		}

		// Access to the shared index. Const only — guarantees nobody
		// mutates it during the parallel phase.
		const GridType& spatial_grid() const noexcept {
			return grid;
		}

		// ------------------------------------------------------------------
		//  Build the spatial index. Must run BEFORE the parallel cell
		//  construction: during that phase the grid is read-only.
		//
		//  points_frac  — fractional atom coordinates in [0,1).
		//  use_pbc_flag — passed to the SG constructor; SG has no setter,
		//                 so the grid is recreated on every rebuild.
		// ------------------------------------------------------------------
		void rebuildGrid(const std::vector<PointType>& points_frac,
						 const UnitCell& cell_def,
						 InternalFloatType cutoff,
						 bool use_pbc_flag)
		{
			grid = GridType(cell_def, use_pbc_flag);
			grid.updateCellAndPoints(points_frac, cutoff, cell_def);
			use_pbc = use_pbc_flag;
		}

		// ------------------------------------------------------------------
		//  Recompute metrics for a single cell. Fully local, safe to call
		//  from a parallel loop: touches only the passed cell.
		// ------------------------------------------------------------------
		static void recompute_cell_metrics(VoronoiCell& c) noexcept {
			c.meta.volume = MathEngine::cell_volume(c.topo, c.meta.center);
		}

		struct Stats {
			size_t n_cells = 0;
			size_t n_vertices = 0;
			size_t n_edges = 0;
			size_t n_faces = 0;
			InternalFloatType total_volume = InternalFloatType(0);
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
// ============================================================================

namespace cpplib::voronoi {

	class VoronoiFused;   // forward declaration

// ============================================================================
//  WorkerPool
//
//  Portable worker pool built on std::jthread. Used instead of
//  std::execution::par because:
//
//    - <execution> is not available or not usable on every toolchain
//      (GCC without TBB, libc++ historically). On those platforms
//      std::execution::par silently degrades to sequential execution,
//      which is worse than a compile-time error.
//    - jthread gives explicit control over the thread count.
//    - Per-thread scratch state can be handed to each worker for the
//      whole lifetime of one parallel_for call — critical for the
//      Voronoi per-cell construction, which uses ~10 scratch vectors
//      per cell and would otherwise thrash the allocator.
//
//  The pool does not attempt work stealing. Work is distributed via a
//  single atomic counter; for the near-uniform per-cell cost in this
//  codebase that is both simpler and faster than a stealing queue.
// ============================================================================
	class WorkerPool {
	public:
		// threads == 0 -> use hardware_concurrency (min 1).
		explicit WorkerPool(size_t threads = 0)
			: n_threads_(threads?threads
						 :std::max<size_t>(1, std::thread::hardware_concurrency())) {
		}

		size_t size() const noexcept {
			return n_threads_;
		}

		// Parallel loop over [0, n).
		//   init  -> produces a fresh Scratch value for one worker
		//   fn    -> callback (index, scratch&), must be thread-safe across indices
		//
		// Falls back to a sequential pass if the pool has one thread or
		// the workload is small — the cost of spawning jthreads is not
		// worth it below a few hundred microseconds of total work.
		template <typename Scratch, typename Init, typename Fn>
		void parallel_for(size_t n, const Init& init, const Fn& fn) const {
			if (n == 0) return;

			if (n_threads_ <= 1 || n < kMinParallelWork) {
				Scratch scratch = init();
				for (size_t i = 0; i < n; ++i)
					fn(i, scratch);
				return;
			}

			// Per-worker scratch: allocated once, reused across indices.
			std::vector<Scratch> scratch;
			scratch.reserve(n_threads_);
			for (size_t t = 0; t < n_threads_; ++t)
				scratch.emplace_back(init());

			std::atomic<size_t> next{0};

			auto worker = [&](size_t tid) noexcept {
				Scratch& local = scratch[tid];
				for (;;) {
					const size_t i = next.fetch_add(1, std::memory_order_relaxed);
					if (i >= n) return;
					fn(i, local);
				}
				};

				// Spawn n_threads_ - 1 workers; the calling thread also runs.
			std::vector<std::jthread> workers;
			workers.reserve(n_threads_ - 1);
			for (size_t t = 1; t < n_threads_; ++t)
				workers.emplace_back(worker, t);

			worker(0);

			// jthread destructors join.
		}

	private:
		// Below this iteration count, a sequential pass wins.
		static constexpr size_t kMinParallelWork = 64;

		size_t n_threads_;
	};


// ============================================================================
//  VoronoiPipeline
//
//  Orchestrator. The cutoff is supplied externally: the caller knows
//  the physics (6 A typical, 10-30 A large complexes, up to the cell
//  size for low-density structures).
//
//  Parallelism:
//    - the spatial index is built once, then read-only;
//    - per-cell construction writes only into cells[i];
//    - metrics are recomputed per-cell, also independently.
//
//  The parallel loops go through a caller-supplied WorkerPool. If pool
//  is null, everything runs sequentially — useful for debugging and for
//  small inputs where thread startup would dominate.
//
//  Per-cell construction uses WorkerScratch, a bundle of reusable
//  vectors. Each worker gets one instance per parallel_for call, so the
//  hot path performs zero allocations after the first few cells warm up.
// ============================================================================
	struct VoronoiPipeline {

		struct Config {
			// Reserved for future tuning knobs.
		};

		Config         config;
		VoronoiContext context;

		// Optional. When non-null, build()/update() distribute work
		// across the pool. Owned by the caller; the pipeline does not
		// take ownership.
		WorkerPool* pool = nullptr;

		// ==================================================================
		//  Public entry points
		// ==================================================================

		// ------------------------------------------------------------------
		//  Full build. Caller provides the cutoff in Cartesian units.
		// ------------------------------------------------------------------
		template <typename ExtT = ExternalFloatType>
		inline void build(const std::vector<geometry::Point<ExtT>>& points_frac,
				   const geometry::Cell<ExtT>& cell_def,
				   ExtT cutoff,
				   bool use_pbc_flag = true)
		{
			using ExtPoint = geometry::Point<ExtT>;
			using ExtCell = geometry::Cell<ExtT>;

			if constexpr (std::is_same_v<ExtT, InternalFloatType>) {
				// Same type — forward directly, no copy.
				build_internal(points_frac, cell_def, cutoff, use_pbc_flag);
			} else {
				// Convert once, then run the internal pipeline.
				std::vector<PointType> internal_pts(points_frac.size());
				for (size_t i = 0; i < points_frac.size(); ++i)
					internal_pts[i] = PointType(points_frac[i]);

				UnitCell internal_cell = to_internal_cell(cell_def);
				build_internal(internal_pts, internal_cell, cutoff, use_pbc_flag);
			}
		}

		// Rebuild a Cell<double> from a Cell<ExtT> by going through the public
		// parameters: lattice lengths and angles. The frac->cart matrix is
		// recomputed from these, so precision matches the internal type.
		template <typename ExtT>
		static UnitCell to_internal_cell(const geometry::Cell<ExtT>& ext) {
			return UnitCell(
				static_cast<InternalFloatType>(ext.lat_dir(0)),
				static_cast<InternalFloatType>(ext.lat_dir(1)),
				static_cast<InternalFloatType>(ext.lat_dir(2)),
				static_cast<InternalFloatType>(ext.getAngleGrad(0)),
				static_cast<InternalFloatType>(ext.getAngleGrad(1)),
				static_cast<InternalFloatType>(ext.getAngleGrad(2)),
				/*is_grad=*/true);
		}

		void build_internal(const std::vector<PointType>& points_frac,
				   const UnitCell& cell_def,
				   InternalFloatType cutoff,
				   bool use_pbc_flag = true)
		{
			assert(!points_frac.empty());
			assert(points_frac.size() <= static_cast<size_t>(INVALID_ATOM_ID));

			context.clear();
			const size_t n = points_frac.size();
			context.cells.resize(n);

			// Spatial index: single-threaded, sequential.
			context.rebuildGrid(points_frac, cell_def, cutoff, use_pbc_flag);

			// Per-cell construction: parallel over independent cells.
			for_each_cell(n, [&](size_t i, WorkerScratch& scratch) {
				build_single_cell(context.cells[i], static_cast<AtomId>(i),
								  points_frac, cutoff, scratch);
						  });

						  // Metric pass: parallel over independent cells.
			if (pool) {
				pool->parallel_for<EmptyScratch>(
					n,
					[] { return EmptyScratch{}; },
					[&](size_t i, EmptyScratch&) {
						VoronoiContext::recompute_cell_metrics(context.cells[i]);
					});
			} else {
				EmptyScratch scratch;
				for (size_t i = 0; i < n; ++i)
					VoronoiContext::recompute_cell_metrics(context.cells[i]);
			}
		}

		// ------------------------------------------------------------------
		//  Update: only atom positions changed. Topology kept.
		// ------------------------------------------------------------------
		template <typename ExtT = ExternalFloatType>
		void update(const std::vector<geometry::Point<ExtT>>& points_frac) {
			if constexpr (std::is_same_v<ExtT, InternalFloatType>) {
				update_internal(points_frac);
			} else {
				std::vector<PointType> internal_pts(points_frac.size());
				for (size_t i = 0; i < points_frac.size(); ++i)
					internal_pts[i] = PointType(points_frac[i]);
				update_internal(internal_pts);
			}
		}
		
		void update_internal(const std::vector<PointType>& points_frac) {
			const size_t n = points_frac.size();
			const MatrixType& FtoC = context.spatial_grid().frac_to_cart();

			if (pool) {
				pool->parallel_for<EmptyScratch>(
					n,
					[] { return EmptyScratch{}; },
					[&](size_t i, EmptyScratch&) {
						context.cells[i].meta.center = FtoC * points_frac[i];
						VoronoiContext::recompute_cell_metrics(context.cells[i]);
					});
			} else {
				for (size_t i = 0; i < n; ++i) {
					context.cells[i].meta.center = FtoC * points_frac[i];
					VoronoiContext::recompute_cell_metrics(context.cells[i]);
				}
			}
		}

		size_t                n_cells() const noexcept {
			return context.n_cells();
		}
		VoronoiContext::Stats stats()   const noexcept {
			return context.get_stats();
		}
		const VoronoiCell& cell(size_t i) const noexcept {
			return context.cell(i);
		}
		const VoronoiContext& ctx()     const noexcept {
			return context;
		}

		VoronoiFused fuse() const;

	private:
		// ==================================================================
		//  Worker scratch
		//
		//  One instance per worker thread. All vectors keep their capacity
		//  across cells, so the hot path performs zero allocations after
		//  the first few cells warm up. Never shared between threads.
		// ==================================================================
		struct WorkerScratch {
			std::vector<geometry::SG<InternalFloatType>::Neighbor> cand;
			std::vector<uint8_t> state;              // VertexSide as uint8_t
			std::vector<LocalId> vmap;               // old vertex -> new vertex
			std::vector<LocalId> cross_vert;         // one entry per old half-edge
			std::vector<LocalId> plane_verts;        // intersection vertex IDs
			std::vector<size_t>  he_li;              // half-edge local indices
			std::vector<size_t>  v_li;               // vertex local indices
			std::vector<LocalId> new_cycle;          // scratch for face rebuild
			std::vector<LocalId> he_ids;             // scratch for add_face_cycle
			std::vector<LocalId> twin_a;             // rebuild_twins: origin per half-edge
			std::vector<LocalId> twin_b;             // rebuild_twins: next-origin per half-edge
		};

		// Placeholder for parallel loops that do not need scratch.
		struct EmptyScratch {
		};

// ==================================================================
//  Dispatch one parallel_for through the pool, or run sequentially.
// ==================================================================
		template <typename Fn>
		void for_each_cell(size_t n, const Fn& fn) const {
			if (pool) {
				pool->parallel_for<WorkerScratch>(
					n,
					[] { return WorkerScratch{}; },
					fn);
			} else {
				WorkerScratch scratch;
				for (size_t i = 0; i < n; ++i)
					fn(i, scratch);
			}
		}

		// ==================================================================
		//  Vertex classification result for one clipping pass.
		// ==================================================================
		enum class VertexSide : uint8_t {
			Delete = 0,  // plane.side < -eps — will be removed
			Keep = 1,  // plane.side >  eps — kept as-is
			OnPlane = 2,  // |plane.side| <= eps — kept and snapped
		};

		// ==================================================================
		//  Per-cell construction — top-level driver.
		//
		//  Steps:
		//    1. Query the shared index for candidates within cutoff.
		//    2. Sort candidates by squared distance, nearest first.
		//    3. Seed with an axis-aligned box.
		//    4. Clip iteratively:
		//         early-exit if nb.dist_sq > 4 * max_vertex_dist_sq;
		//         clip by the perpendicular bisector plane;
		//         if the clip removed any vertex, refresh max_vertex_dist_sq.
		//
		//  The bisector between the cell center c and neighbour n sits at
		//  distance |c - n| / 2 from c. If that distance exceeds the
		//  current cell's farthest-vertex radius, the plane cannot cut
		//  the cell. Solving |c - n| / 2 < r_max gives the bound used.
		// ==================================================================
		void build_single_cell(VoronoiCell& cell, AtomId atom_id,
							   const std::vector<PointType>& points_frac,
							   InternalFloatType cutoff,
							   WorkerScratch& scratch) const
		{
			cell.clear();
			cell.reserve(/*nv*/ 64, /*ne*/ 192, /*nf*/ 64);

			const auto& grid = context.spatial_grid();
			const MatrixType& FtoC = grid.frac_to_cart();
			cell.meta.center = FtoC * points_frac[atom_id];

			// --- 1. Collect candidates (reuse scratch buffer) -----------
			auto& cand = scratch.cand;
			cand.clear();
			grid.for_each_neighbor(static_cast<uint32_t>(atom_id),
								   [&](const geometry::SG<InternalFloatType>::Neighbor& nb) {
									   cand.push_back(nb);
								   });

							   // --- 2. Sort by squared distance -----------------------------
			std::sort(cand.begin(), cand.end(),
					  [](const auto& a, const auto& b) { return a.dist_sq < b.dist_sq; });

			// --- 3. Seed with an axis-aligned box ------------------------
			seed_box(cell, cutoff, scratch);
			float max_vd_sq = max_vertex_dist_sq(cell);

			// --- 4. Iterative clipping -----------------------------------
			const auto& geom = grid.geometry();
			for (const auto& nb : cand) {
				if (nb.dist_sq > InternalFloatType(4) * static_cast<InternalFloatType>(max_vd_sq)) break;

				const PointType nbr_cart =
					FtoC * points_frac[nb.idx] +
					geom.shift_cart[nb.shift.get_code()];

				// Bisector plane: perpendicular to (center - nbr), through
				// the midpoint, oriented so the center is on the +side.
				const PointType mid = (cell.meta.center + nbr_cart) * InternalFloatType(0.5);
				const PointType normal = cell.meta.center - nbr_cart;
				const geometry::Plane<InternalFloatType> bisector(mid, normal);

				if (clip_by_plane(cell, bisector,
								  static_cast<AtomId>(nb.idx),
								  nb.shift, scratch))
					max_vd_sq = max_vertex_dist_sq(cell);
			}
		}

		// ==================================================================
		//  clip_by_plane — pipeline entry point.
		//
		//  Delegates each stage to a step helper. All scratch vectors come
		//  from WorkerScratch; no allocations on the hot path.
		// ==================================================================
		bool clip_by_plane(VoronoiCell& cell,
						   const geometry::Plane<InternalFloatType>& plane,
						   AtomId other_cell,
						   ShiftCode other_shift,
						   WorkerScratch& scratch) const
		{
			const PolyhedronData& old = cell.topo;

			// Step 1 — classify; early exit if the plane cuts nothing.
			if (!classify_vertices(old, plane, scratch.state))
				return false;

			// Step 2 — fresh topology; copy surviving vertices.
			PolyhedronData fresh;
			fresh.reserve(old.vert.size() + 16,
						  old.edge.size() + 32,
						  old.face.size() + 4);

			copy_surviving_vertices(old, scratch.state, fresh, scratch.vmap);

			// Step 3 — compute intersections on crossing half-edges.
			find_plane_crossings(old, scratch.state, plane, cell.meta.center,
								 fresh, scratch.cross_vert, scratch.plane_verts);

			// Step 4 — rebuild surviving faces.
			rebuild_surviving_faces(old, scratch.state, scratch.vmap,
									scratch.cross_vert, fresh, scratch);

			// Step 5 — add the cutting face.
			add_cutting_face(fresh, plane, cell.meta.center,
							 scratch.plane_verts,
							 other_cell, other_shift, scratch);

			// Step 6 — restore twin links.
			rebuild_twins(fresh, scratch);

			// Step 7 — commit.
			cell.topo = std::move(fresh);
			return true;
		}

		// ==================================================================
		//  clip_by_plane — step 1
		//
		//  Classify every vertex relative to the plane. Vertices whose
		//  |side| <= ON_PLANE_EPS are treated as ON_PLANE and survive;
		//  this reduces the chance that the same physical vertex is
		//  deleted in one cell and kept in its twin (numerical drift).
		//
		//  Returns true if at least one vertex was classified as Delete.
		// ==================================================================
		static bool classify_vertices(const PolyhedronData& old,
									  const geometry::Plane<InternalFloatType>& plane,
									  std::vector<uint8_t>& state)
		{
			constexpr InternalFloatType ON_PLANE_EPS = InternalFloatType(1e-8);

			const size_t nv = old.vert.size();
			state.resize(nv);

			bool any_deleted = false;
			for (size_t i = 0; i < nv; ++i) {
				const InternalFloatType s = plane.side(old.vert.data()[i].pos);
				if (s < -ON_PLANE_EPS) {
					state[i] = static_cast<uint8_t>(VertexSide::Delete);
					any_deleted = true;
				} else if (s < ON_PLANE_EPS) {
					state[i] = static_cast<uint8_t>(VertexSide::OnPlane);
				} else {
					state[i] = static_cast<uint8_t>(VertexSide::Keep);
				}
			}
			return any_deleted;
		}

		// ==================================================================
		//  clip_by_plane — step 2
		//
		//  Copy all surviving vertices from `old` into `fresh`, building
		//  the old -> new vertex ID map along the way. Deleted vertices
		//  leave INVALID_LOCAL_ID in the map.
		//
		//  The whole Vertex is copied verbatim, including the cached
		//  dist_sq field — no recomputation needed since the cell center
		//  is unchanged across a clip.
		// ==================================================================
		static void copy_surviving_vertices(const PolyhedronData& old,
											const std::vector<uint8_t>& state,
											PolyhedronData& fresh,
											std::vector<LocalId>& vmap)
		{
			const size_t nv = old.vert.size();
			vmap.assign(nv, INVALID_LOCAL_ID);

			const auto del = static_cast<uint8_t>(VertexSide::Delete);
			for (size_t i = 0; i < nv; ++i) {
				if (state[i] != del)
					vmap[i] = fresh.copy_vertex(old.vert.data()[i]);
			}
		}

		// ==================================================================
		//  clip_by_plane — step 3
		//
		//  Find every crossing half-edge (KEEP/ON_PLANE -> DELETE), compute
		//  the intersection with the plane once per undirected edge, and
		//  record the fresh vertex in both directed slots of cross_vert.
		//
		//  plane_verts collects the unique intersection vertex IDs — they
		//  are the boundary of the new cutting face.
		// ==================================================================
		static void find_plane_crossings(const PolyhedronData& old,
										 const std::vector<uint8_t>& state,
										 const geometry::Plane<InternalFloatType>& plane,
										 const PointType& cell_center,
										 PolyhedronData& fresh,
										 std::vector<LocalId>& cross_vert,
										 std::vector<LocalId>& plane_verts)
		{
			const size_t ne = old.edge.size();
			cross_vert.assign(ne, INVALID_LOCAL_ID);
			plane_verts.clear();
			plane_verts.reserve(16);

			const auto del = static_cast<uint8_t>(VertexSide::Delete);

			for (size_t hi = 0; hi < ne; ++hi) {
				if (cross_vert[hi] != INVALID_LOCAL_ID) continue;

				const HalfEdge& h = old.edge.data()[hi];
				const size_t o_li = old.vert.index_of(h.origin_vertex_id);
				const size_t n_li = old.vert.index_of(
					old.edge.get_by_id(h.next_edge_id).origin_vertex_id);

				if (state[o_li] == del) continue;
				if (state[n_li] != del) continue;

				const PointType& pa = old.vert.data()[o_li].pos;
				const PointType& pb = old.vert.data()[n_li].pos;
				const InternalFloatType sa = plane.side(pa);
				const InternalFloatType sb = plane.side(pb);
				const InternalFloatType t = sa / (sa - sb);
				const PointType P = pa + (pb - pa) * t;

				const LocalId vid = fresh.add_vertex(P, cell_center);
				cross_vert[hi] = vid;
				plane_verts.push_back(vid);

				const LocalId twin = h.twin_edge_id;
				if (twin != INVALID_LOCAL_ID) {
					const size_t twin_li = old.edge.index_of(twin);
					cross_vert[twin_li] = vid;
				}
			}
		}

		// ==================================================================
		//  clip_by_plane — step 4
		//
		//  Rebuild every surviving face. For each face, walk its old
		//  half-edge cycle and emit a fresh one:
		//    - KEEP / ON_PLANE vertices are copied through vmap;
		//    - at each KEEP <-> DELETE transition, the precomputed
		//      intersection vertex is inserted.
		//
		//  Faces entirely on the DELETE side are dropped.
		//
		//  Scratch buffers are reused across faces to avoid per-face
		//  allocations on the hot path.
		// ==================================================================
		static void rebuild_surviving_faces(const PolyhedronData& old,
											const std::vector<uint8_t>& state,
											const std::vector<LocalId>& vmap,
											const std::vector<LocalId>& cross_vert,
											PolyhedronData& fresh,
											WorkerScratch& scratch)
		{
			auto& he_li = scratch.he_li;
			auto& v_li = scratch.v_li;
			auto& new_cycle = scratch.new_cycle;

			const auto del = static_cast<uint8_t>(VertexSide::Delete);

			for (size_t fi = 0; fi < old.face.size(); ++fi) {
				const Polygon& old_f = old.face.data()[fi];

				he_li.clear();
				{
					LocalId cur = old_f.first_edge_id;
					do {
						he_li.push_back(old.edge.index_of(cur));
						cur = old.edge.data()[he_li.back()].next_edge_id;
					} while (cur != old_f.first_edge_id);
				}

				v_li.resize(he_li.size());
				for (size_t k = 0; k < he_li.size(); ++k)
					v_li[k] = old.vert.index_of(
						old.edge.data()[he_li[k]].origin_vertex_id);

				bool any_survivor = false;
				for (size_t v:v_li)
					if (state[v] != del) {
						any_survivor = true; break;
					}
				if (!any_survivor) continue;

				const size_t n = v_li.size();
				new_cycle.clear();
				new_cycle.reserve(n + 2);

				for (size_t k = 0; k < n; ++k) {
					const size_t v_here = v_li[k];
					const size_t v_prev = v_li[(k + n - 1) % n];
					const size_t v_next = v_li[(k + 1) % n];

					if (state[v_here] == del) continue;

					if (state[v_prev] == del) {
						const size_t he_prev = he_li[(k + n - 1) % n];
						const LocalId cv = cross_vert[he_prev];
						if (cv != INVALID_LOCAL_ID) new_cycle.push_back(cv);
					}

					new_cycle.push_back(vmap[v_here]);

					if (state[v_next] == del) {
						const LocalId cv = cross_vert[he_li[k]];
						if (cv != INVALID_LOCAL_ID) new_cycle.push_back(cv);
					}
				}

				if (new_cycle.size() < 3) continue;
				add_face_cycle(fresh, old_f.other_cell, old_f.other_shift,
							   new_cycle, scratch.he_ids);
			}
		}

		// ==================================================================
		//  clip_by_plane — step 5
		//
		//  Add the cutting face itself. Its vertices are the intersection
		//  points collected in step 3; they are sorted by polar angle
		//  around the reversed plane normal so the face is traversed CCW
		//  as seen from outside the cell (valid for the convex cross-section).
		//
		//  The vertices of the cutting face are owned by the same fresh
		//  topology that already holds them — plane_verts is a list of
		//  references, not a copy.
		// ==================================================================
		static void add_cutting_face(PolyhedronData& fresh,
									 const geometry::Plane<InternalFloatType>& plane,
									 const PointType& cell_center,
									 std::vector<LocalId>& plane_verts,
									 AtomId other_cell,
									 ShiftCode other_shift,
									 WorkerScratch& scratch)
		{
			if (plane_verts.size() < 3) return;

			PointType centroid(0, 0, 0);
			for (LocalId v:plane_verts)
				centroid = centroid + fresh.vert.get_by_id(v).pos;
			centroid = centroid / static_cast<InternalFloatType>(plane_verts.size());

			PointType nrm(-plane.a[0], -plane.a[1], -plane.a[2]);
			const InternalFloatType nl = nrm.r();
			if (nl > 1e-12) nrm = nrm / nl;

			PointType up(0, 0, 1);
			if (std::abs(nrm[2]) > InternalFloatType(0.9)) up = PointType(1, 0, 0);
			PointType ex = PointType::Vector(up, nrm);
			const InternalFloatType exl = ex.r();
			if (exl > 1e-12) ex = ex / exl;
			const PointType ey = PointType::Vector(nrm, ex);

			std::sort(plane_verts.begin(), plane_verts.end(),
					  [&](LocalId a, LocalId b) {
						  const PointType da = fresh.vert.get_by_id(a).pos - centroid;
						  const PointType db = fresh.vert.get_by_id(b).pos - centroid;
						  const InternalFloatType aa = std::atan2(
							  PointType::Scalar(da, ey), PointType::Scalar(da, ex));
						  const InternalFloatType ab = std::atan2(
							  PointType::Scalar(db, ey), PointType::Scalar(db, ex));
						  return aa < ab;
					  });

			add_face_cycle(fresh, other_cell, other_shift,
						   plane_verts, scratch.he_ids);
		}

		// ==================================================================
		//  seed_box
		//
		//  Initial polyhedron: an axis-aligned cube centered on the atom,
		//  six faces oriented CCW viewed from outside. Each shared edge
		//  appears in opposite directions in the two faces that share it —
		//  verified by the table below.
		//
		//  Half-size must be >= cutoff so the box contains the final
		//  Voronoi cell. Face shift codes are placeholders — every face
		//  is expected to be cut at least once by a real bisector plane.
		// ==================================================================
		void seed_box(VoronoiCell& cell, InternalFloatType half_size,
					  WorkerScratch& scratch) const
		{
			const PointType c = cell.meta.center;
			const InternalFloatType h = half_size;

			const PointType pts[8] = {
				c + PointType(-h, -h, -h),  // 0
				c + PointType(+h, -h, -h),  // 1
				c + PointType(+h, +h, -h),  // 2
				c + PointType(-h, +h, -h),  // 3
				c + PointType(-h, -h, +h),  // 4
				c + PointType(+h, -h, +h),  // 5
				c + PointType(+h, +h, +h),  // 6
				c + PointType(-h, +h, +h)   // 7
			};

			LocalId vids[8];
			for (int i = 0; i < 8; ++i)
				vids[i] = cell.topo.add_vertex(pts[i], c);

			static constexpr int face_verts[6][4] = {
				{1, 5, 6, 2}, {0, 3, 7, 4}, {2, 6, 7, 3},
				{0, 4, 5, 1}, {4, 7, 6, 5}, {0, 1, 2, 3}
			};
			static constexpr uint8_t face_shifts[6] = {14, 12, 16, 10, 22, 4};

			std::vector<LocalId> cycle(4);
			for (int f = 0; f < 6; ++f) {
				cycle[0] = vids[face_verts[f][0]];
				cycle[1] = vids[face_verts[f][1]];
				cycle[2] = vids[face_verts[f][2]];
				cycle[3] = vids[face_verts[f][3]];
				add_face_cycle(cell.topo, /*other_cell*/ 0,
							   ShiftCode(face_shifts[f]), cycle, scratch.he_ids);
			}

			rebuild_twins(cell.topo, scratch);
		}

		// ==================================================================
		//  max_vertex_dist_sq
		//
		//  Maximum squared distance from the cell center to any vertex.
		//  Drives the early-exit bound in build_single_cell.
		//
		//  Reads the cached dist_sq from each Vertex instead of recomputing
		//  (pos - center).rSq(): the scan becomes a 4-byte-per-element
		//  comparison against a flat array, roughly 6x less memory traffic
		//  and 10x fewer arithmetic ops than the recomputing version.
		// ==================================================================
		float max_vertex_dist_sq(const VoronoiCell& cell) const {
			const Vertex* vd = cell.topo.vert.data();
			const size_t  nv = cell.topo.vert.size();

			float best = 0.0f;
			for (size_t i = 0; i < nv; ++i)
				if (vd[i].dist_sq > best) best = vd[i].dist_sq;
			return best;
		}

		// ==================================================================
		//  Helpers — half-edge primitives
		// ==================================================================

		// ------------------------------------------------------------------
		//  add_face_cycle
		//
		//  Emit a closed polygon cycle. next_edge_id and polygon_id are
		//  set for every half-edge. twin_edge_id stays INVALID —
		//  rebuild_twins() fixes the twin links once the whole topology
		//  is stable.
		// ------------------------------------------------------------------
		static LocalId add_face_cycle(PolyhedronData& topo,
									  AtomId other_cell,
									  ShiftCode shift,
									  const std::vector<LocalId>& verts,
									  std::vector<LocalId>& he_ids_scratch)
		{
			const size_t n = verts.size();
			assert(n >= 3);

			he_ids_scratch.resize(n);
			for (size_t i = 0; i < n; ++i)
				he_ids_scratch[i] = topo.add_halfedge(verts[i], /*poly*/ 0,
													  /*next*/ 0, INVALID_LOCAL_ID);

			for (size_t i = 0; i < n; ++i)
				topo.edge.get_by_id(he_ids_scratch[i]).next_edge_id =
				he_ids_scratch[(i + 1) % n];

			const LocalId face_id = topo.add_face(he_ids_scratch[0], other_cell, shift);
			for (LocalId h:he_ids_scratch)
				topo.edge.get_by_id(h).polygon_id = face_id;

			return face_id;
		}

		// ------------------------------------------------------------------
		//  rebuild_twins
		//
		//  Pair every half-edge (a -> b) with its twin (b -> a). Linear
		//  O(E^2) scan: for E in the 40-400 range (one Voronoi cell),
		//  this beats sort + binary search — no sort, no allocation for
		//  scratch, no unpredictable branches inside lower_bound.
		// ------------------------------------------------------------------
		static void rebuild_twins(PolyhedronData& topo, WorkerScratch& scratch) {
			const size_t ne = topo.edge.size();
			if (ne == 0) return;

			auto& ep_a = scratch.twin_a;
			auto& ep_b = scratch.twin_b;
			ep_a.resize(ne);
			ep_b.resize(ne);

			for (size_t i = 0; i < ne; ++i) {
				const HalfEdge& h = topo.edge.data()[i];
				ep_a[i] = h.origin_vertex_id;
				ep_b[i] = topo.edge.get_by_id(h.next_edge_id).origin_vertex_id;
			}

			for (size_t i = 0; i < ne; ++i)
				topo.edge.data()[i].twin_edge_id = INVALID_LOCAL_ID;

			for (size_t i = 0; i < ne; ++i) {
				HalfEdge& hi = topo.edge.data()[i];
				if (hi.twin_edge_id != INVALID_LOCAL_ID) continue;

				const LocalId a = ep_a[i];
				const LocalId b = ep_b[i];

				for (size_t j = i + 1; j < ne; ++j) {
					if (ep_a[j] != b || ep_b[j] != a) continue;
					HalfEdge& hj = topo.edge.data()[j];
					if (hj.twin_edge_id != INVALID_LOCAL_ID) continue;

					hi.twin_edge_id = topo.edge.get_id(static_cast<LocalIdx>(j));
					hj.twin_edge_id = topo.edge.get_id(static_cast<LocalIdx>(i));
					break;
				}
			}
		}
	};

} // namespace cpplib::voronoi

// ============================================================================
//  VoronoiFused — flattened mesh built from all per-cell topologies.
//
//  Layout is fully CSR: no nested vectors inside POD records, no maps.
//  The build is decomposed into five top-level steps, each of which is
//  further split into small single-purpose helpers.
//
//  Twin faces: numerical drift during clipping may give two records of
//  the same physical face slightly different vertex sets (a vertex
//  survives in one cell but is clipped in the other). Records are grouped
//  by the exact canonical key (min_atom, max_atom, shift) and their
//  vertex sets are merged within tolerance in fractional space.
// ============================================================================

namespace cpplib::voronoi {

	class VoronoiFused {
	public:
		struct EdgeFused {
			uint32_t v0 = 0;   // global vertex ID, v0 < v1
			uint32_t v1 = 0;
		};

		struct PolygonFused {
			uint32_t          vert_offset;
			uint16_t          vert_count;
			uint32_t          edge_offset;
			uint16_t          edge_count;
			AtomId            owner_atom = INVALID_ATOM_ID;
			AtomId            other_atom = INVALID_ATOM_ID;
			ShiftCode         other_shift;
			InternalFloatType area = 0;
			InternalFloatType solid_angle = 0;
		};

		struct PolyhedronFused {
			PointType         center;
			InternalFloatType volume = 0;
			uint32_t vert_offset; uint16_t vert_count;
			uint32_t edge_offset; uint16_t edge_count;
			uint32_t poly_offset; uint16_t poly_count;
		};

		std::vector<PointType>       vertices;
		std::vector<EdgeFused>       edges;
		std::vector<PolygonFused>    polygons;
		std::vector<PolyhedronFused> polyhedra;

		// Flat buffers referenced by the CSR offsets above.
		std::vector<uint32_t> poly_verts;
		std::vector<uint32_t> poly_edges;
		std::vector<uint32_t> ph_verts;
		std::vector<uint32_t> ph_edges;
		std::vector<uint32_t> ph_polys;

		struct ValidationReport {
			size_t cells = 0;
			size_t face_records = 0;
			size_t physical_faces = 0;
			size_t twin_pairs = 0;
			size_t self_faces = 0;
			size_t vertex_count_bad = 0;
			size_t vertex_set_bad = 0;
			size_t atom_pair_bad = 0;
			size_t shift_bad = 0;
			size_t twin_missing = 0;
			std::vector<std::string> messages;

			bool ok() const noexcept {
				return vertex_count_bad == 0
					&& vertex_set_bad == 0
					&& atom_pair_bad == 0
					&& shift_bad == 0
					&& twin_missing == 0;
			}
		};
		ValidationReport report;

		// ------------------------------------------------------------------
		//  Public entry point — pipeline driver.
		// ------------------------------------------------------------------
		static VoronoiFused build(const VoronoiContext& ctx);

	private:
		// ==================================================================
		//  Internal data records
		// ==================================================================

		struct FaceRec {
			AtomId   a;           // canonical min atom ID
			AtomId   b;           // canonical max atom ID
			uint8_t  s;           // canonical shift (a -> b)
			AtomId   src_cell;    // which cell this record came from
			LocalIdx src_local;   // local face index inside src_cell
		};

		struct MergedVert {
			PointType pos_cart;
			PointType pos_frac;   // wrapped to [0, 1)
			uint32_t  gid;        // global vertex ID
		};

		// ==================================================================
		//  Step 1 — merge vertices
		//
		//  Fractional coordinates are computed once per vertex and reused
		//  for both the sort key and the merge test. The merge test uses
		//  the minimum-image metric, so boundary vertices living near
		//  opposite ends of an axis still match.
		// ==================================================================
		static void merge_vertices(
			const VoronoiContext& ctx,
			const MatrixType& CtoF,
			std::vector<PointType>& out_vertices,
			std::vector<std::vector<uint32_t>>& cell_vert_gid);

		// ==================================================================
		//  Step 2 — merge edges
		//
		//  Flat array of (v0, v1, cell, local) with v0 < v1, sorted by
		//  (v0, v1). A two-pointer sweep assigns global IDs; the global
		//  edge array ends up sorted by (v0, v1).
		// ==================================================================
		static void merge_edges(
			const VoronoiContext& ctx,
			const std::vector<std::vector<uint32_t>>& cell_vert_gid,
			std::vector<EdgeFused>& out_edges,
			std::vector<std::vector<uint32_t>>& cell_edge_gid);

		// ==================================================================
		//  Step 3 — merge faces
		//
		//  Collect all face records, sort by canonical key (a, b, s), and
		//  emit one polygon per physical face. Per-group work is delegated
		//  to emit_face_group().
		// ==================================================================
		static void merge_faces(
			const VoronoiContext& ctx,
			const MatrixType& CtoF,
			const std::vector<std::vector<uint32_t>>& cell_vert_gid,
			const std::vector<std::vector<uint32_t>>& cell_edge_gid,
			VoronoiFused& out,
			std::vector<std::vector<uint32_t>>& cell_face_gid);

		// ==================================================================
		//  Step 3a — collect one FaceRec per (cell, face) pair.
		// ==================================================================
		static std::vector<FaceRec> collect_face_records(
			const VoronoiContext& ctx,
			std::vector<std::vector<uint32_t>>& cell_face_gid);

		// ==================================================================
		//  Step 3b — group scan: find the end of the (a, b, s) run.
		// ==================================================================
		static size_t find_group_end(const std::vector<FaceRec>& recs, size_t i);

		// ==================================================================
		//  Step 3c — emit one canonical polygon for a group of records that
		//  all share the same (a, b, s) key.
		//
		//  This is itself a small pipeline:
		//    1. Merge the group's vertex sets with tolerance.
		//    2. Order the merged set by polar angle around the centroid.
		//    3. Emit vertices, edges, metadata, metrics.
		//    4. Register the polygon with every source record.
		//    5. Validate group consistency.
		// ==================================================================
		static void emit_face_group(
			const VoronoiContext& ctx,
			const MatrixType& CtoF,
			const std::vector<std::vector<uint32_t>>& cell_vert_gid,
			const std::vector<std::vector<uint32_t>>& cell_edge_gid,
			VoronoiFused& out,
			std::vector<std::vector<uint32_t>>& cell_face_gid,
			const std::vector<FaceRec>& recs,
			size_t i, size_t j);

		// ==================================================================
		//  Step 3c.i — collect merged vertex set across the group.
		//
		//  Deduplication is by fractional-position matching within
		//  merge_eps, using the minimum-image metric. min_single_count /
		//  max_single_count capture the spread in per-record vertex counts
		//  for validation.
		// ==================================================================
		static void merge_group_vertices(
			const VoronoiContext& ctx,
			const MatrixType& CtoF,
			const std::vector<std::vector<uint32_t>>& cell_vert_gid,
			const std::vector<PointType>& global_vertices,
			const std::vector<FaceRec>& recs,
			size_t i, size_t j,
			std::vector<MergedVert>& merged,
			size_t& min_single_count,
			size_t& max_single_count);

		// ==================================================================
		//  Step 3c.ii — polar sort of the merged vertex set.
		//
		//  Voronoi faces are convex, so the polar order around the
		//  centroid equals the boundary order. The face normal is
		//  computed via Newell's method over the unordered set.
		// ==================================================================
		static void order_by_polar_angle(std::vector<MergedVert>& merged);

		// ==================================================================
		//  Step 3c.iii — emit edges for a face group.
		//
		//  Edges are taken from the record with the largest vertex count.
		//  If drift dropped one vertex in the twin record, the larger
		//  record still carries every edge of the merged face.
		// ==================================================================
		static void emit_polygon_edges(
			const VoronoiContext& ctx,
			const std::vector<std::vector<uint32_t>>& cell_edge_gid,
			const std::vector<FaceRec>& recs,
			size_t i, size_t j,
			PolygonFused& pg,
			std::vector<uint32_t>& poly_edges);

		// ==================================================================
		//  Step 3c.iv — emit owner/other/shift canonical metadata.
		//
		//  The (a, b, s) key already encodes a < b, so the polygon always
		//  stores the pair in that order.
		// ==================================================================
		static void emit_polygon_metadata(
			const VoronoiContext& ctx,
			const std::vector<FaceRec>& recs,
			size_t i,
			PolygonFused& pg);

		// ==================================================================
		//  Step 3c.v — face metrics.
		// ==================================================================
		static InternalFloatType compute_face_area(
			const std::vector<MergedVert>& merged);

		static InternalFloatType compute_face_solid_angle(
			const std::vector<MergedVert>& merged,
			const PointType& origin);

		// ==================================================================
		//  Step 3c.vi — group validation.
		//
		//  Records all anomalies in out.report without aborting the build.
		// ==================================================================
		static void validate_face_group(
			const VoronoiContext& ctx,
			const std::vector<FaceRec>& recs,
			size_t i, size_t j,
			const std::vector<MergedVert>& merged,
			size_t min_single_count,
			size_t max_single_count,
			uint32_t gp,
			VoronoiFused& out);

		// ==================================================================
		//  Step 4 — assemble per-cell polyhedra (CSR).
		// ==================================================================
		static void assemble_polyhedra(
			const VoronoiContext& ctx,
			const std::vector<std::vector<uint32_t>>& cell_vert_gid,
			const std::vector<std::vector<uint32_t>>& cell_edge_gid,
			const std::vector<std::vector<uint32_t>>& cell_face_gid,
			VoronoiFused& out);

		// ==================================================================
		//  Step 5 — volumes via the pyramid formula.
		//     V = (1/3) * sum_f A_f * h_f
		//  h_f is the perpendicular distance from the polyhedron centre
		//  to the plane of face f. No sign convention needed.
		// ==================================================================
		static void compute_volumes(VoronoiFused& out);

		// ==================================================================
		//  PBC-aware comparison of two fractional positions.
		// ==================================================================
		static inline bool same_frac_pbc(const PointType& a,
										 const PointType& b,
										 InternalFloatType eps) noexcept {
			for (int k = 0; k < 3; ++k) {
				InternalFloatType d = std::abs(a[k] - b[k]);
				d -= std::floor(d);
				d = std::min(d, InternalFloatType(1) - d);
				if (d > eps) return false;
			}
			return true;
		}
	};


	// ===========================================================================
	//  VoronoiFused — implementation
	// ===========================================================================

	inline VoronoiFused VoronoiFused::build(const VoronoiContext& ctx) {
		VoronoiFused out;

		const size_t n_cells = ctx.n_cells();
		out.report.cells = n_cells;

		const MatrixType& CtoF = ctx.grid.cart_to_frac();

		// Per-cell local -> global index maps, one entry per local element.
		std::vector<std::vector<uint32_t>> cell_vert_gid;
		std::vector<std::vector<uint32_t>> cell_edge_gid;
		std::vector<std::vector<uint32_t>> cell_face_gid;

		// Pipeline: five independent steps, each one self-contained.
		merge_vertices(ctx, CtoF, out.vertices, cell_vert_gid);
		merge_edges(ctx, cell_vert_gid, out.edges, cell_edge_gid);
		merge_faces(ctx, CtoF, cell_vert_gid, cell_edge_gid, out, cell_face_gid);
		assemble_polyhedra(ctx, cell_vert_gid, cell_edge_gid, cell_face_gid, out);
		compute_volumes(out);

		return out;
	}

	// ---------------------------------------------------------------------------
	//  Step 1 — merge vertices
	//
	//  Two-phase approach:
	//
	//    Phase A (sliding window): sort all vertex entries by the wrapped
	//    fractional coordinate sum, then collapse adjacent entries within a
	//    tolerance. This assigns an initial global vertex ID to every entry.
	//
	//    Phase B (boundary wrap-up): entries near key ~ 0 and entries near
	//    key ~ 3 may describe the same physical vertex across a periodic
	//    boundary, but the sliding window misses those pairs. Instead of
	//    merging them by erasing out_vertices (O(N) per match, quadratic in
	//    the worst case), we collect the matches and resolve them with a
	//    union-find over the global IDs. The final compaction of
	//    out_vertices is a single O(N) pass.
	//
	//  Complexity: O(N log N) for the sort + O(N * alpha(N)) for the
	//  union-find, versus O(N * K) for the erase-based version, where K is
	//  the number of boundary matches.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::merge_vertices(
		const VoronoiContext& ctx,
		const MatrixType& CtoF,
		std::vector<PointType>& out_vertices,
		std::vector<std::vector<uint32_t>>& cell_vert_gid)
	{
		// Tolerance for matching vertices across cells, in fractional units.
		// 1e-8 is ~1e-6 A for a 100 A cell: large enough to absorb numerical
		// drift, small enough that distinct physical vertices never collide.
		const InternalFloatType merge_eps = InternalFloatType(1e-8);

		struct VEntry {
			AtomId            cell;
			LocalIdx          local;
			PointType         pos_cart;
			PointType         pos_frac;   // wrapped to [0, 1)
			InternalFloatType key;        // pos_frac[0] + pos_frac[1] + pos_frac[2]
		};

		const size_t n_cells = ctx.n_cells();

		// ---- Collect all entries ------------------------------------------
		size_t total_v = 0;
		for (const auto& c : ctx.cells) total_v += c.topo.vert.size();

		std::vector<VEntry> v;
		v.reserve(total_v);
		cell_vert_gid.assign(n_cells, {});

		for (AtomId ci = 0; ci < n_cells; ++ci) {
			const auto& topo = ctx.cells[ci].topo;
			const Vertex* vd = topo.vert.data();
			const size_t  nv = topo.vert.size();
			cell_vert_gid[ci].assign(nv, std::numeric_limits<uint32_t>::max());

			for (size_t li = 0; li < nv; ++li) {
				const PointType& pc = vd[li].pos;
				PointType pf = CtoF * pc;
				pf[0] -= std::floor(pf[0]);
				pf[1] -= std::floor(pf[1]);
				pf[2] -= std::floor(pf[2]);

				v.push_back({ci, static_cast<LocalIdx>(li), pc, pf,
							  pf[0] + pf[1] + pf[2]});
			}
		}

		std::sort(v.begin(), v.end(),
				  [](const VEntry& a, const VEntry& b) { return a.key < b.key; });

		// ---- Phase A: sliding window to assign initial global IDs ---------
		constexpr uint32_t NO_ID = std::numeric_limits<uint32_t>::max();

		std::vector<uint32_t> assigned(v.size(), NO_ID);
		uint32_t num_globals = 0;

		for (size_t i = 0; i < v.size(); ++i) {
			if (assigned[i] != NO_ID) continue;

			const uint32_t gid = num_globals++;
			out_vertices.push_back(v[i].pos_cart);
			assigned[i] = gid;

			for (size_t j = i + 1; j < v.size(); ++j) {
				if (v[j].key - v[i].key > 3 * merge_eps) break;
				if (assigned[j] != NO_ID) continue;
				if (same_frac_pbc(v[j].pos_frac, v[i].pos_frac, merge_eps))
					assigned[j] = gid;
			}
		}

		// ---- Phase B: boundary wrap-up via union-find --------------------
		//
		//  The sliding window cannot see pairs that straddle the wrap
		//  boundary: entries near key ~ 0 and entries near key ~ 3. Scan
		//  those two short ranges and record every match as a union.
		//
		//  Nothing is erased or shifted during the scan; the union-find
		//  resolves all chains of merges at once, and a single compaction
		//  pass rebuilds out_vertices.
		{
			// Left edge: entries with key close to 0.
			size_t lo_end = 0;
			while (lo_end < v.size() && v[lo_end].key < 3 * merge_eps) ++lo_end;

			// Right edge: entries with key close to 3.
			size_t hi_beg = v.size();
			while (hi_beg > 0 && v[hi_beg - 1].key > 3 - 3 * merge_eps) --hi_beg;

			// Union-find over global IDs. Path halving for cheap lookup.
			std::vector<uint32_t> parent(num_globals);
			std::iota(parent.begin(), parent.end(), 0u);

			auto find_root = [&parent](uint32_t x) noexcept {
				while (parent[x] != x) {
					parent[x] = parent[parent[x]];   // path halving
					x = parent[x];
				}
				return x;
				};

			for (size_t i = 0; i < lo_end; ++i) {
				for (size_t j = hi_beg; j < v.size(); ++j) {
					const uint32_t gi = assigned[i];
					const uint32_t gj = assigned[j];
					if (gi == NO_ID || gj == NO_ID || gi == gj) continue;

					if (same_frac_pbc(v[i].pos_frac, v[j].pos_frac, merge_eps)) {
						const uint32_t ri = find_root(gi);
						const uint32_t rj = find_root(gj);
						if (ri != rj) {
							// Always attach the larger root to the smaller one.
							// This keeps the smallest global ID as the
							// representative, which preserves insertion order
							// in the compacted array.
							const uint32_t lo = std::min(ri, rj);
							const uint32_t hi = std::max(ri, rj);
							parent[hi] = lo;
						}
					}
				}
			}

			// ---- Compact out_vertices in a single pass -------------------
			//
			//  Build new_id[old_global] -> new_global. Only roots survive;
			//  their positions in out_vertices become contiguous indices in
			//  the compacted array.
			std::vector<uint32_t> new_id(num_globals, NO_ID);
			uint32_t next_id = 0;
			for (uint32_t gid = 0; gid < num_globals; ++gid) {
				if (find_root(gid) == gid)
					new_id[gid] = next_id++;
			}
			for (uint32_t gid = 0; gid < num_globals; ++gid) {
				if (new_id[gid] == NO_ID)
					new_id[gid] = new_id[find_root(gid)];
			}

			// Move surviving vertices into a fresh buffer. Non-roots are
			// dropped; their position is not needed by anyone.
			std::vector<PointType> compacted(next_id);
			for (uint32_t gid = 0; gid < num_globals; ++gid) {
				if (find_root(gid) == gid)
					compacted[new_id[gid]] = out_vertices[gid];
			}
			out_vertices = std::move(compacted);

			// Remap every assigned entry through new_id.
			for (auto& a : assigned) {
				if (a != NO_ID) a = new_id[a];
			}
		}

		// ---- Store global IDs into per-cell maps --------------------------
		for (size_t e = 0; e < v.size(); ++e)
			cell_vert_gid[v[e].cell][v[e].local] = assigned[e];
	}

	// ---------------------------------------------------------------------------
	//  Step 2 — merge edges
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::merge_edges(
		const VoronoiContext& ctx,
		const std::vector<std::vector<uint32_t>>& cell_vert_gid,
		std::vector<EdgeFused>& out_edges,
		std::vector<std::vector<uint32_t>>& cell_edge_gid)
	{
		struct EEntry {
			uint32_t v0, v1;
			AtomId   cell;
			LocalIdx local;
		};

		const size_t n_cells = ctx.n_cells();

		size_t total_e = 0;
		for (const auto& c : ctx.cells) total_e += c.topo.edge.size();

		std::vector<EEntry> e;
		e.reserve(total_e);
		cell_edge_gid.assign(n_cells, {});

		for (AtomId ci = 0; ci < n_cells; ++ci) {
			const auto& topo = ctx.cells[ci].topo;
			const HalfEdge* ed = topo.edge.data();
			const size_t    ne = topo.edge.size();
			cell_edge_gid[ci].assign(ne, std::numeric_limits<uint32_t>::max());

			for (size_t li = 0; li < ne; ++li) {
				const HalfEdge& he = ed[li];
				const HalfEdge& nxt = topo.edge.get_by_id(he.next_edge_id);

				const LocalIdx lv0 = topo.vert.index_of(he.origin_vertex_id);
				const LocalIdx lv1 = topo.vert.index_of(nxt.origin_vertex_id);

				const uint32_t gv0 = cell_vert_gid[ci][lv0];
				const uint32_t gv1 = cell_vert_gid[ci][lv1];

				e.push_back({std::min(gv0, gv1), std::max(gv0, gv1),
							  ci, static_cast<LocalIdx>(li)});
			}
		}

		std::sort(e.begin(), e.end(), [](const EEntry& a, const EEntry& b) {
			if (a.v0 != b.v0) return a.v0 < b.v0;
			return a.v1 < b.v1;
				  });

		for (size_t i = 0; i < e.size(); ) {
			const uint32_t gid = static_cast<uint32_t>(out_edges.size());
			out_edges.push_back(EdgeFused{e[i].v0, e[i].v1});
			size_t j = i;
			while (j < e.size() && e[j].v0 == e[i].v0 && e[j].v1 == e[i].v1) {
				cell_edge_gid[e[j].cell][e[j].local] = gid;
				++j;
			}
			i = j;
		}
	}

	// ---------------------------------------------------------------------------
	//  Step 3 — merge faces
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::merge_faces(
		const VoronoiContext& ctx,
		const MatrixType& CtoF,
		const std::vector<std::vector<uint32_t>>& cell_vert_gid,
		const std::vector<std::vector<uint32_t>>& cell_edge_gid,
		VoronoiFused& out,
		std::vector<std::vector<uint32_t>>& cell_face_gid)
	{
		std::vector<FaceRec> recs = collect_face_records(ctx, cell_face_gid);

		std::sort(recs.begin(), recs.end(),
				  [](const FaceRec& x, const FaceRec& y) {
					  if (x.a != y.a) return x.a < y.a;
					  if (x.b != y.b) return x.b < y.b;
					  return x.s < y.s;
				  });

		out.report.face_records = recs.size();

		size_t i = 0;
		while (i < recs.size()) {
			const size_t j = find_group_end(recs, i);
			emit_face_group(ctx, CtoF, cell_vert_gid, cell_edge_gid,
							out, cell_face_gid, recs, i, j);
			i = j;
		}
	}

	// ---------------------------------------------------------------------------
	//  Step 3a — collect one FaceRec per (cell, face) pair.
	// ---------------------------------------------------------------------------
	inline std::vector<VoronoiFused::FaceRec> VoronoiFused::collect_face_records(
		const VoronoiContext& ctx,
		std::vector<std::vector<uint32_t>>& cell_face_gid)
	{
		const size_t n_cells = ctx.n_cells();

		size_t total_f = 0;
		for (const auto& c : ctx.cells) total_f += c.topo.face.size();

		std::vector<FaceRec> recs;
		recs.reserve(total_f);
		cell_face_gid.assign(n_cells, {});

		for (AtomId ci = 0; ci < n_cells; ++ci) {
			const auto& topo = ctx.cells[ci].topo;
			const Polygon* fd = topo.face.data();
			const size_t   nf = topo.face.size();
			cell_face_gid[ci].assign(nf, std::numeric_limits<uint32_t>::max());

			for (size_t li = 0; li < nf; ++li) {
				const AtomId  oa = fd[li].other_cell;
				const AtomId  ow = ci;
				const uint8_t sc = fd[li].other_shift.get_code();

				// Canonical orientation: (a, b, s) with a < b and s the
				// shift from a to b. If the record already has owner <
				// other, keep it as-is; otherwise invert both.
				AtomId  ca, cb;
				uint8_t cs;
				if (ow <= oa) {
					ca = ow; cb = oa; cs = sc;
				} else {
					ca = oa; cb = ow;
					cs = ShiftCode::inverse(static_cast<uint8_t>(sc));
				}
				recs.push_back({ca, cb, cs, ci, static_cast<LocalIdx>(li)});
			}
		}
		return recs;
	}

	// ---------------------------------------------------------------------------
	//  Step 3b — group scan: end index of the run with the same (a, b, s).
	// ---------------------------------------------------------------------------
	inline size_t VoronoiFused::find_group_end(const std::vector<FaceRec>& recs,
											   size_t i)
	{
		size_t j = i + 1;
		while (j < recs.size() &&
			   recs[j].a == recs[i].a &&
			   recs[j].b == recs[i].b &&
			   recs[j].s == recs[i].s) ++j;
		return j;
	}

	// ---------------------------------------------------------------------------
	//  Step 3c — emit one canonical polygon for a group.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::emit_face_group(
		const VoronoiContext& ctx,
		const MatrixType& CtoF,
		const std::vector<std::vector<uint32_t>>& cell_vert_gid,
		const std::vector<std::vector<uint32_t>>& cell_edge_gid,
		VoronoiFused& out,
		std::vector<std::vector<uint32_t>>& cell_face_gid,
		const std::vector<FaceRec>& recs,
		size_t i, size_t j)
	{
		const size_t group_size = j - i;

		// (1) Merge vertex sets of every record in the group.
		std::vector<MergedVert> merged;
		size_t min_single_count = 0, max_single_count = 0;

		merge_group_vertices(ctx, CtoF, cell_vert_gid, out.vertices,
							 recs, i, j, merged,
							 min_single_count, max_single_count);

		if (group_size > 2) {
			out.report.vertex_set_bad += group_size - 2;
			out.report.messages.emplace_back(
				"duplicate records for face ("
				+ std::to_string(recs[i].a) + ","
				+ std::to_string(recs[i].b) + ",s="
				+ std::to_string(recs[i].s) + ")");
		}

		if (merged.size() < 3) {
			++out.report.vertex_set_bad;
			return;
		}

		// (2) Order the merged set by polar angle around the centroid.
		order_by_polar_angle(merged);

		// (3) Emit the canonical polygon.
		const uint32_t gp = static_cast<uint32_t>(out.polygons.size());
		out.polygons.push_back(PolygonFused{});
		PolygonFused& pg = out.polygons.back();

		// Vertices — cyclic order.
		pg.vert_offset = static_cast<uint32_t>(out.poly_verts.size());
		pg.vert_count = static_cast<uint16_t>(merged.size());
		for (const auto& m : merged)
			out.poly_verts.push_back(m.gid);

		// Edges — taken from the most complete source record.
		emit_polygon_edges(ctx, cell_edge_gid, recs, i, j, pg, out.poly_edges);

		// Owner / other / shift — canonical metadata.
		emit_polygon_metadata(ctx, recs, i, pg);

		// Metrics.
		pg.area = compute_face_area(merged);
		pg.solid_angle = compute_face_solid_angle(merged,
												  ctx.cells[pg.owner_atom].meta.center);

											  // (4) Register the polygon with every source record.
		for (size_t k = i; k < j; ++k)
			cell_face_gid[recs[k].src_cell][recs[k].src_local] = gp;

		// (5) Validation.
		validate_face_group(ctx, recs, i, j, merged,
							min_single_count, max_single_count, gp, out);

		++out.report.physical_faces;
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.i — collect merged vertex set across the group.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::merge_group_vertices(
		const VoronoiContext& ctx,
		const MatrixType& CtoF,
		const std::vector<std::vector<uint32_t>>& cell_vert_gid,
		const std::vector<PointType>& global_vertices,
		const std::vector<FaceRec>& recs,
		size_t i, size_t j,
		std::vector<MergedVert>& merged,
		size_t& min_single_count,
		size_t& max_single_count)
	{
		const InternalFloatType merge_eps = InternalFloatType(1e-8);

		merged.clear();
		max_single_count = 0;
		min_single_count = std::numeric_limits<size_t>::max();

		for (size_t k = i; k < j; ++k) {
			const auto& topo = ctx.cells[recs[k].src_cell].topo;
			const Polygon& poly = topo.face.data()[recs[k].src_local];

			size_t cnt = 0;
			LocalId cur = poly.first_edge_id;
			do {
				const HalfEdge& he = topo.edge.get_by_id(cur);
				const LocalIdx  lv = topo.vert.index_of(he.origin_vertex_id);
				const uint32_t  gid = cell_vert_gid[recs[k].src_cell][lv];

				const PointType& pc = global_vertices[gid];
				PointType pf = CtoF * pc;
				pf[0] -= std::floor(pf[0]);
				pf[1] -= std::floor(pf[1]);
				pf[2] -= std::floor(pf[2]);

				// Deduplicate against already-collected vertices.
				bool dup = false;
				for (const auto& m : merged) {
					if (same_frac_pbc(m.pos_frac, pf, merge_eps)) {
						dup = true;
						break;
					}
				}
				if (!dup) merged.push_back({pc, pf, gid});

				++cnt;
				cur = he.next_edge_id;
			} while (cur != poly.first_edge_id);

			max_single_count = std::max(max_single_count, cnt);
			min_single_count = std::min(min_single_count, cnt);
		}
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.ii — polar sort of the merged vertex set.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::order_by_polar_angle(std::vector<MergedVert>& merged) {
		// Centroid — origin for polar sorting.
		PointType centroid(0, 0, 0);
		for (const auto& m : merged) centroid = centroid + m.pos_cart;
		centroid = centroid / static_cast<InternalFloatType>(merged.size());

		// Face normal via Newell's method over the unordered set.
		InternalFloatType nx = 0, ny = 0, nz = 0;
		for (size_t k = 0; k < merged.size(); ++k) {
			const PointType& a = merged[k].pos_cart;
			const PointType& b = merged[(k + 1) % merged.size()].pos_cart;
			nx += (a[1] - b[1]) * (a[2] + b[2]);
			ny += (a[2] - b[2]) * (a[0] + b[0]);
			nz += (a[0] - b[0]) * (a[1] + b[1]);
		}
		PointType nrm(nx, ny, nz);
		const InternalFloatType nlen = nrm.r();
		if (nlen > 1e-12) nrm = nrm / nlen;

		// In-plane orthonormal basis.
		PointType up(0, 0, 1);
		if (std::abs(nrm[2]) > InternalFloatType(0.9)) up = PointType(1, 0, 0);
		PointType ex = PointType::Vector(up, nrm);
		const InternalFloatType exlen = ex.r();
		if (exlen > 1e-12) ex = ex / exlen;
		const PointType ey = PointType::Vector(nrm, ex);

		std::sort(merged.begin(), merged.end(),
				  [&](const MergedVert& A, const MergedVert& B) {
					  const PointType da = A.pos_cart - centroid;
					  const PointType db = B.pos_cart - centroid;
					  const InternalFloatType aa = std::atan2(
						  PointType::Scalar(da, ey), PointType::Scalar(da, ex));
					  const InternalFloatType ab = std::atan2(
						  PointType::Scalar(db, ey), PointType::Scalar(db, ex));
					  return aa < ab;
				  });
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.iii — emit edges for a face group.
	//
	//  Take them from the record with the largest vertex count. If drift
	//  dropped one vertex in the twin record, the larger record still
	//  carries every edge of the merged face.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::emit_polygon_edges(
		const VoronoiContext& ctx,
		const std::vector<std::vector<uint32_t>>& cell_edge_gid,
		const std::vector<FaceRec>& recs,
		size_t i, size_t j,
		PolygonFused& pg,
		std::vector<uint32_t>& poly_edges)
	{
		// Pick the record with the most edges.
		size_t best_k = i;
		size_t best_cnt = 0;
		for (size_t k = i; k < j; ++k) {
			const auto& topo = ctx.cells[recs[k].src_cell].topo;
			const Polygon& poly = topo.face.data()[recs[k].src_local];
			size_t cnt = 0;
			LocalId cur = poly.first_edge_id;
			do {
				++cnt;
				cur = topo.edge.get_by_id(cur).next_edge_id;
			} while (cur != poly.first_edge_id);
			if (cnt > best_cnt) {
				best_cnt = cnt; best_k = k;
			}
		}

		// Emit that record's edges as global IDs.
		const auto& src_topo = ctx.cells[recs[best_k].src_cell].topo;
		const Polygon& src_poly = src_topo.face.data()[recs[best_k].src_local];

		pg.edge_offset = static_cast<uint32_t>(poly_edges.size());
		pg.edge_count = 0;

		LocalId cur = src_poly.first_edge_id;
		do {
			const HalfEdge& he = src_topo.edge.get_by_id(cur);
			const LocalIdx  le = src_topo.edge.index_of(cur);
			poly_edges.push_back(cell_edge_gid[recs[best_k].src_cell][le]);
			++pg.edge_count;
			cur = he.next_edge_id;
		} while (cur != src_poly.first_edge_id);
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.iv — emit canonical owner / other / shift metadata.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::emit_polygon_metadata(
		const VoronoiContext& ctx,
		const std::vector<FaceRec>& recs,
		size_t i,
		PolygonFused& pg)
	{
		const auto& src_topo = ctx.cells[recs[i].src_cell].topo;
		const Polygon& src_poly = src_topo.face.data()[recs[i].src_local];

		AtomId    owner = recs[i].src_cell;
		AtomId    other = src_poly.other_cell;
		ShiftCode shift = src_poly.other_shift;

		// Canonical orientation: owner < other, shift points owner -> other.
		if (owner > other) {
			std::swap(owner, other);
			shift = ShiftCode::inverse(shift.get_code());
		}
		pg.owner_atom = owner;
		pg.other_atom = other;
		pg.other_shift = shift;
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.v — face metrics.
	// ---------------------------------------------------------------------------
	inline InternalFloatType VoronoiFused::compute_face_area(
		const std::vector<MergedVert>& merged)
	{
		if (merged.size() < 3) return InternalFloatType(0);

		InternalFloatType area = 0;
		const PointType p0 = merged[0].pos_cart;
		PointType v1 = merged[1].pos_cart - p0;
		for (size_t k = 2; k < merged.size(); ++k) {
			const PointType v2 = merged[k].pos_cart - p0;
			area += MathEngine::triangle_area(v1, v2);
			v1 = v2;
		}
		return area;
	}

	inline InternalFloatType VoronoiFused::compute_face_solid_angle(
		const std::vector<MergedVert>& merged,
		const PointType& origin)
	{
		if (merged.size() < 3) return InternalFloatType(0);

		const PointType v0 = merged[0].pos_cart - origin;
		const InternalFloatType r0 = v0.r();
		PointType         v_prev = merged[1].pos_cart - origin;
		InternalFloatType r_prev = v_prev.r();
		InternalFloatType dot0_prev = PointType::Scalar(v0, v_prev);

		InternalFloatType solid = 0;
		for (size_t k = 2; k < merged.size(); ++k) {
			const PointType v_curr = merged[k].pos_cart - origin;
			const InternalFloatType r_curr = v_curr.r();
			const InternalFloatType dot0_curr = PointType::Scalar(v0, v_curr);

			const InternalFloatType cross = std::abs(
				PointType::Scalar(v0, PointType::Vector(v_prev, v_curr)));
			const InternalFloatType denom =
				r0 * r_prev * r_curr
				+ dot0_prev * r_curr
				+ dot0_curr * r_prev
				+ PointType::Scalar(v_prev, v_curr) * r0;
			if (cross > MathEngine::EPS)
				solid += std::atan2(cross, denom);

			v_prev = v_curr;
			r_prev = r_curr;
			dot0_prev = dot0_curr;
		}
		return solid * InternalFloatType(2);
	}

	// ---------------------------------------------------------------------------
	//  Step 3c.vi — group validation. Records anomalies without aborting.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::validate_face_group(
		const VoronoiContext& ctx,
		const std::vector<FaceRec>& recs,
		size_t i, size_t j,
		const std::vector<MergedVert>& merged,
		size_t min_single_count,
		size_t max_single_count,
		uint32_t gp,
		VoronoiFused& out)
	{
		const size_t group_size = j - i;

		// Per-record vertex count drift.
		if (group_size > 1 && max_single_count != min_single_count) {
			++out.report.vertex_count_bad;
			out.report.messages.emplace_back(
				"vertex count drift on face ("
				+ std::to_string(recs[i].a) + ","
				+ std::to_string(recs[i].b) + ",s="
				+ std::to_string(recs[i].s) + "): records ["
				+ std::to_string(min_single_count) + ".."
				+ std::to_string(max_single_count) + "], merged "
				+ std::to_string(merged.size()));
		}

		if (group_size == 1) {
			const auto& src_topo = ctx.cells[recs[i].src_cell].topo;
			const Polygon& src_poly = src_topo.face.data()[recs[i].src_local];
			if (src_poly.other_cell == recs[i].src_cell
				&& src_poly.other_shift.get_code() != 13)
				++out.report.self_faces;
			else
				++out.report.twin_missing;
			return;
		}

		if (group_size == 2) {
			++out.report.twin_pairs;

			const auto& A = recs[i];
			const auto& B = recs[i + 1];

			const bool atom_ok =
				A.src_cell != B.src_cell &&
				ctx.cells[A.src_cell].topo.face.data()[A.src_local].other_cell
				== B.src_cell &&
				ctx.cells[B.src_cell].topo.face.data()[B.src_local].other_cell
				== A.src_cell;
			if (!atom_ok) {
				++out.report.atom_pair_bad;
				out.report.messages.emplace_back(
					"atom pair mismatch on face " + std::to_string(gp));
			}
			const uint8_t inv_a = ShiftCode::inverse(
				ctx.cells[A.src_cell].topo.face.data()[A.src_local]
				.other_shift.get_code());
			const uint8_t sb =
				ctx.cells[B.src_cell].topo.face.data()[B.src_local]
				.other_shift.get_code();
			if (sb != inv_a) {
				++out.report.shift_bad;
				out.report.messages.emplace_back(
					"shift mismatch on face " + std::to_string(gp));
			}
		}
	}

	// ---------------------------------------------------------------------------
	//  Step 4 — assemble per-cell polyhedra (CSR).
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::assemble_polyhedra(
		const VoronoiContext& ctx,
		const std::vector<std::vector<uint32_t>>& cell_vert_gid,
		const std::vector<std::vector<uint32_t>>& cell_edge_gid,
		const std::vector<std::vector<uint32_t>>& cell_face_gid,
		VoronoiFused& out)
	{
		const size_t n_cells = ctx.n_cells();
		out.polyhedra.resize(n_cells);

		// Reserve the flat buffers up front.
		size_t tot_v = 0, tot_e = 0, tot_p = 0;
		for (AtomId ci = 0; ci < n_cells; ++ci) {
			tot_v += ctx.cells[ci].topo.vert.size();
			tot_e += ctx.cells[ci].topo.edge.size();
			tot_p += ctx.cells[ci].topo.face.size();
		}
		out.ph_verts.reserve(tot_v);
		out.ph_edges.reserve(tot_e);
		out.ph_polys.reserve(tot_p);

		for (AtomId ci = 0; ci < n_cells; ++ci) {
			const auto& src = ctx.cells[ci];
			auto& dst = out.polyhedra[ci];
			dst.center = src.meta.center;

			dst.vert_offset = static_cast<uint32_t>(out.ph_verts.size());
			dst.vert_count = static_cast<uint16_t>(src.topo.vert.size());
			for (size_t li = 0; li < src.topo.vert.size(); ++li)
				out.ph_verts.push_back(cell_vert_gid[ci][li]);

			dst.edge_offset = static_cast<uint32_t>(out.ph_edges.size());
			dst.edge_count = static_cast<uint16_t>(src.topo.edge.size());
			for (size_t li = 0; li < src.topo.edge.size(); ++li)
				out.ph_edges.push_back(cell_edge_gid[ci][li]);

			dst.poly_offset = static_cast<uint32_t>(out.ph_polys.size());
			dst.poly_count = static_cast<uint16_t>(src.topo.face.size());
			for (size_t li = 0; li < src.topo.face.size(); ++li)
				out.ph_polys.push_back(cell_face_gid[ci][li]);
		}
	}

	// ---------------------------------------------------------------------------
	//  Step 5 — volumes via the pyramid formula.
	// ---------------------------------------------------------------------------
	inline void VoronoiFused::compute_volumes(VoronoiFused& out) {
		for (auto& ph : out.polyhedra) {
			InternalFloatType vol = 0;

			for (uint16_t k = 0; k < ph.poly_count; ++k) {
				const uint32_t pid = out.ph_polys[ph.poly_offset + k];
				const PolygonFused& poly = out.polygons[pid];
				if (poly.vert_count < 3) continue;

				const PointType& p0 = out.vertices[out.poly_verts[poly.vert_offset + 0]];
				const PointType& p1 = out.vertices[out.poly_verts[poly.vert_offset + 1]];
				const PointType& p2 = out.vertices[out.poly_verts[poly.vert_offset + 2]];

				geometry::Plane<InternalFloatType> plane(p0, p1, p2);
				const InternalFloatType h = std::abs(plane.distance(ph.center));
				vol += poly.area * h;
			}
			ph.volume = vol / InternalFloatType(3);
		}
	}

	// ---------------------------------------------------------------------------
	//  VoronoiPipeline::fuse — convenience wrapper.
	// ---------------------------------------------------------------------------
	inline VoronoiFused VoronoiPipeline::fuse() const {
		return VoronoiFused::build(context);
	}

} // namespace cpplib::voronoi