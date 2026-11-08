// ============================================================================
//  Surface Construction Types and Constants Documentation
// ============================================================================
//
//  This header-only file documents the type system and constants used by
//  the Bader surface construction module. Include this along with
//  SurfaceConstruction.h for reference.
// ============================================================================

#pragma once

#include "SurfaceConstruction.h"

/**
 * @file SurfaceConstruction_types.h
 * @brief Type documentation for Bader surface construction module
 *
 * This module implements Stages A-D of the Bader zero-flux surface construction:
 * - Stage A: Define point density constant and data structures
 * - Stage B: Initialize surfaces from Voronoi faces using Poisson disk sampling
 * - Stage C: rtsafe iteration for zero-flux boundary finding
 * - Stage D: Build boundaries by connecting surface points
 */

namespace cpplib {

/**
 * @brief Constants for surface construction configuration
 *
 * These constants control the density and quality of the initial surface points.
 * They can be adjusted to balance memory usage vs. accuracy.
 */
namespace surface_constants {

    /**
     * Point density on surface (Å² per point).
     *
     * @note This is a compile-time constant for fast parameter tuning in tests.
     */
    constexpr double POINT_DENSITY_A2 = 0.1;

    /**
     * Minimum number of points per polygon.
     * Ensures adequate sampling even for very small faces.
     */
    static constexpr uint32_t MIN_POINTS_PER_POLYGON = 4;

    /**
     * Maximum iterations for Poisson disk sampling.
     * If target count not reached, returns whatever was generated.
     */
    static constexpr uint32_t MAX_POISSON_ITERATIONS = 1024;

}  // namespace surface_constants


/**
 * @brief Surface point structure used throughout the pipeline
 *
 * Each SurfacePoint represents a seed location on the boundary between two atoms.
 * Points are initialized from Voronoi shared faces and refined via rtsafe iteration.
 */
struct SurfacePoint {
    PointType pos;           ///< Position on the surface (Cartesian coordinates)
    PointType normal;        ///< Outward normal (unit vector), updated by Hessian analysis
    uint32_t basin_a;        ///< Owner atom index
    uint32_t basin_b;        ///< Other atom index (shared face partner)

    /**
     * Create initial surface point with given position and normal.
     *
     * @param r Position on the surface
     * @param normal_in Incoming normal direction (will be normalized)
     * @param b_a Basin A identifier (owner atom index)
     * @param b_b Basin B identifier (other atom index - shared)
     * @return Initialized SurfacePoint
     */
    static SurfacePoint make_initial(const PointType& r,
                                    const PointType& normal_in,
                                    uint32_t b_a, uint32_t b_b);
};


/**
 * @brief Polygon area computation utilities
 */
struct PolygonAreaCalculator {

    /**
     * Compute target point count for uniform distribution.
     *
     * Formula: target = ceil(area / POINT_DENSITY_A2)
     * Minimum of MIN_POINTS_PER_POLYGON is enforced.
     *
     * @param area Polygon area in Å²
     * @param density Target density (points per Å²), defaults to 1/POINT_DENSITY_A2
     * @return Target point count
     */
    static uint32_t compute_target_points(InternalFloatType area,
                                          double density = 1.0 / surface_constants::POINT_DENSITY_A2);

    /**
     * Compute target point count with MIN_POINTS_PER_POLYGON enforcement.
     * Public wrapper that ensures minimum sampling for small polygons.
     */
    static uint32_t compute_target_points_enforced(InternalFloatType area,
                                                   double density = 1.0 / surface_constants::POINT_DENSITY_A2) {
        // FIX: Enforce minimum point count even for tiny polygons
        const auto target = std::ceil(area / density);
        return std::max<uint32_t>(target, surface_constants::MIN_POINTS_PER_POLYGON);
    }

    /**
     * Compute approximate spacing between points.
     *
     * Formula: spacing ≈ sqrt(area / n_points)
     *
     * @param area Polygon area in Å²
     * @param n_points Number of points
     * @return Approximate average spacing in Å
     */
    static InternalFloatType compute_spacing(InternalFloatType area, uint32_t n_points);

};


/**
 * @brief Surface initialization from Voronoi data
 */
struct SurfaceInitialization {

    /**
     * Create initial surface points from Voronoi shared faces.
     *
     * For each shared face between two atoms:
     * 1. Uses precomputed area from Voronoi pipeline
     * 2. Computes target point count based on POINT_DENSITY_A2
     * 3. Creates seed points at centroid with Voronoi normal
     * 4. Produces TWO surface points per face (one for each basin side)
     *
     * @param vf Fused Voronoi data with precomputed polygon areas
     * @param psoa Spatial grid for electron density evaluation
     * @param ed ElectronDensity accessor (defaults to singleton)
     * @return Vector of initial SurfacePoints ready for rtsafe iteration
     */
    static std::vector<SurfacePoint> from_voronoi(
        const voronoi::VoronoiFused& vf,
        const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
        const ElectronDensity& ed = ElectronDensity::instance());

};


/**
 * @brief Poisson disk sampler for uniform surface coverage
 *
 * Generates points with minimum separation while attempting to reach target count.
 */
struct PoissonDiskSampler {

    /**
     * Generate uniformly distributed points on polygon.
     *
     * Algorithm: Stochastic rejection sampling with Gaussian spread around
     * accepted points, minimum separation enforcement via distance buffer.
     *
     * @param verts Ordered polygon vertices (CCW from outside view)
     * @param target_count Desired number of points
     * @param min_separation Minimum distance between any two points (Å)
     * @return Sampled point positions (may be fewer than target if max iterations reached)
     */
    static std::vector<PointType> sample(
        const std::vector<PointType>& verts,
        uint32_t target_count,
        double min_separation = 0.5);

};

}  // namespace cpplib


/**
 * @brief Quick reference for common operations
 */

// ============================================================================
//  RTSafe ITERATION (Stage C) - Zero-flux surface refinement
// ============================================================================

namespace cpplib {

/**
 * Refine surface points using r-squared safe iteration.
 *
 * Uses gradient and Hessian from CriticalPoint to iteratively move points toward
 * zero-flux condition where ∇ρ is perpendicular to the surface normal.
 *
 * @param sp Points to refine (modified in place)
 * @param psoa Spatial grid for electron density evaluation
 * @param max_iterations Maximum iterations per point
 * @param tolerance Convergence threshold for |∇ρ ⋅ tangent|
 * @return Number of successfully converged points
 */
inline std::vector<SurfacePoint> rtsafe_iteration(
    std::vector<SurfacePoint>& sp,
    const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
    uint32_t max_iterations = 50,
    double tolerance = 1e-8,
    std::vector<SurfacePointEigenAnalysis>* return_analysis = nullptr
) {
    return SurfaceConstruction::rtsafe_iteration_with_lambda_criteria(
        sp, psoa, max_iterations, tolerance, return_analysis);
}

}  // namespace cpplib

// ============================================================================
//  EXPORTED FUNCTIONS (convenience wrappers)
// ============================================================================

/**
 * Get polygon area from precomputed Voronoi data.
 *
 * @param poly Polygon from VoronoiFused
 * @return Area in Å²
 */
inline InternalFloatType get_polygon_area(const voronoi::PolygonFused& poly) {
    return poly.area;
}

/**
 * Compute target point count for polygon area.
 *
 * @param area Polygon area in Å²
 * @return Number of points needed
 */
inline uint32_t compute_surface_points_for_area(InternalFloatType area,
                                                double density = 1.0 / surface_constants::POINT_DENSITY_A2) {
    return PolygonAreaCalculator::compute_target_points(area, density);
}

/**
 * Initialize surface from Voronoi decomposition.
 *
 * Main entry point for Bader boundary construction starting from Voronoi faces.
 *
 * @param vf Fused Voronoi data
 * @param psoa Spatial grid for density evaluation
 * @return Initial surface points for rtsafe/boundary pipeline
 */
inline std::vector<SurfacePoint> initialize_surface_from_voronoi(
    const voronoi::VoronoiFused& vf,
    const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
    const ElectronDensity& ed = ElectronDensity::instance()) {
    return SurfaceInitialization::from_voronoi(vf, psoa, ed);
}

/**
 * Generate uniformly distributed points using Poisson disk sampling.
 *
 * @param verts Polygon vertices (CCW ordered)
 * @param n_points Target count
 * @param min_separation Minimum distance (Å), defaults to 0.5 Å
 * @return Sampled points
 */
inline std::vector<PointType> poisson_disk_sample(
    const std::vector<PointType>& verts,
    uint32_t n_points,
    double min_separation = 0.5) {
    return PoissonDiskSampler::sample(verts, n_points, min_separation);
}

/**
 * Compute target point count with MIN_POINTS_PER_POLYGON enforcement.
 */
inline uint32_t compute_surface_points_for_area_enforced(InternalFloatType area,
                                                         double density = 1.0 / surface_constants::POINT_DENSITY_A2) {
    const auto target = PolygonAreaCalculator::compute_target_points(area, density);
    return std::max(target, surface_constants::MIN_POINTS_PER_POLYGON);
}

// ============================================================================
//  RTSafe ITERATION (Stage C) - Zero-flux surface refinement
// ============================================================================

/**
 * Refine surface points using r-squared safe iteration.
 *
 * Uses gradient and Hessian from CriticalPoint to iteratively move points toward
 * zero-flux condition where ∇ρ is perpendicular to surface normal.
 *
 * @param sp Points to refine (modified in place)
 * @param psoa Spatial grid for electron density evaluation
 * @param max_iterations Maximum iterations per point
 * @param tolerance Convergence threshold for |∇ρ ⋅ tangent|
 * @return Number of successfully converged points
 */
// Removed: duplicate rtsafe_iteration wrapper - replaced with unified version above

// ============================================================================
//  BOUNDARY CONSTRUCTION (Stage D) - Connecting surface points into boundaries
// ============================================================================

namespace cpplib {

/**
 * Surface boundary for Bader zero-flux surfaces.
 *
 * Represents a closed polygon that forms part of the atom basin boundary.
 * Vertices should be ordered counter-clockwise when viewed from outside.
 */
struct SurfaceBoundary {
    uint32_t basin_a;       ///< Owner atom index
    uint32_t basin_b;       ///< Neighbor atom index
    std::vector<geometry::Point<double>> verts;   ///< Vertices in order (CCW from outside)

    /**
     * Create a new boundary with given basin pair.
     * @param ba Owner atom identifier
     * @param bb Neighbor atom identifier
     */
    SurfaceBoundary(uint32_t ba, uint32_t bb) : basin_a(ba), basin_b(bb) {}

    /**
     * Add a vertex to the boundary.
     */
    void add_vertex(const geometry::Point<double>& v) { verts.push_back(v); }

    /**
     * Compute area using cross-product method.
     * @return Surface area in Å²
     */
    double compute_area() const;

    /**
     * Check if boundary has valid polygon (at least 3 vertices).
     */
    bool is_valid() const { return verts.size() >= 3; }

    /**
     * Update area field after vertex list changes.
     */
    void update_area() { compute_area(); }
};

/**
 * Connect surface points into boundaries using basin pairing.
 *
 * Groups refined surface points by (basin_a, basin_b) pairs and connects
 * adjacent points along the zero-flux boundary between atoms.
 *
 * @param sp Refined surface points from rtsafe iteration (Stage C)
 * @param point_spacing Minimum expected spacing between points (Å)
 * @return Vector of SurfaceBoundary objects
 */
inline std::vector<SurfaceBoundary> build_boundaries_from_points(
    const std::vector<SurfacePoint>& sp,
    double point_spacing = 0.5
) {
    return cpplib::surface::build_boundaries_from_points(sp, point_spacing);
}

/**
 * Compute area of surface boundary using cross-product method.
 */
inline double SurfaceBoundary::compute_area() const {
    if (verts.size() < 3) return 0.0;

    const auto& v0 = verts[0];
    double sum_x(0.0), sum_y(0.0), sum_z(0.0);

    for (size_t i = 1, n = verts.size(); i < n; ++i) {
        const auto& vi = verts[i];
        const auto& vi1 = verts[(i + 1) % n];

        sum_x += v0[1] * vi[2] - v0[2] * vi[1];
        sum_y += v0[2] * vi1[0] - v0[0] * vi[2];
        sum_z += v0[0] * vi1[1] - v0[1] * vi[0];
    }

    return 0.5 * std::sqrt(sum_x * sum_x + sum_y * sum_y + sum_z * sum_z);
}

}  // namespace cpplib


// ============================================================================
//  EXAMPLE USAGE PATTERNS
// ============================================================================

/**
 * Typical workflow pattern for surface construction:
 *
 * ```cpp
 * // 1. Build Voronoi tessellation
 * voronoi::VoronoiFused vf = voronoi::VoronoiFused::build(cell_list, frac_to_cart);
 *
 * // 2. Initialize surface from Voronoi faces
 * std::vector<SurfacePoint> seeds = initialize_surface_from_voronoi(vf, psoa, ed);
 *
 * // 3. Perform rtsafe iteration to find zero-flux boundary (Stage C)
 * auto refined = cpplib::rtsafe_iteration(seeds, psoa, 50, 1e-8);
 *
 * // 4. Build boundaries by connecting points (Stage D)
 * std::vector<SurfaceBoundary> boundaries = build_boundaries_from_points(refined);
 * ```
 */

/**
 * Point density configuration example:
 *
 * ```cpp
 * constexpr double density_1 = cpplib::surface_constants::POINT_DENSITY_A2;  // 0.5 Å²/point
 * constexpr double density_2 = cpplib::surface_constants::POINT_DENSITY_A2 * 2.0;  // finer (0.25)
 * constexpr double density_3 = cpplib::surface_constants::POINT_DENSITY_A2 / 2.0;  // coarser (1.0)
 *
 * auto points = poisson_disk_sample(vertices, target_count,
 *                                   std::sqrt(surface_constants::POINT_DENSITY_A2));
 * ```
 */


// ============================================================================
//  PERFORMANCE NOTES
// ============================================================================

/**
 * @file Performance considerations for surface construction
 *
 * Memory:
 * - Each SurfacePoint requires: 6 × sizeof(PointType) + 2 × sizeof(uint32_t) ≈ 48 bytes (for double)
 * - For N shared faces, expect ~2N points with density=0.5 Å²
 * - Typical crystal (~100 atoms): 50-200 faces → 100-400 surface points
 * - SurfaceBoundary overhead: ~32 bytes + vertex storage per boundary
 *
 * Time complexity:
 * - Stage A (constants/structures): O(1)
 * - Stage B (Poisson sampling): O(n_points × max_iterations) where n_points = target per polygon
 * - Stage C (rtsafe iteration): O(points × grid_accesses × max_iters) per point
 *   - Grid access cost depends on SG spatial lookup implementation
 * - Stage D (boundary construction): O(points log points) due to map grouping
 *
 * Typical memory usage:
 * - 200 surface points with SurfacePoint struct: ~10 KB
 * - 50 boundaries with avg 10 vertices each: ~160 bytes overhead + vertex data
 *
 * Optimization opportunities:
 * - Spatial hashing for neighbor lookup in Poisson sampling
 * - SIMD vectorization of distance calculations in sampling
 * - Parallel processing of independent polygon groups (Stage D)
 * - Early termination in rtsafe when convergence achieved
 * - Cache-friendly SoA layout for electron density grid access
 *
 */
