// Surface Construction - QTAIM zero-flux surface generation using Voronoi polyhedra approach
// Implements:
//   Stage A: Constants and data structures
//   Stage B: Initialize from Voronoi faces (Poisson disk sampling)
//   Stage C: rtsafe iteration with lambda-criteria for zero-flux boundary
//   Stage D: Build boundaries (placeholder, using surface:: module)

#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <array>
#include <iomanip>

#include "../Classes/Geometry.h"
#include "../Classes/Splines.h"
#include "../Classes/BaderOperator.h"
#include "../Classes/SurfaceConstruction_types.h"

namespace cpplib {

// ========================================================================
// Hessian Eigenvalue Analysis Structure (Stage 2)
// ========================================================================

struct SurfacePointEigenAnalysis {
	TripleDouble lambda_max;    ///< Largest eigenvalue λ₃
	TripleDouble lambda_mid;    ///< Middle eigenvalue λ₂
	TripleDouble lambda_min;    ///< Smallest eigenvalue λ₁
	PointType eigenvector_max;  ///< Dominant eigenvector e₃ (normal candidate)

	bool is_reliable = false;
	double lambda_ratio = 0.0;
	std::string flag_type;

	void analyze(const TripleDouble::value_type& ev0, const TripleDouble::value_type& ev1,
				            const TripleDouble::value_type& ev2, const PointType& eigenvector_max_val,
				            double epsilon = 1e-30) {

		std::array<TripleDouble, 3> values = {TripleDouble(ev0), TripleDouble(ev1), TripleDouble(ev2)};
		std::sort(values.begin(), values.end(), [](const TripleDouble& a, const TripleDouble& b) {
			return a.value_ > b.value_;
		});

		lambda_max = values[0];
		lambda_mid = values[1];
		lambda_min = values[2];
		eigenvector_max = eigenvector_max_val;

		const auto abs_max = std::max(std::abs(lambda_max),
		                              std::max(std::abs(lambda_mid), std::abs(lambda_min)));

		if (lambda_max > 0.0 && lambda_max >= abs_max * 0.99) {
			lambda_ratio = lambda_max / (std::max(std::abs(lambda_mid), epsilon));

			if (lambda_ratio > 1.5 && (lambda_mid < -epsilon || std::abs(lambda_mid) < abs_max * 0.1)) {
				is_reliable = true;
				flag_type = "ok";
			} else if (lambda_ratio < 1.3) {
				is_reliable = false;
				flag_type = "degenerate";
			} else {
				if ((lambda_max > epsilon && lambda_mid > -epsilon) ||
				    (lambda_max > epsilon && lambda_min > -epsilon)) {
					is_reliable = false;
					flag_type = "two_positive";
				} else if (lambda_min > -epsilon) {
					is_reliable = false;
					flag_type = "all_positive";
				} else {
					is_reliable = true;
					flag_type = "ok";
				}
			}
		} else if (lambda_max > epsilon && lambda_mid > -epsilon) {
			is_reliable = false;
			flag_type = "two_positive";
		} else if (lambda_max <= -epsilon) {
			is_reliable = false;
			flag_type = "all_negative";
		} else {
			is_reliable = false;
			flag_type = "ambiguous";
		}
	};

};

// ========================================================================
// Helper Structure for Iteration Data (Stage 3)
// ========================================================================

struct RtsafeIterationData {
	PointType pos;
	PointType normal;
	BaderOperator::CriticalPoint crit;
	PointType pos_prev;
	auto hess_eig = std::make_tuple(
		std::optional<TripleDouble>(), std::optional<TripleDouble>(),
		std::optional<TripleDouble>(), std::optional<PointType>()
	);

	uint32_t iteration_count = 0;
	bool converged = false;
	SurfacePointEigenAnalysis analysis;

private:
	template <typename Tuple>
	void set_hess_eig(const Tuple& t) {
		std::get<0>(hess_eig) = std::make_optional<TripleDouble>(std::get<0>(t));
		std::get<1>(hess_eig) = std::make_optional<TripleDouble>(std::get<1>(t));
		std::get<2>(hess_eig) = std::make_optional<TripleDouble>(std::get<2>(t));
		std::get<3>(hess_eig) = std::make_optional<PointType>(std::get<3>(t));
	}

public:
	void update_hess_eig(const TripleDouble& v0, const TripleDouble& v1,
	                     const TripleDouble& v2, const PointType& vec2) {
		std::get<0>(hess_eig) = v0;
		std::get<1>(hess_eig) = v1;
		std::get<2>(hess_eig) = v2;
		std::get<3>(hess_eig) = vec2;
	}

	TripleDouble get_ev(size_t idx) const { return *std::get<idx>(hess_eig); }
	PointType get_vec2() const { return std::get<3>(hess_eig); }
};

static void prepare_iteration_initialization(
	std::vector<SurfacePoint>& sp,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	uint32_t max_iterations,
	double tolerance,
	RtsafeIterationData& iter_data,
	std::unordered_map<std::string, uint32_t>* flag_counts = nullptr
) {
	if (sp.empty()) return;

	iter_data.pos = sp[0].pos;
	iter_data.normal = sp[0].normal;
	iter_data.crit = BaderOperator::CriticalPoint(iter_data.pos, psoa);
	iter_data.crit.CalculateEDinPoint();

	flag_counts["total"] = 0;
	flag_counts["ok"] = 0;
	flag_counts["degenerate"] = 0;
	flag_counts["two_positive"] = 0;
	flag_counts["all_negative"] = 0;
	flag_counts["ambiguous"] = 0;
	flag_counts["critical_point"] = 0;
}

static SurfacePointEigenAnalysis evaluate_hessian_eigenvalues(
	const RtsafeIterationData& data,
	std::vector<SurfacePointEigenAnalysis>* return_analysis = nullptr
) {
	auto hess_eig = data.crit.value_.hess.EigenvaluesAndVectors();

	SurfacePointEigenAnalysis analysis;
	analysis.analyze(hess_eig.values[0], hess_eig.values[1], hess_eig.values[2], hess_eig.vectors[2]);

	if (return_analysis) {
		return_analysis->clear();
		return_analysis->push_back(analysis);
	} else {
		return analysis;
	}
}

static bool build_initial_bracket(
	RtsafeIterationData& data,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	TripleDouble& t_low, TripleDouble& t_high,
	std::vector<std::pair<TripleDouble, TripleDouble>>& f_low, std::vector<std::pair<TripleDouble, TripleDouble>>& f_high
) {
	const auto& grad = data.crit.value_.grad;
	TripleDouble f0(0.0);
	f0[0] = grad[0] * data.normal[0];
	f0[1] = grad[1] * data.normal[1];
	f0[2] = grad[2] * data.normal[2];

	constexpr double BRACKET_TOLERANCE = 1e-6;
	if (std::abs(f0.value_) < BRACKET_TOLERANCE * 10) {
		t_low = -BRACKET_TOLERANCE;
		t_high = BRACKET_TOLERANCE;
		return true;
	}

	TripleDouble hess_normal(0.0);
	hess_normal[0] = data.crit.value_.hess[0] * data.normal[0] * data.normal[0] +
		               data.crit.value_.hess[1] * data.normal[1] * data.normal[1] +
		               data.crit.value_.hess[2] * data.normal[2] * data.normal[2] +
		               2.0 * (data.crit.value_.hess[0]*data.normal[1]*data.normal[0] +
	                          data.crit.value_.hess[0]*data.normal[2]*data.normal[0] +
	                          data.crit.value_.hess[1]*data.normal[2]*data.normal[1]);

	const bool increasing = hess_normal.value_ > 0;

	if (f0.value_ > 0.0) {
		t_low = (increasing) ? -1.0 : 0.0;
		t_high = (increasing) ? 0.0 : 1.0;
	} else {
		t_low = (increasing) ? 0.0 : -1.0;
		t_high = (increasing) ? 1.0 : 0.0;
	}

	const PointType r_low_est = data.pos + data.normal * t_low.value_;
	const PointType r_high_est = data.pos + data.normal * t_high.value_;

	BaderOperator::CriticalPoint crit_low(r_low_est, psoa);
	crit_low.CalculateEDinPoint();
	f_low.emplace_back(TripleDouble(crit_low.value_.grad[0] * data.normal[0]),
	                   TripleDouble(crit_low.value_.grad[1] * data.normal[1]),
	                   TripleDouble(crit_low.value_.grad[2] * data.normal[2]));

	BaderOperator::CriticalPoint crit_high(r_high_est, psoa);
	crit_high.CalculateEDinPoint();
	f_high.emplace_back(TripleDouble(crit_high.value_.grad[0] * data.normal[0]),
	                    TripleDouble(crit_high.value_.grad[1] * data.normal[1]),
	                    TripleDouble(crit_high.value_.grad[2] * data.normal[2]));

	return false;
}

static bool perform_newton_search(
	RtsafeIterationData& data,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	TripleDouble t_low, TripleDouble t_high,
	std::pair<TripleDouble, TripleDouble> f_low, std::pair<TripleDouble, TripleDouble> f_high
) {
	constexpr double MAX_NEWTON_STEP = 3.0;
	for (uint32_t niter = 0; niter < 10; ++niter) {
		if (std::abs(t_high - t_low) < 1e-12) break;

		const TripleDouble t_mid = 0.5 * (t_low + t_high);
		const PointType r_trial = data.pos + data.normal * t_mid.value_;

		BaderOperator::CriticalPoint crit_trial(r_trial, psoa);
		crit_trial.CalculateEDinPoint();

		TripleDouble f_val(0.0);
		f_val[0] = crit_trial.value_.grad[0] * data.normal[0];
		f_val[1] = crit_trial.value_.grad[1] * data.normal[1];
		f_val[2] = crit_trial.value_.grad[2] * data.normal[2];

		const double f_norm = std::sqrt(f_val[0]*f_val[0] + f_val[1]*f_val[1] + f_val[2]*f_val[2]);
		if (f_norm < tolerance) {
			data.pos = r_trial;
			data.crit = crit_trial;
			data.normal = hess_eig.vectors[2];

			SurfacePointEigenAnalysis analysis_final;
			const auto& hess_eig_final = data.crit.value_.hess.EigenvaluesAndVectors();
			analysis_final.analyze(hess_eig_final.values[0], hess_eig_final.values[1],
			                      hess_eig_final.values[2], hess_eig_final.vectors[2]);

			if (return_analysis) {
				return_analysis->clear();
				return_analysis->push_back(analysis_final);
			} else {
				data.analysis = analysis_final;
			}

			break;
		}

		TripleDouble df(0.0);
		df[0] = data.crit.value_.hess[0] * data.normal[0] * data.normal[0] +
		       data.crit.value_.hess[1] * data.normal[1] * data.normal[1] +
		       data.crit.value_.hess[2] * data.normal[2] * data.normal[2] +
		       2.0 * (data.crit.value_.hess[0]*data.normal[1]*data.normal[0] +
	                  data.crit.value_.hess[0]*data.normal[2]*data.normal[0] +
	                  data.crit.value_.hess[1]*data.normal[2]*data.normal[1]);

		if (std::abs(df.value_) < 1e-30) {
			df = TripleDouble(0.1);
		}

		TripleDouble dt(-f_val / df);
		TripleDouble t_new = t_mid + dt;

		if ((t_new.value_ > t_low && t_new.value_ < t_high) ||
		    (std::abs(dt.value_) < 0.5)) {
			t_low = std::min(t_low, t_new);
			t_high = std::max(t_high, t_new);
		} else {
			if (f_val.value_ * f_low.first.value_ > 0) {
				t_low = t_mid;
				f_low = std::make_pair(f_val, f_low.second);
			} else {
				t_high = t_mid;
			}
		}

		if (std::abs(t_mid.value_) > MAX_NEWTON_STEP) {
			break;
		}
	}

	return true;
}

static RtsafeIterationData finalize_iteration(
	RtsafeIterationData& data,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	bool converged,
	TripleDouble t_low, TripleDouble t_high,
	double tolerance
) {
	constexpr double MAX_NEWTON_STEP = 3.0;

	const auto& grad = data.crit.value_.grad;

	data.pos = data.pos + data.normal * 0.5 * (t_low + t_high);

	double disp_sq(0);
	for(int d=0; d<3; ++d) {
		disp_sq += (data.pos[d] - data.pos_prev[d]) * (data.pos[d] - data.pos_prev[d]);
	}

	if (std::sqrt(disp_sq) > MAX_NEWTON_STEP) {
		TripleDouble f_current(0.0);
		f_current[0] = grad[0] * data.normal[0];
		f_current[1] = grad[1] * data.normal[1];
		f_current[2] = grad[2] * data.normal[2];

		PointType tangent_grad = grad;
		tangent_grad[0] -= f_current.value_ * data.normal[0];
		tangent_grad[1] -= f_current.value_ * data.normal[1];
		tangent_grad[2] -= f_current.value_ * data.normal[2];

		const double tg_norm_sq = tangent_grad[0]*tangent_grad[0] +
		                           tangent_grad[1]*tangent_grad[1] +
		                           tangent_grad[2]*tangent_grad[2];

		if (tg_norm_sq > 1e-30) {
			tangent_grad[0] /= std::sqrt(tg_norm_sq);
			tangent_grad[1] /= std::sqrt(tg_norm_sq);
			tangent_grad[2] /= std::sqrt(tg_norm_sq);
			data.pos[0] += tangent_grad[0] * tolerance;
			data.pos[1] += tangent_grad[1] * tolerance;
			data.pos[2] += tangent_grad[2] * tolerance;
		}
	}

	data.crit = BaderOperator::CriticalPoint(data.pos, psoa);
	data.crit.CalculateEDinPoint();
	auto hess_eig_new = data.crit.value_.hess.EigenvaluesAndVectors();

	data.update_hess_eig(hess_eig_new.values[0], hess_eig_new.values[1],
		                 hess_eig_new.values[2], hess_eig_new.vectors[2]);

	SurfacePointEigenAnalysis analysis_new;
	const auto& ev0_new = hess_eig_new.values[0];
	const auto& ev1_new = hess_eig_new.values[1];
	const auto& ev2_new = hess_eig_new.values[2];
	const auto& ve2_new = hess_eig_new.vectors[2];

	analysis_new.analyze(ev0_new, ev1_new, ev2_new, ve2_new);

	data.analysis = analysis_new;

	PointType new_normal = data.normal;
	if (analysis_new.is_reliable) {
		new_normal = ve2_new;
	} else {
		const auto& grad_new = data.crit.value_.grad;
		const double gnorm_sq = grad_new[0]*grad_new[0] +
		                       grad_new[1]*grad_new[1] +
		                       grad_new[2]*grad_new[2];
		if (gnorm_sq > 1e-24) {
			new_normal[0] = grad_new[0] / std::sqrt(gnorm_sq);
			new_normal[1] = grad_new[1] / std::sqrt(gnorm_sq);
			new_normal[2] = grad_new[2] / std::sqrt(gnorm_sq);
		} else {
			new_normal = PointType(0.0, 0.0, 1.0);
		}
	}

	data.normal[0] = new_normal[0];
	data.normal[1] = new_normal[1];
	data.normal[2] = new_normal[2];

	const auto& grad_new = data.crit.value_.grad;
	TripleDouble f_tangent(0.0);
	f_tangent[0] = grad_new[0] * data.normal[0];
	f_tangent[1] = grad_new[1] * data.normal[1];
	f_tangent[2] = grad_new[2] * data.normal[2];

	if (std::sqrt(f_tangent[0]*f_tangent[0] + f_tangent[1]*f_tangent[1] + f_tangent[2]*f_tangent[2]) < tolerance) {
		data.converged = true;
	}

	return data;
}

static RtsafeIterationData perform_newton_iteration(
	RtsafeIterationData& data,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	uint32_t max_iterations,
	double tolerance,
	std::vector<SurfacePointEigenAnalysis>* return_analysis = nullptr
) {
	for (uint32_t iter = 0; iter < max_iterations; ++iter) {
		data.pos_prev[0] = data.pos[0];
		data.pos_prev[1] = data.pos[1];
		data.pos_prev[2] = data.pos[2];

		TripleDouble t_low(-0.5), t_high(0.5);
		std::vector<std::pair<TripleDouble, TripleDouble>> f_low_storage;
		std::vector<std::pair<TripleDouble, TripleDouble>> f_high_storage;

		if (!build_initial_bracket(data, psoa, t_low, t_high, f_low_storage, f_high_storage)) {
			continue;
		}

		perform_newton_search(data, psoa, t_low, t_high,
		                    f_low_storage.front(), f_high_storage.front());

		TripleDouble t_mid = 0.5 * (t_low + t_high);
		f_low = std::make_pair(TripleDouble(f_low_storage.front().first.value_),
		                       TripleDouble(f_low_storage.front().second.value_));
		f_high = std::make_pair(TripleDouble(f_high_storage.front().first.value_),
		                        TripleDouble(f_high_storage.front().second.value_));

		perform_newton_search(data, psoa, t_low, t_high, f_low, f_high);

		data.pos = data.pos + data.normal * 0.5 * (t_low + t_high);

		// Check convergence after each iteration
		const auto& grad = data.crit.value_.grad;
		TripleDouble f_tangent(0.0);
		f_tangent[0] = grad[0] * data.normal[0];
		f_tangent[1] = grad[1] * data.normal[1];
		f_tangent[2] = grad[2] * data.normal[2];

		if (std::sqrt(f_tangent[0]*f_tangent[0] + f_tangent[1]*f_tangent[1] + f_tangent[2]*f_tangent[2]) < tolerance) {
			return finalize_iteration(data, psoa, true, t_low, t_high, tolerance);
		}
	}

	return finalize_iteration(data, psoa, false, TripleDouble(0.0), TripleDouble(0.0), tolerance);

}

</function>
</tool_call>
<parameter=old_string>

static std::vector<SurfacePoint> rtsafe_iteration_with_lambda_criteria(
	std::vector<SurfacePoint>& sp,
	const BaderOperator::PointsSoA<TripleDouble::value_type>& psoa,
	uint32_t max_iterations = 50,
	double tolerance = 1e-8,
	std::vector<SurfacePointEigenAnalysis>* return_analysis = nullptr
) {
	if (return_analysis && !return_analysis->empty()) {
		return_analysis->clear();
		return_analysis->resize(sp.size());
	}

	for (size_t i = 0; i < sp.size(); ++i) {
		RtsafeIterationData iter_data;

		iter_data.pos = sp[i].pos;
		iter_data.normal = sp[i].normal;

		BaderOperator::CriticalPoint crit(iter_data.pos, psoa);
		crit.CalculateEDinPoint();

		iter_data.crit = crit;

		perform_newton_iteration(iter_data, psoa, max_iterations, tolerance, return_analysis ? &(*return_analysis)[i] : nullptr);

		const auto& grad_final = iter_data.crit.value_.grad;
		TripleDouble f_tangent(0.0);
		f_tangent[0] = grad_final[0] * iter_data.normal[0];
		f_tangent[1] = grad_final[1] * iter_data.normal[1];
		f_tangent[2] = grad_final[2] * iter_data.normal[2];

		if (std::sqrt(f_tangent[0]*f_tangent[0] + f_tangent[1]*f_tangent[1] + f_tangent[2]*f_tangent[2]) < tolerance) {
			break;
		}
	}

	std::cout << "\n=== rtsafe Iteration Complete ===" << std::endl;

	return sp;
}

} // namespace cpplib
