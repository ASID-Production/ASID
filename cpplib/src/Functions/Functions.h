#pragma once
#include <vector>
#include <string>
#include "AllInOneAndCurrent.h"
bool CompareGraph(const char* search1,
					const char* search2,
					const bool exact);

std::vector<int> SearchMain(const char* search,
							std::vector <const char*>&& data,
							const int np,
							const bool exact);

std::pair<std::string, std::vector<std::pair<std::vector<std::tuple<cpplib::currents::PointType, cpplib::currents::AtomIndex>>, int>>>
	FindMoleculesInCell(const std::array<float, 6>& unit_cell,
						std::vector<const char*>& symm,
						std::vector<int>& types,
						std::vector<float>& xyz);

std::string FindMoleculesWithoutCell(const std::vector<int>& types,
										std::vector<float>& xyz);


std::string FindDistanceWC(std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 2>& type,
							const std::pair<float, float>& value);

std::string FindDistanceIC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 2>& type,
							const std::pair<float, float>& value);

std::string FindAngleWC(std::vector<int>& types,
						std::vector<float>& xyz,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a);

std::string FindAngleIC(const std::array<float, 6>& unit_cell,
						std::vector<const char*>& symm,
						std::vector<int>& types,
						std::vector<float>& xyz,
						const std::array<int, 3>& type,
						const std::array<std::pair<float, float>, 2>& value_d,
						const std::pair<float, float>& value_a);

std::string FindTorsionWC(std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 4>& type,
							const std::array<std::pair<float, float>, 3>& value_d,
							const std::array<std::pair<float, float>, 2>& value_a,
							const std::pair<float, float>& value_t);

std::string FindTorsionIC(const std::array<float, 6>& unit_cell,
							std::vector<const char*>& symm,
							std::vector<int>& types,
							std::vector<float>& xyz,
							const std::array<int, 4>& type,
							const std::array<std::pair<float, float>, 3>& value_d,
							const std::array<std::pair<float, float>, 2>& value_a,
							const std::pair<float, float>& value_t);
