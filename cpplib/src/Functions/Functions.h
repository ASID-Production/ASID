#pragma once
#include <vector>
#include <string>

bool CompareGraph(const char* search1, const char* search2, const bool exact);
std::vector<int> SearchMain(const char* search, std::vector <const char*>&& data, const int np, const bool exact);
std::string FindMoleculesInCell(const std::array<float, 6>& unit_cell, std::vector<const char*>& symm, std::vector<int>& types, std::vector<float>& xyz);
std::string FindMoleculesWithoutCell(const std::vector<int>& types, std::vector<float>& xyz);

const char* FindDistanceWC(const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value);
const char* FindDistanceIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value);
const char* FindAngleWC(const int* types, const float* xyz, const int types_s,
						const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123);
const char* FindAngleIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
						const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123);
const char* FindTorsionWC(const int* types, const float* xyz, const int types_s,
						  const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
						  const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234);
const char* FindTorsionIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
						  const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
						  const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234);
