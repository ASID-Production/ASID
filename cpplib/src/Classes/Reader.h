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

#include <iostream>
#include <vector>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <system_error>

#include "../Classes/Splines.h"

namespace cpplib {

    class DensityParser {
    public:
        struct AtomData {
            std::vector<double> values;
            double r_min;
            double r_max;
            int element;
            uint32_t n_knots;
        };
    private:
        [[nodiscard]] static constexpr const char* skip_spaces(const char* ptr, const char* end) noexcept {
            while (ptr < end && (*ptr == ' ' || *ptr == '\t' || *ptr == '\n' || *ptr == '\r')) {
                ++ptr;
            }
            return ptr;
        }

        // Universal numeric parsing template. Optimized perfectly by the compiler.
        template <typename T>
        [[nodiscard]] static const char* parse_value(const char* ptr, const char* end, T& value) {
            ptr = skip_spaces(ptr, end);
            if (ptr == end) {
                throw std::runtime_error("Parser error: Unexpected End of File (EOF).");
            }

            // std::from_chars is locale-independent, working directly with ASCII characters in memory
            auto [next_ptr, ec] = std::from_chars(ptr, end, value);
            if (ec != std::errc{}) {
                throw std::runtime_error("Parser error: Failed to convert numerical value.");
            }
            return next_ptr;
        }

    public:
        [[nodiscard]] static std::vector<RadialSpline> parse_file(const std::filesystem::path& filepath) {
            if (!std::filesystem::is_regular_file(filepath)) {
                throw std::runtime_error("File not found: " + filepath.string());
            }

            // 1. Block Binary IO: Determine file size
            const std::uintmax_t file_size = std::filesystem::file_size(filepath);

            // Allocate a contiguous buffer (2-3 MB size fits perfectly into CPU L3 cache)
            std::vector<char> buffer(file_size);

            // Open strictly in binary mode to prevent cross-platform issues with \r\n on Windows
            std::ifstream file(filepath, std::ios::binary);
            if (!file.is_open()) {
                throw std::runtime_error("Failed to open file: " + filepath.string());
            }

            // Single system call to read the entire file into memory
            file.read(buffer.data(), file_size);

            // Setup Zero-Copy Sliding Pointers
            const char* ptr = buffer.data();
            const char* end = ptr + buffer.size();

            // 2. Parse global header
            int max_element = 0;
            ptr = parse_value(ptr, end, max_element);

            // Pre-allocate the top-level vector (eliminates redundant heap reallocations)
            std::vector<RadialSpline> rs_ret;
            rs_ret.resize(max_element + 1);

            AtomData atom;
            atom.values.reserve(RadialSpline::SIZE);

            // 3. Atom block parsing pipeline
            while (true) {
                ptr = skip_spaces(ptr, end);
                if (ptr == end) break; // Parsing successfully completed

                // Construct the object directly within the vector's allocated memory (In-place)

                // Read block metadata
                ptr = parse_value(ptr, end, atom.element);
                ptr = parse_value(ptr, end, atom.r_min);
                ptr = parse_value(ptr, end, atom.r_max);
                ptr = parse_value(ptr, end, atom.n_knots);

                // Validation against specification boundaries
                if (atom.n_knots > RadialSpline::SIZE) {
                    throw std::runtime_error("Specification violated: N_KNOTS > 576 for atom ID " + std::to_string(atom.element));
                }
                if (atom.element < 1 || atom.element > max_element) {
                    throw std::runtime_error("Invalid element_no: " + std::to_string(atom.element) + " (expected 1-" + std::to_string(max_element) + ")");
                }

                // Read the values array (N_KNOTS elements)
                const size_t total_values = atom.n_knots;
                atom.values.resize(atom.n_knots + 1);
                for (size_t i = 0; i < total_values; ++i) {
                    ptr = parse_value(ptr, end, atom.values[i]);
                }
                atom.values.back() = 0;
                rs_ret[atom.element] = RadialSpline(atom.r_min, atom.r_max, atom.n_knots, atom.values);
            }

            return rs_ret;
        }
    };
}
