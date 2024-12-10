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
#include <mutex>
#include <vector>
#include <bitset>
#include "../Classes/Geometry.h"
#include "../Classes/FindMolecules.h"
#include "../Functions/AllInOneAndCurrent.h"
namespace cpplib {
	class SearchDataInterface {
	public:
		using MoleculeIndex = currents::MoleculeIndex;
		using RawVector = std::vector<const char*>;
		using MultiflagType = std::bitset<mend_size>;
		using size_type = RawVector::size_type;
	private:
		size_type iterator_ = 0;
		const size_type size_;
		const RawVector rawdata_;
		const MultiflagType multiflag_;

		std::mutex mutexIN_;
		std::mutex mutexOUT_;

		std::vector<int> ret_;

	public:
		SearchDataInterface() = delete;
		explicit SearchDataInterface(std::vector<const char*>&& rawdata, MultiflagType&& multiflag) noexcept
			: rawdata_(rawdata), size_(rawdata.size()), multiflag_(multiflag) {
			if (size_ >= 1024)
				ret_.reserve(1024);
			else
				ret_.reserve(static_cast<size_t>(size_));
		};
		inline size_type size() const noexcept {
			return size_;
		}
		const char* getNext() {
			size_type iter;
			do {
				iter = getNextIterator();
				if (iter == size_type(-1))
					return nullptr;
			} while (rawdata_[iter][0] == '\0');

			return rawdata_[iter];
		}
		inline const MultiflagType& getMulty() const noexcept {
			return multiflag_;
		}
		void push_result(const MoleculeIndex molecularID) {
			std::lock_guard<std::mutex> lock(mutexOUT_);
			ret_.push_back(molecularID);
		}
		std::vector<int>&& getAllResults() noexcept {
			std::lock_guard<std::mutex> lock(mutexOUT_);
			return std::move(ret_);
		}
	private:
		size_type getNextIterator() {
			std::lock_guard<std::mutex> lock(mutexIN_);
			if (iterator_ == size_)
				return size_type(-1);
			size_type i = iterator_;
			iterator_++;
			return i;
		}
	};

	struct ParseData {
		using FAMStructType = FAM_Struct;
		using FAMCellType = FAM_Cell;

		using AtomType = FAMStructType::AtomType;
		using PointType = FAMStructType::PointType;
		using SymmType = FAMCellType::SymmType;

		ParseData(FAMStructType& fs, const std::vector<int>& types, const std::vector<float>& xyz) {
			const auto types_s = types.size();
			fs.points.reserve(types_s);
			fs.types.reserve(types_s);
			fs.sizePoints = types_s;
			fs.sizeUnique = types_s;
			for (int i = 0, i3 = 0; i < types_s; i++, i3 += 3) {
				fs.types.emplace_back(types[i]);
				fs.points.emplace_back(xyz[i3], xyz[i3 + 1], xyz[i3 + 2]);
			}
			fs.parseIndex.resize(types_s);
			std::iota(fs.parseIndex.begin(), fs.parseIndex.end(), 0); // Fill with 0, 1...
		}
		ParseData(FAMStructType& fs, FAMStructType::AtomContainerType && types, FAMStructType::PointConteinerType && points) {
			const auto types_s = types.size();
			fs.types = std::move(types);
			fs.points = std::move(points);
			fs.sizePoints = types_s;
			fs.sizeUnique = types_s;

			fs.parseIndex.resize(types_s);
			std::iota(fs.parseIndex.begin(), fs.parseIndex.end(), 0); // Fill with 0, 1...
		}

		ParseData(FAMStructType& fs, FAMCellType& fc, const std::vector<const char*>& symm, FAMStructType::AtomContainerType&& types, FAMStructType::PointConteinerType&& points)
			: ParseData(fs, std::move(types), std::move(points))
		{
			const int symm_s = symm.size();
			std::vector<SymmType> symmv;
			symmv.reserve(symm_s - 1);
			for (int i = 1; i < symm_s; i++) {
				symmv.emplace_back(symm[i]);
			}

			fc.GenerateSymm(fs, symmv, true);
			fs.sizePoints = fs.points.size();
			fs.types.reserve(fs.sizePoints);
			for (decltype(fs.sizeUnique) i = fs.sizeUnique; i < fs.sizePoints; i++) {
				fs.types.emplace_back(fs.types[fs.parseIndex[i]]);
			}
		}
	};
}