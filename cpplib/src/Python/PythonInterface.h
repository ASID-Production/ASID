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

#define Py_LIMITED_API 0x030A0000
#include <Python.h>

#include <array>
#include <concepts>
#include <list>
#include <utility>
#include <vector>

#include "../Functions/Functions.h"
#include "../BaseHeaders/BaseTypes.h"
#include "../Classes/Bond.h"
#include "../Classes/Cluster.h"
#include "../Classes/Geometry.h"
#include "../Classes/Voronoi.h"

namespace py_util {

	template<std::integral I> inline PyObject* convert(I v);
	inline PyObject* convert(cpplib::basic_types::FloatingPointType v);
	inline PyObject* convert(const std::string& v);

	template<typename T> PyObject* convert(const std::vector<T>& vec);
	template<typename T, std::size_t N> inline PyObject* convert(const std::array<T, N>& arr);

	template <typename T> inline PyObject* convert(const cpplib::geometry::Point<T>& p);

	inline PyObject* convert(const cpplib::Bond& b);
	inline PyObject* convert(const cpplib::BondEx& b);

	inline PyObject* convert(const cpplib::cluster_detail::ClusterData& cd);

	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::EdgeIn& e);
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::PolygonIn& p);
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::Polyhedron& p);
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused& vf);

	inline bool add_to_dict(const char* key, PyObject* val, PyObject* dict) {
		if (!val) return false;
		if (PyDict_SetItemString(dict, key, val) != 0) [[unlikely]] {
			Py_DECREF(val);
			return false;
		}
		Py_DECREF(val);
		return true;
	}

	inline PyObject* convert(cpplib::basic_types::FloatingPointType v) {
		return PyFloat_FromDouble(static_cast<double>(v));
	}
	inline PyObject* convert(const std::string& v) {
		return PyUnicode_FromString(v.c_str());
	}
	template<std::integral I>
	inline PyObject* convert(I v) {
		if constexpr (std::is_same_v<I, bool>) {
			return PyBool_FromLong(static_cast<long>(v));
		} else if constexpr (std::is_unsigned_v<I>) {
			// (size_t, uint64_t)
			return PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(v));
		} else {
			// (int, int64_t)
			return PyLong_FromLongLong(static_cast<long long>(v));
		}
	}

	template <typename T>
	inline PyObject* convert(const cpplib::geometry::Point<T>& p) {
		PyObject* tuple = PyTuple_New(3);
		if (!tuple) return nullptr;

		for (Py_ssize_t i = 0; i < 3; i++) {
			PyObject* val = convert(p[i]);
			if (!val) [[unlikely]] {
				Py_DECREF(tuple);
				return nullptr;
			}

			if (PyTuple_SetItem(tuple, i, val) != 0) [[unlikely]] {
				Py_DECREF(val);
				Py_DECREF(tuple);
				return nullptr;
			}
		}
		return tuple;
	}

	// (std::vector -> List)
	template<typename T>
	PyObject* convert(const std::vector<T>& vec) {
		Py_ssize_t n = static_cast<Py_ssize_t>(vec.size());
		PyObject* list = PyList_New(n);
		if (!list) return nullptr;
		for (Py_ssize_t i = 0; i < n; ++i) {
			PyObject* item = convert(vec[i]);
			if (!item) [[unlikely]] {
				Py_DECREF(list); 
				return nullptr;
			}
			if (PyList_SetItem(list, i, item) != 0) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}
		}
		return list;
	}
	// (std::list -> List)
	template<typename T>
	PyObject* convert(const std::list<T>& l) {
		Py_ssize_t n = static_cast<Py_ssize_t>(l.size());
		PyObject* list = PyList_New(n);
		if (!list) return nullptr;
		int i = 0;
		for (const auto& e : l) {
			PyObject* item = convert(e);
			if (!item) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}
			if (PyList_SetItem(list, i, item) != 0) [[unlikely]] {

				Py_DECREF(list);
				return nullptr;
			}
			i++;
		}
		return list;
	}
	// (std::array -> List)
	template<typename T, std::size_t N>
	inline PyObject* convert(const std::array<T, N>& arr) {
		PyObject* list = PyList_New(static_cast<Py_ssize_t>(N));
		if (!list) return nullptr;

		for (std::size_t i = 0; i < N; ++i) {
			PyObject* item = convert(arr[i]);

			if (!item) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}

			if (PyList_SetItem(list, static_cast<Py_ssize_t>(i), item) != 0) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}
		}
		return list;
	}

	// (Bond -> Tuple(Int,Int))
	inline PyObject* convert(const cpplib::Bond& b) {
		return Py_BuildValue("(NN)", convert(b.first), convert(b.second));
	}
	// (BondEx -> Tuple(Int,Int))
	inline PyObject* convert(const cpplib::BondEx& b) {
		return Py_BuildValue("(NNN)", convert(b.first), convert(b.second), convert(b.length));
	}

	// (cluster_detail::ClusterData -> Dict)
	inline PyObject* convert(const cpplib::cluster_detail::ClusterData& cd) {
		auto N = static_cast<Py_ssize_t>(cd.indexes.size());
		PyObject* list = PyList_New(N);
		if (!list) return nullptr;

		for (Py_ssize_t i = 0; i < N; ++i) {
			PyObject* dict = PyDict_New();
			if (!dict) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}

			if (!add_to_dict("index", convert(cd.indexes[i]), dict) ||
				!add_to_dict("type", convert(cd.types[i]), dict) ||
				!add_to_dict("point_frac", convert(cd.points[i]), dict) ||
				!add_to_dict("symmref", convert(cd.symm_indexes[i]), dict) ||
				!add_to_dict("shift", convert(cd.shifts[i]), dict)) [[unlikely]] {
				Py_DECREF(list);
				Py_DECREF(dict);
				return nullptr;
			}

			if (PyList_SetItem(list, i, dict) != 0) [[unlikely]] {
				Py_DECREF(list);
				return nullptr;
			}
		}
		return list;
	}

	// (VoronoiFused::EdgeIn -> Dict)
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::EdgeIn& e) {
		return convert(e.vert_ids);
	}
	// (VoronoiFused::PolygonIn -> Dict)
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::PolygonIn& p) {
		PyObject* dict = PyDict_New();
		if (!dict) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		if (!add_to_dict("area",        convert(p.area), dict) ||
			!add_to_dict("vertices",    convert(p.vert_ids), dict) ||
			!add_to_dict("shift",       convert(p.second_shift), dict) ||
			!add_to_dict("edges",       convert(p.edge_ids), dict) ||
			!add_to_dict("solid_angle", convert(p.solidangle), dict) ||
			!add_to_dict("atoms",       convert(p.atom_ids), dict)) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		return dict;
	}
	// (VoronoiFused::Polyhedron -> Dict)
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused::Polyhedron& p) {
		PyObject* dict = PyDict_New();
		if (!dict) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		if (!add_to_dict("vertices", convert(p.vert_ids), dict) ||
			!add_to_dict("center",   convert(p.center), dict) ||
			!add_to_dict("edges",    convert(p.edge_ids), dict) ||
			!add_to_dict("volume",   convert(p.volume), dict) ||
			!add_to_dict("polygons", convert(p.poly_ids), dict)) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		return dict;
	}

	// (VoronoiFused)
	inline PyObject* convert(const cpplib::voronoi::VoronoiFused& vf) {
		PyObject* dict = PyDict_New();
		if (!dict) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		if (!add_to_dict("vertices",  convert(vf.vertices), dict) ||
			!add_to_dict("edges",     convert(vf.edges), dict) ||
			!add_to_dict("polygons",  convert(vf.polygons), dict) ||
			!add_to_dict("polyhedra", convert(vf.polyhedra), dict)) [[unlikely]] {
			Py_DECREF(dict);
			return nullptr;
		}

		return dict;
	}
}