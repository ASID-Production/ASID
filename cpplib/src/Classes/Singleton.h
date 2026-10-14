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
//
// ******************************************************************************************
#pragma once

#include <atomic>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <mutex>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace cpplib {

	/// <summary>
	/// Program-wide result/state code of cpplib operations. Shared by ALL classes
	/// (not only singletons): singleton initialization, data-file loading and any
	/// other API that must report a status without throwing exceptions.
	/// Zero value (OK) means success; every non-zero value describes what exactly
	/// went wrong during the operation.
	/// </summary>
	enum class STATE : unsigned char {
		OK = 0,             ///< Operation completed successfully / object fully loaded and valid.
		NOT_INITIALIZED,    ///< Object was not created yet ("несобран"): lazy singleton has no instance.
		FILE_NOT_FOUND,     ///< The data file does not exist or cannot be opened.
		BAD_SIZE,           ///< File content is malformed: wrong header count, truncated data, size mismatch.
		BAD_CONTENT,        ///< File exists with plausible size but contains invalid values / parse errors.
		IO_ERROR,           ///< Underlying stream reported a failure that is not covered above.
		INVALID_ARGUMENT,   ///< Caller passed an unusable argument (e.g. empty path).
		ERROR,              ///< Any other unspecified failure.
	};

	// ---- Aliases for other program-wide error codes: existing consumers may
	// refer to these names; they denote the same STATE values. ----
	inline constexpr auto WRONG_INPUT_DATA   = STATE::BAD_CONTENT;
	inline constexpr auto PATH_DOESNT_EXISTS = STATE::FILE_NOT_FOUND;
	inline constexpr auto NO_SUCH_FILE       = STATE::FILE_NOT_FOUND;

	[[nodiscard]] inline constexpr const char* to_string(STATE s) noexcept {
		switch (s) {
			case STATE::OK:            return "OK";
			case STATE::NOT_INITIALIZED: return "NOT_INITIALIZED";
			case STATE::FILE_NOT_FOUND:  return "FILE_NOT_FOUND";
			case STATE::BAD_SIZE:        return "BAD_SIZE";
			case STATE::BAD_CONTENT:     return "BAD_CONTENT";
			case STATE::IO_ERROR:        return "IO_ERROR";
			case STATE::INVALID_ARGUMENT: return "INVALID_ARGUMENT";
			case STATE::ERROR:           return "ERROR";
		}
		return "UNKNOWN";
	}

	inline constexpr bool succeeded(STATE s) noexcept { return s == STATE::OK; }

	/// <summary>
	/// Translates a standard exception thrown while reading/parsing a data file
	/// into a program-wide STATE code (used instead of propagating exceptions out
	/// of the singleton interface).
	/// </summary>
	[[nodiscard]] inline STATE state_from_exception(const std::exception& e) noexcept {
		const std::string_view what = e.what();
		// DensityParser / ifstream messages produced by cpplib readers:
		if (what.find("File not found") != std::string_view::npos ||
		    what.find("Failed to open file") != std::string_view::npos) {
			return STATE::FILE_NOT_FOUND;
		}
		if (what.find("N_KNOTS") != std::string_view::npos ||
		    what.find("Unexpected End of File") != std::string_view::npos) {
			return STATE::BAD_SIZE;
		}
		if (what.find("Invalid element_no") != std::string_view::npos ||
		    what.find("Failed to convert numerical value") != std::string_view::npos ||
		    what.find("incorrect data") != std::string_view::npos) {
			return STATE::BAD_CONTENT;
		}
		return STATE::ERROR;
	}

	/// <summary>
	/// catch(...) variant: re-dispatches the active exception through
	/// std::rethrow_exception and maps it to a STATE code.
	/// </summary>
	[[nodiscard]] inline STATE state_from_current_exception() noexcept {
		try {
			throw;
		} catch (const std::exception& e) {
			return state_from_exception(e);
		} catch (...) {
			return STATE::ERROR;
		}
	}

	/// <summary>
	/// Absolute directory of the loaded executable or shared library hosting
	/// cpplib. Uses the native loader at runtime, never a source path or a CWD
	/// fallback. Failure to locate the module is reported as an exception.
	/// </summary>
	inline const std::filesystem::path& get_module_directory() {
		static const std::filesystem::path dir = [] {
			static const int anchor = 0;
#ifdef _WIN32
			HMODULE module = nullptr;
			if (!GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
				GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
				reinterpret_cast<LPCWSTR>(&anchor), &module)) {
				throw std::runtime_error("Could not locate cpplib module");
			}
			std::vector<wchar_t> buffer(256);
			for (;;) {
				const auto size = GetModuleFileNameW(module, buffer.data(), static_cast<DWORD>(buffer.size()));
				if (size == 0)
					throw std::runtime_error("Could not locate cpplib module");
				if (size < buffer.size())
					return std::filesystem::path(std::wstring(buffer.data(), size)).parent_path();
				buffer.resize(buffer.size() * 2);
			}
#else
			Dl_info info{};
			if (!dladdr(&anchor, &info) || !info.dli_fname || !*info.dli_fname)
				throw std::runtime_error("Could not locate cpplib module");
			return std::filesystem::canonical(info.dli_fname).parent_path();
#endif
		}();
		return dir;
	}

	/// <summary>
	/// Resolve data next to the loaded module. The canonical distribution paths
	/// module/BondLength.ini and module/ED.txt name files in that directory (the
	/// extension itself is installed in module/). Absolute paths are used as-is.
	/// </summary>
	inline std::filesystem::path get_module_file(const std::filesystem::path& relpath) {
		if (relpath.empty())
			return {};
		if (relpath.is_absolute())
			return std::filesystem::weakly_canonical(relpath);
		const auto local = *relpath.begin() == "module" ? relpath.lexically_relative("module") : relpath;
		return std::filesystem::weakly_canonical(get_module_directory() / local);
	}

	/// <summary>
	/// Per-singleton default relative data path. Each concrete singleton provides
	/// its own explicit specialization (see Distances / ElectronDensity), e.g.:
	///   template <> inline constexpr std::string_view default_relpath<Distances>
	///       = "module/BondLength.ini";
	/// The fallback keeps the standardized interface usable even without one.
	/// </summary>
	template <typename Derived>
	inline constexpr std::string_view default_relpath = "BondLength.ini";

	/// <summary>
	/// Standardized interface (CRTP) for every lazy static singleton in cpplib
	/// that owns read-only data loaded from a file (Distances, ElectronDensity, ...).
	///
	/// Contract provided to all derived singletons:
	///   - Derived::try_call(relpath = default_relpath) -> std::pair<const Derived&, STATE>
	///       THE primary entry point: exception-free. Creates the instance lazily
	///       on the FIRST call only (never eagerly at program startup); the
	///       RELATIVE path is resolved against the starting module directory
	///       (get_module_file()) exactly like Distances did before. Every
	///       subsequent call must resolve to the same file, otherwise it returns
	///       INVALID_ARGUMENT. If loading fails, the
	///       state is FILE_NOT_FOUND / BAD_SIZE / BAD_CONTENT / IO_ERROR instead
	///       of an exception being thrown. Thread-safe (mutex-guarded).
	///   - Derived::instance(relpath = ...) -> const Derived&
	///       Convenience accessor: same semantics, returns the reference only;
	///       throws if loading fails or a different path is requested.
	///       Use try_call() to receive loading errors as STATE values.
	///   - Derived::state()
	///       Current STATE of the singleton: OK after a successful load,
	///       NOT_INITIALIZED before the very first call ("несобран"), or the
	///       concrete failure code recorded during initialization.
	///   - Derived::exists()
	///       true if the instance has been successfully loaded; never triggers
	///       initialization.
	///   - Derived::reload(new_relpath) -> STATE
	///       Test/diagnostic hook: forces the singleton to re-read the data from
	///       a NEW location (resolved the same way as the initial relpath) and
	///       updates the stored STATE. Returns the resulting state. On failure
	///       the previously loaded data stay intact.
	///   - filename()
	///       absolute path the data were actually loaded from.
	///
	/// Objects are IMMUTABLE after loading: the instance is stored as const and
	/// all accessors hand out const references only — Distances and ED data are
	/// read-only by design.
	///
	/// Uniqueness guarantees: default construction, copy and move are deleted in
	/// this base, so the ONLY way to construct the object is via this interface
	/// (the protected constructor is reachable from Derived alone).
	///
	/// How to add a new data singleton:
	///   class Foo : public cpplib::SingletonBase<Foo> {
	///       explicit Foo(std::filesystem::path absfile, STATE& st); // loads data,
	///           // sets st = OK / FILE_NOT_FOUND / BAD_SIZE / ... (must not throw)
	///       ...data members and CONST accessors...
	///   };
	///   template <> inline constexpr std::string_view cpplib::default_relpath<Foo>
	///       = "module/foo.txt";
	/// </summary>
	template <typename Derived>
	class SingletonBase {
	public:
		SingletonBase() = delete;
		SingletonBase(const SingletonBase&) = delete;
		SingletonBase(SingletonBase&&) = delete;
		SingletonBase& operator=(const SingletonBase&) = delete;
		SingletonBase& operator=(SingletonBase&&) = delete;

		/// <summary>
		/// Exception-free standardized access to the unique static instance.
		/// See the class documentation above for the full contract.
		/// </summary>
		/// <param name="relpath">Relative path to the data file, resolved against
		/// the loaded module directory (get_module_directory()). The path
		/// must match the loaded file on later calls.</param>
		/// <returns>Reference to the (immutable) singleton and the STATE describing
		/// whether the data are loaded (OK) or what prevented initialization.</returns>
		static ::std::pair<const Derived&, STATE> try_call(
			const std::filesystem::path& relpath = std::filesystem::path(default_relpath<Derived>)) {
			std::lock_guard lk(mutex_());
			const auto st = ensure_(relpath);
			if (st != STATE::OK)
				return {Empty_(), st};
			return {*storage_(), st};
		}

		/// <summary>
		/// Return the immutable singleton; report load errors and conflicting
		/// paths as exceptions. Failed initialization can be retried.
		/// </summary>
		static const Derived& instance(
			const std::filesystem::path& relpath = std::filesystem::path(default_relpath<Derived>)) {
			auto [object, st] = try_call(relpath);
			if (st != STATE::OK)
				throw std::runtime_error(std::string("Could not load singleton data: ") + to_string(st));
			return object;
		}

		/// <summary>
		/// Current state of the singleton: OK after a successful load,
		/// NOT_INITIALIZED before the first call, or the failure code recorded
		/// during initialization. A failed reload preserves an existing OK state.
		/// Never triggers initialization.
		/// </summary>
		static STATE state() noexcept {
			return state_().load(::std::memory_order_acquire);
		}

		/// <summary>
		/// Check whether the singleton has already been constructed, without
		/// triggering initialization.
		/// </summary>
		static bool exists() noexcept {
			return state_().load(::std::memory_order_acquire) == STATE::OK;
		}

		/// <summary>
		/// Test/diagnostic helper: force the singleton to re-read its data from a
		/// NEW relative location (resolved against the module directory like the
		/// initial path does) and refresh the stored STATE. Returns the resulting
		/// state; on failure the previously loaded data remain untouched.
		/// </summary>
		static STATE reload(const std::filesystem::path& new_relpath) {
			std::lock_guard lk(mutex_());
			return load_(new_relpath);
		}

		/// <summary>
		/// Absolute path of the data file the singleton was created from.
		/// </summary>
		const std::filesystem::path& filename() const noexcept {
			return filename_;
		}

	protected:
		// The derived constructor receives the resolved absolute file path and an
		// output STATE reference it must fill in (OK on success, otherwise the
		// concrete failure code). It must not throw for expected failures.
		explicit SingletonBase(std::filesystem::path absfile) : filename_(::std::move(absfile)) {}

	private:
		// One storage + one state + one mutex per Derived instantiation.
		static ::std::unique_ptr<const Derived>& storage_() {
			static ::std::unique_ptr<const Derived> obj;
			return obj;
		}
		static ::std::atomic<STATE>& state_() {
			static ::std::atomic<STATE> st{STATE::NOT_INITIALIZED};
			return st;
		}
		static ::std::mutex& mutex_() {
			static ::std::mutex mtx;
			return mtx;
		}

		/// <summary>
		/// Fallback object handed out by try_call() when no data could
		/// be loaded (e.g. INVALID_ARGUMENT with an empty path). A valid
		/// reference is kept alive for the whole process, so callers that ignore
		/// the STATE never hit a dangling reference; properly written callers
		/// check the returned state instead. The derived type must expose a
		/// default constructor for this to compile (Distances, ElectronDensity
		/// do — their data members are value-initialized/empty).
		/// </summary>
		struct EmptyHolder {
			const Derived value{};
		};
		static const Derived& Empty_() {
			static const EmptyHolder holder;
			return holder.value;
		}

		// Called under mutex_(). Publish only a fully parsed object; failures
		// leave an existing object intact and permit a later initialization retry.
		static STATE load_(const std::filesystem::path& relpath) {
			STATE st = STATE::ERROR;
			try {
				const auto absfile = get_module_file(relpath);
				if (absfile.empty()) {
					st = STATE::INVALID_ARGUMENT;
				} else {
					std::unique_ptr<const Derived> fresh(new Derived(absfile, st));
					if (st == STATE::OK)
						storage_() = std::move(fresh);
				}
			} catch (...) {
				st = state_from_current_exception();
			}
			if (!storage_())
				state_().store(st, std::memory_order_release);
			else
				state_().store(STATE::OK, std::memory_order_release);
			return st;
		}

		static STATE ensure_(const std::filesystem::path& relpath) {
			if (!storage_())
				return load_(relpath);
			try {
				return get_module_file(relpath) == storage_()->filename()
					? STATE::OK : STATE::INVALID_ARGUMENT;
			} catch (...) {
				return state_from_current_exception();
			}
		}

		std::filesystem::path filename_;
	};
}
