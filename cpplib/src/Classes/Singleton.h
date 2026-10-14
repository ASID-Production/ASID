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
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

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
	/// Platform-independent location of the starting module (executable / shared
	/// library) that hosts cpplib code, e.g. .../site-packages/cpplib.
	/// Uses the C++20 standard reflection std::source_location::current(): it yields
	/// the path of this very header at compile time - no WinAPI, dladdr or
	/// compiler-specific intrinsics involved. Data files shipped next to the module
	/// (module/BondLength.ini, module/ED.txt) are resolved relative to it, not to
	/// the process CWD. If the location cannot be determined, falls back to the
	/// current working directory ("."). Computed once on first use (thread-safe
	/// magic static).
	/// </summary>
	inline const std::filesystem::path& get_module_directory() {
		static const std::filesystem::path dir = [] {
			const auto refl = std::source_location::current();
			const std::filesystem::path self(refl.file_name());
			// Absolute path to this very header proves the module location is known.
			if (self.is_absolute()) {
				// <module_dir>/src/Classes/Singleton.h -> <module_dir>
				return self.parent_path().parent_path().parent_path();
			}
			return std::filesystem::path(".");
		}();
		return dir;
	}

	/// <summary>
	/// Resolve a relative data-file path against the starting module directory,
	/// exactly the way SingletonBase::instance() does.
	/// </summary>
	inline std::filesystem::path get_module_file(const std::filesystem::path& relpath) {
		return get_module_directory() / relpath;
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
	///       subsequent call IGNORES the argument and returns the previously
	///       loaded version together with its stored state. If loading fails, the
	///       state is FILE_NOT_FOUND / BAD_SIZE / BAD_CONTENT / IO_ERROR instead
	///       of an exception being thrown. Thread-safe (mutex-guarded).
	///   - Derived::instance(relpath = ...) -> const Derived&
	///       Convenience accessor: same semantics, returns the reference only;
	///       when the load failed the returned object simply holds no data
	///       (state() reports why). Use try_call() when the status is needed.
	///   - Derived::state()
	///       Current STATE of the singleton: OK after a successful load,
	///       NOT_INITIALIZED before the very first call ("несобран"), or the
	///       concrete failure code recorded during initialization.
	///   - Derived::exists()
	///       true if the instance has already been constructed; never triggers
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
		/// the starting module directory (get_module_directory()). Significant
		/// only for the very first call.</param>
		/// <returns>Reference to the (immutable) singleton and the STATE describing
		/// whether the data are loaded (OK) or what prevented initialization.</returns>
		static ::std::pair<const Derived&, STATE> try_call(
			const std::filesystem::path& relpath = std::filesystem::path(default_relpath<Derived>)) {
			ensure_(relpath);
			if (!storage_().has_value()) // INVALID_ARGUMENT: nothing could be loaded
				return {Empty_<>(), state_()};
			return {*storage_(), state_()};
		}

		/// <summary>
		/// Convenience wrapper over try_call(): returns the reference to the
		/// immutable singleton. The argument matters only on the first call;
		/// every later call ignores it and yields the previously loaded version.
		/// Check state() (or use try_call()) if the load status matters.
		/// </summary>
		static const Derived& instance(
			const std::filesystem::path& relpath = std::filesystem::path(default_relpath<Derived>)) {
			ensure_(relpath);
			if (!storage_().has_value())
				return Empty_<>();
			return *storage_();
		}

		/// <summary>
		/// Current state of the singleton: OK after a successful load,
		/// NOT_INITIALIZED before the first call, or the failure code recorded
		/// during the last (re)load attempt. Never triggers initialization.
		/// </summary>
		static STATE state() noexcept {
			return state_().load(::std::memory_order_acquire);
		}

		/// <summary>
		/// Check whether the singleton has already been constructed, without
		/// triggering initialization.
		/// </summary>
		static bool exists() noexcept {
			return state_().load(::std::memory_order_acquire) != STATE::NOT_INITIALIZED;
		}

		/// <summary>
		/// Test/diagnostic helper: force the singleton to re-read its data from a
		/// NEW relative location (resolved against the module directory like the
		/// initial path does) and refresh the stored STATE. Returns the resulting
		/// state; on failure the previously loaded data remain untouched.
		/// </summary>
		static STATE reload(const std::filesystem::path& new_relpath) {
			const std::filesystem::path absfile = get_module_file(new_relpath);
			if (absfile.empty())
				return STATE::INVALID_ARGUMENT;
			// Scratch state filled in by the fresh object's constructor. It is
			// copied into the shared atomic while still alive, because the derived
			// ctor keeps a reference to it only for the duration of construction.
			STATE st = STATE::ERROR;
			::std::lock_guard lk(mutex_());
			try {
				storage_().emplace(absfile, ::std::ref(st));
				storage_()->swap_state_target(state_());
			} catch (...) {
				// Failed reload must not destroy the previously loaded data:
				// ensure_() recreates the old instance from the default location
				// if emplace() left the storage empty.
				st = state_from_current_exception();
				ensure_(default_relpath<Derived>);
				return state_().load(::std::memory_order_acquire);
			}
			if (st == STATE::NOT_INITIALIZED)
				st = STATE::OK; // defensive: ctor forgot to report
			state_().store(st, ::std::memory_order_release);
			return st;
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
		static ::std::optional<Derived>& storage_() {
			static ::std::optional<Derived> obj;
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
		/// Fallback object handed out by instance()/try_call() when no data could
		/// ever be loaded (e.g. INVALID_ARGUMENT with an empty path). A valid
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

		static void ensure_(const std::filesystem::path& relpath) {
			// Fast path: already initialized (success OR recorded failure) —
			// subsequent calls ignore the argument and reuse the loaded version.
			if (state_().load(::std::memory_order_acquire) != STATE::NOT_INITIALIZED)
				return;
			::std::lock_guard lk(mutex_());
			if (state_().load(::std::memory_order_relaxed) != STATE::NOT_INITIALIZED)
				return; // another thread won the race
			STATE st = STATE::ERROR;
			const std::filesystem::path absfile = relpath.empty() ? std::filesystem::path{}
			                                                      : get_module_file(relpath);
			if (absfile.empty()) {
				st = STATE::INVALID_ARGUMENT;
			} else {
				// Scratch state filled in by the constructor; copied into the
				// shared atomic while still alive (the derived ctor keeps only a
				// reference to it for the duration of construction).
				try {
					storage_().emplace(absfile, ::std::ref(st));
					// Rebind the object's internal reference from the dead local
					// scratch to the long-lived shared atomic.
					storage_()->swap_state_target(state_());
					if (st == STATE::NOT_INITIALIZED)
						st = STATE::OK; // defensive: ctor forgot to report
				} catch (...) {
					st = state_from_current_exception();
					storage_().reset();
				}
			}
			state_().store(st, ::std::memory_order_release);
		}

		std::filesystem::path filename_;
	};
}
