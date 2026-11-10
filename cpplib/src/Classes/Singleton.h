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

namespace cpplib {
	template <typename T>
	inline constexpr std::string_view default_relpath{};

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
		ERROR_STATE,              ///< Any other unspecified failure.
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
			case STATE::ERROR_STATE:           return "ERROR_STATE";
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
        if (what.find("read failure") != std::string_view::npos) return STATE::IO_ERROR;
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
		    what.find("incorrect data") != std::string_view::npos ||
            dynamic_cast<const std::invalid_argument*>(&e)) {
			return STATE::BAD_CONTENT;
		}
		return STATE::ERROR_STATE;
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
			return STATE::ERROR_STATE;
		}
	}

	struct ModuleDirectoryState {
		std::mutex mutex;
		std::filesystem::path directory;
		bool initialized = false;
	};

	inline ModuleDirectoryState& module_directory_state() {
		static ModuleDirectoryState state;
		return state;
	}

	/// <summary>
	/// Use the directory reported by Python for its imported extension module.
	/// It must be configured before singleton paths are first resolved.
	/// </summary>
	inline void set_module_directory(std::filesystem::path directory) {
		if (directory.empty())
			throw std::invalid_argument("Python module directory is empty");
		directory = std::filesystem::weakly_canonical(std::filesystem::absolute(directory));

		auto& state = module_directory_state();
		std::lock_guard lock(state.mutex);
		if (state.initialized) {
			if (state.directory != directory)
				throw std::logic_error("cpplib module directory is already in use");
			return;
		}
		state.directory = std::move(directory);
		state.initialized = true;
	}

	/// <summary>
	/// Return the configured Python module directory, or an uncached current
	/// working directory for standalone C++ use.
	/// </summary>
	inline std::filesystem::path get_module_directory() {
		auto& state = module_directory_state();
		std::lock_guard lock(state.mutex);
		return state.initialized ? state.directory : std::filesystem::current_path();
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
			const std::filesystem::path& relpath) {
			std::lock_guard lk(mutex_());
			const auto st = ensure_(relpath);
			if (st != STATE::OK)
				return {Empty_(), st};
			return {*storage_(), st};
		}

        // An omitted path means the current data, independent of cwd changes.
        static std::pair<const Derived&, STATE> try_call() {
            std::lock_guard lk(mutex_());
            if (storage_()) return {*storage_(), STATE::OK};
            const auto st = load_(std::filesystem::path(default_relpath<Derived>));
            return {st == STATE::OK ? *storage_() : Empty_(), st};
        }
        static const Derived& instance() {
            auto [object, st] = try_call();
            if (st != STATE::OK)
                throw std::runtime_error(std::string("Could not load singleton data: ") + to_string(st));
            return object;
        }

		/// <summary>
		/// Return the immutable singleton; report load errors and conflicting
		/// paths as exceptions. Failed initialization can be retried.
		/// </summary>
		static const Derived& instance(
			const std::filesystem::path& relpath) {
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
        // Raw references escape the API: no bounded safe reclamation period exists.
        // Keep previous immutable generations alive until process shutdown.
        static std::vector<std::unique_ptr<const Derived>>& retired_() {
            static std::vector<std::unique_ptr<const Derived>> objects;
            return objects;
        }
        static std::filesystem::path& module_key_() {
            static std::filesystem::path directory;
            return directory;
        }
        static std::filesystem::path& request_key_() {
            static std::filesystem::path request;
            return request;
        }
        static std::filesystem::path normalized_(const std::filesystem::path& path) {
            auto local = path.lexically_normal();
            if (!local.empty() && *local.begin() == "module")
                local = local.lexically_relative("module");
            return local;
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
            if (relpath.empty()) {
                if (!storage_()) state_().store(STATE::INVALID_ARGUMENT, std::memory_order_release);
                return STATE::INVALID_ARGUMENT;
            }
			STATE st = STATE::ERROR_STATE;
			try {
				auto directory = get_module_directory();
                auto request = normalized_(relpath);
                const auto absfile = relpath.empty() ? std::filesystem::path{} :
                    std::filesystem::weakly_canonical(directory / request);
				if (absfile.empty()) {
					st = STATE::INVALID_ARGUMENT;
				} else {
					std::unique_ptr<const Derived> fresh(new Derived(absfile, st));
					if (st == STATE::OK) {
                        if (storage_()) retired_().push_back(std::move(storage_()));
                        module_key_().swap(directory);
                        request_key_().swap(request);
                        storage_() = std::move(fresh);
                    }
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
				if (relpath.empty()) return STATE::INVALID_ARGUMENT;
                const auto request = normalized_(relpath);
                return (request == request_key_() ||
                    (module_key_() / request).lexically_normal() == storage_()->filename())
					? STATE::OK : STATE::INVALID_ARGUMENT;
			} catch (...) {
				return state_from_current_exception();
			}
		}

		std::filesystem::path filename_;
	};
}
