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

    class DataLoadError : public std::runtime_error {
    public:
        DataLoadError(STATE code, const std::filesystem::path& file, const std::string& detail = {})
            : std::runtime_error(std::string("Could not load singleton data: ") + to_string(code) +
                " [" + path_text(file) + "]" + (detail.empty() ? "" : ": " + detail)), code_(code) {}
        STATE code() const noexcept { return code_; }
    private:
        static std::string path_text(const std::filesystem::path& file) {
            const auto text = file.u8string();
            return {reinterpret_cast<const char*>(text.data()), text.size()};
        }
        STATE code_;
    };

	inline constexpr bool succeeded(STATE s) noexcept { return s == STATE::OK; }

	/// <summary>
	/// Translates a standard exception thrown while reading/parsing a data file
	/// into a program-wide STATE code (used instead of propagating exceptions out
	/// of the singleton interface).
	/// </summary>
	[[nodiscard]] inline STATE state_from_exception(const std::exception& e) noexcept {
		if (const auto* error = dynamic_cast<const DataLoadError*>(&e)) return error->code();
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

    // Immutable, lazily loaded data. Production access throws on
    // critical loading failures. try_call is a compatibility/status adapter.
    // Published once per process lifetime. Initialized access is one acquire
    // load, and the immutable object lives until process shutdown.
    template <typename Derived>
    class SingletonBase {
    public:
        SingletonBase() = delete;
        SingletonBase(const SingletonBase&) = delete;
        SingletonBase(SingletonBase&&) = delete;
        SingletonBase& operator=(const SingletonBase&) = delete;
        SingletonBase& operator=(SingletonBase&&) = delete;

        static const Derived& instance() {
            if (const auto* object = published_().load(std::memory_order_acquire)) return *object;
            std::lock_guard lock(mutex_());
            if (!storage_()) load_(std::filesystem::path(default_relpath<Derived>));
            return *storage_();
        }
        static const Derived& instance(const std::filesystem::path& relpath) {
            std::lock_guard lock(mutex_());
            if (!storage_()) load_(relpath);
            else {
                const auto request = normalized_(relpath);
                if (relpath.empty() || (request != request_key_() &&
                    (module_key_() / request).lexically_normal() != storage_()->filename()))
                    throw DataLoadError(STATE::INVALID_ARGUMENT, relpath);
            }
            return *storage_();
        }
        static std::pair<const Derived&, STATE> try_call() {
            try { return {instance(), STATE::OK}; }
            catch (...) { return {Empty_(), state_from_current_exception()}; }
        }
        static std::pair<const Derived&, STATE> try_call(const std::filesystem::path& relpath) {
            try { return {instance(relpath), STATE::OK}; }
            catch (...) { return {Empty_(), state_from_current_exception()}; }
        }
        static STATE state() noexcept { return state_().load(std::memory_order_acquire); }
        static bool exists() noexcept { return published_().load(std::memory_order_acquire) != nullptr; }


		/// <summary>
		/// Absolute path of the data file the singleton was created from.
		/// </summary>
		const std::filesystem::path& filename() const noexcept {
			return filename_;
		}

	protected:
		// Derived loads the resolved absolute file, throwing on critical failures.
		explicit SingletonBase(std::filesystem::path absfile) : filename_(::std::move(absfile)) {}

	private:
		// One storage + one state + one mutex per Derived instantiation.
		static ::std::unique_ptr<const Derived>& storage_() {
			static ::std::unique_ptr<const Derived> obj;
			return obj;
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

        static std::atomic<const Derived*>& published_() {
            static std::atomic<const Derived*> object{nullptr};
            return object;
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

        // Called under mutex_(), only until the first successful load. Failures
        // leave no object published and can be retried on the next call.
        static void load_(const std::filesystem::path& relpath) {
            try {
                if (relpath.empty()) throw DataLoadError(STATE::INVALID_ARGUMENT, relpath);
                auto directory = get_module_directory();
                auto request = normalized_(relpath);
                const auto absfile = std::filesystem::weakly_canonical(directory / request);
                std::unique_ptr<const Derived> fresh(new Derived(absfile));
                module_key_().swap(directory);
                request_key_().swap(request);
                storage_() = std::move(fresh);
                state_().store(STATE::OK, std::memory_order_release);
                published_().store(storage_().get(), std::memory_order_release);
            } catch (const DataLoadError& error) {
                if (!published_().load(std::memory_order_relaxed)) state_().store(error.code(), std::memory_order_release);
                throw;
            } catch (const std::exception& error) {
                const auto code = state_from_exception(error);
                if (!published_().load(std::memory_order_relaxed)) state_().store(code, std::memory_order_release);
                throw DataLoadError(code, relpath, error.what());
            }
        }

		std::filesystem::path filename_;
	};
}
