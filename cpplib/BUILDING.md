# Building cpplib

Configure the repository root with CMake and build with a C++20 compiler and Python development libraries.
Use `-DBUILD_TESTING=ON` to register GoogleTest tests on all platforms.
PythonTests is registered only when pytest is installed in the selected Python interpreter **at configuration time**.
After installing pytest, re-run CMake to register PythonTests. CMake emits a warning when pytest is missing.
Python tests use the extension produced by the build through `CPPLIB_TEST_EXTENSION`; set `-DCPPLIB_COPY_MODULE_TO_SOURCE=OFF` to keep test builds from overwriting the source module.

Release builds preserve IEEE floating-point semantics by default. `CPPLIB_NATIVE_CPU=ON` opts into local CPU instructions
(MSVC uses AVX2); leave it OFF for distributed binaries. `CPPLIB_ENABLE_IPO=ON` enables Release IPO through CMake's
compiler/linker support check on all platforms. Both options default to OFF.

Run `ctest --test-dir <build-dir> --output-on-failure` after building.

On Windows the extension uses Python Stable ABI (`python3.lib`) with the selected release interpreter in both C++ Debug and Release configurations. C++ assertions remain enabled in Debug.
