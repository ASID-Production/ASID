# Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
# This file is part of ASID - Atomistic Simulation Instruments and Database
# For more information see <https://github.com/ASID-Production/ASID>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ******************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# ******************************************************************************************

cmake_minimum_required(VERSION 3.12)

# ============================================================================
# Project and basic options
# ============================================================================
project(cpplib LANGUAGES CXX C)

if(MSVC AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
    set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add_compile_definitions(DEBUG_MESS_ON)

# ============================================================================
# Politics
# ============================================================================
if(CMAKE_VERSION VERSION_LESS 3.15)
    message(STATUS "CMake version 3.15+ recommended")
endif()
if(POLICY CMP0141)
    cmake_policy(SET CMP0141 NEW)
    set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
        "$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>")
endif()

# ============================================================================
# Optimizations (interface)
# ============================================================================
set(CMAKE_SHOW_INCLUDES OFF)
add_library(optimizations INTERFACE)
option(CPPLIB_NATIVE_CPU "Optimize Release builds for the local CPU" OFF)
option(CPPLIB_ENABLE_IPO "Enable interprocedural optimization when supported" OFF)
if(MSVC)
    target_compile_options(optimizations INTERFACE /fp:precise
        $<$<CONFIG:Release>:/O2>)
    if(CPPLIB_NATIVE_CPU)
        target_compile_options(optimizations INTERFACE $<$<CONFIG:Release>:/arch:AVX2>)
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(optimizations INTERFACE $<$<CONFIG:Release>:-O3>)
    if(CPPLIB_NATIVE_CPU)
        target_compile_options(optimizations INTERFACE
            "$<$<CONFIG:Release>:-march=native;-mtune=native>")
    endif()
endif()
if(CPPLIB_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT CPPLIB_IPO_SUPPORTED OUTPUT CPPLIB_IPO_ERROR LANGUAGES CXX)
    if(NOT CPPLIB_IPO_SUPPORTED)
        message(FATAL_ERROR "IPO requested but unavailable: ${CPPLIB_IPO_ERROR}")
    endif()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()
message(STATUS "cpplib: precise floating point; native CPU=${CPPLIB_NATIVE_CPU}; IPO=${CPPLIB_ENABLE_IPO}")

# ============================================================================
# Python options
# ============================================================================
if(WIN32)
    find_package(Python COMPONENTS Interpreter Development EXACT 3.10 REQUIRED)
else()
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python EXECUTABLE: ${Python_EXECUTABLE}")

# ============================================================================
# Files
# ============================================================================
file(GLOB_RECURSE CPPLIB_HEADERS
    "cpplib/src/*.h"
)

set(CPPLIB_SOURCES
    cpplib/src/Functions/Functions.cpp
)
set(CPPLIB_PY_SOURCES
    cpplib/src/Python/PythonInterface.cpp
)

set(DATA_FILES
    BondLength.ini
    ED.txt
)
set(DATA_FULL_PATHS "")
foreach(file ${DATA_FILES})
    list(APPEND DATA_FULL_PATHS ${CMAKE_SOURCE_DIR}/module/${file})
endforeach()

# ============================================================================
# CPPLIB extension (all configurations)
# ============================================================================
add_library(cpplib SHARED ${CPPLIB_PY_SOURCES} ${CPPLIB_SOURCES})

# Include dirs
target_include_directories(cpplib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/cpplib/src/BaseHeaders
        ${PROJECT_SOURCE_DIR}/cpplib/src/Classes
        ${PROJECT_SOURCE_DIR}/cpplib/src/Functions
        ${PROJECT_SOURCE_DIR}/cpplib/src/Platform
    PRIVATE
        ${Python_INCLUDE_DIRS}
)

set_target_properties(cpplib PROPERTIES
    PREFIX ""
    OUTPUT_NAME "cpplib"
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(cpplib PRIVATE optimizations ${CMAKE_DL_LIBS})

if(MSVC)
    target_compile_options(cpplib PRIVATE $<$<CONFIG:Release>:/MD>)
endif()

# ============================================================================
# Link with Python (Stable ABI / Limited API)
# ============================================================================

if(WIN32)
    find_library(PYTHON_SABI_LIBRARY_RELEASE
        NAMES python3
        PATHS ${Python_LIBRARY_DIRS}
        REQUIRED
    )

    # Stable ABI targets the selected release interpreter in every C++ config.
    target_link_libraries(cpplib PRIVATE ${PYTHON_SABI_LIBRARY_RELEASE})
else()
    target_link_libraries(cpplib PRIVATE ${Python_LIBRARIES})
endif()

# ============================================================================
# Copy target to module dir
# ============================================================================
if(WIN32)
    set(MODULE_EXTENSION "pyd")
else()
    set(MODULE_EXTENSION "so")
endif()

set_target_properties(cpplib PROPERTIES SUFFIX ".${MODULE_EXTENSION}")
add_custom_command(TARGET cpplib POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${DATA_FULL_PATHS} "$<TARGET_FILE_DIR:cpplib>"
    VERBATIM)
option(CPPLIB_COPY_MODULE_TO_SOURCE "Copy the extension into module/ after building" ON)
if(CPPLIB_COPY_MODULE_TO_SOURCE)
    add_custom_command(TARGET cpplib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:cpplib> "${CMAKE_SOURCE_DIR}/module/cpplib.${MODULE_EXTENSION}"
        VERBATIM)
endif()


# ============================================================================
# Tests
# ============================================================================


include(CTest)

if(BUILD_TESTING)

    # ------------------------------------------------------------------------
    # C++ tests (Google Test)
    # ------------------------------------------------------------------------
    if(EXISTS "${CMAKE_SOURCE_DIR}/cpplib/Tests")
        add_subdirectory(cpplib/Tests)
        message(STATUS "Google Test enabled")
    endif()

    # ------------------------------------------------------------------------
    # Python tests (pytest)
    # ------------------------------------------------------------------------
    set(PYTHON_TEST_DIR "${CMAKE_SOURCE_DIR}/cpplib/Tests/Python")

    if(EXISTS "${PYTHON_TEST_DIR}" AND TARGET cpplib)

        # pytest availability is checked at configure time. Re-run CMake after installing it.
        execute_process(
            COMMAND ${Python_EXECUTABLE} -m pytest --version
            RESULT_VARIABLE PYTEST_RESULT
            OUTPUT_QUIET
            ERROR_QUIET
        )

        if(PYTEST_RESULT EQUAL 0)

            add_test(
                NAME PythonTests
                COMMAND
                    ${CMAKE_COMMAND} -E env "CPPLIB_TEST_EXTENSION=$<TARGET_FILE:cpplib>"
                    ${Python_EXECUTABLE}
                    -m pytest
                    -v
                    ${PYTHON_TEST_DIR}
            )

            # Allows Python tests to import cpplib from module/
            set_tests_properties(PythonTests PROPERTIES
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/module"
            )

            message(STATUS "Python pytest enabled")

        else()
            message(WARNING
                "pytest not found. Install it using: "
                "${Python_EXECUTABLE} -m pip install pytest, then re-run CMake to register PythonTests."
            )
        endif()

    endif()

endif()

# add_compile_definitions(DEBUG_MESS_ON)
