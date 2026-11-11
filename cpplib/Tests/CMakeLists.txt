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
cmake_minimum_required(VERSION 3.14)

project(MyTests LANGUAGES CXX)

# ============================================================================
# Google Test
# ============================================================================

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

# For Windows: Prevent overriding parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# ============================================================================
# Test sources
# ============================================================================

file(GLOB_RECURSE TEST_SRC CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(${PROJECT_NAME}
    ${TEST_SRC}
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/Functions/Functions.cpp"
)

# ============================================================================
# Include directories
# ============================================================================

target_include_directories(${PROJECT_NAME}
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/BaseHeaders"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/Classes"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/Functions"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/Platform"
)

# ============================================================================
# C++ standard
# ============================================================================

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

# ============================================================================
# Linking
# ============================================================================

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        GTest::gtest_main
        optimizations
        ${CMAKE_DL_LIBS}
)

# ============================================================================
# Copy DATA files
# ============================================================================

if(DATA_FULL_PATHS)
    add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${DATA_FULL_PATHS}
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
        COMMENT "Copying DATA files into test directory"
        VERBATIM
    )
endif()

# ============================================================================
# Google Test discovery
# ============================================================================

include(GoogleTest)

gtest_discover_tests(${PROJECT_NAME}
    WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
    DISCOVERY_MODE POST_BUILD
)