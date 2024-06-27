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

project(MyTests)

cmake_minimum_required(VERSION 3.14)

# GoogleTest requires at least C++14

if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif(COMPILER_SUPPORTS_CXX14)

if (POLICY CMP0135)
cmake_policy(SET CMP0135 OLD) 
endif()
include(FetchContent)
FetchContent_Declare(
  gest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gest)

enable_testing()


add_executable(
    ${PROJECT_NAME} 
    MainTest.cpp
    Includes.cpp 
    NodeTest.cpp)
target_link_libraries(
    ${PROJECT_NAME} 
    GTest::gtest_main
    cpplibLib
)
set_property(TARGET MyTests PROPERTY CXX_STANDARD 14)

add_custom_command(
        TARGET MyTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/module/BondLength.ini
                ${CMAKE_CURRENT_BINARY_DIR}/BondLength.ini)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})


