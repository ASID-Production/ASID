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

cmake_minimum_required (VERSION 3.12)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(cpplibLib)

  add_library(cpplibLib "cpplib/src/Functions/Functions.cpp" 
                               "cpplib/src/BaseHeaders/Support.h"
                               "cpplib/src/Classes/Distances.h"
                               "cpplib/src/Classes/FindMolecules.h"
                               "cpplib/src/Classes/Hash.h"
                               "cpplib/src/Classes/Interfaces.h"
                               "cpplib/src/Classes/MoleculeGraph.h"
                               "cpplib/src/Classes/Engine.h"
                               "cpplib/src/Classes/Geometry.h"
                               "cpplib/src/Classes/SearchGraph.h"
                               "cpplib/src/Platform/Definitions.h"
                               "cpplib/src/Functions/AllInOneAndCurrent.h")

target_include_directories(cpplibLib
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/cpplib/src/BaseHeaders
        ${PROJECT_SOURCE_DIR}/cpplib/src/Classes
        ${PROJECT_SOURCE_DIR}/cpplib/src/Platform
)

project(cpplib)

############################################################
# Create a library
############################################################


#Generate the shared library from the library sources
if(CMAKE_BUILD_TYPE MATCHES Release)
  add_library(cpplib SHARED "cpplib/src/Functions/Functions.cpp" 
                               "cpplib/src/BaseHeaders/Support.h"
                               "cpplib/src/Classes/Distances.h"
                               "cpplib/src/Classes/FindMolecules.h"
                               "cpplib/src/Classes/Hash.h"
                               "cpplib/src/Classes/Interfaces.h"
                               "cpplib/src/Classes/MoleculeGraph.h"
                               "cpplib/src/Classes/Engine.h"
                               "cpplib/src/Classes/Geometry.h"
                               "cpplib/src/Classes/SearchGraph.h"
                               "cpplib/src/Platform/Definitions.h"
                               "cpplib/src/Functions/AllInOneAndCurrent.h")

target_include_directories(cpplib
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/cpplib/src/BaseHeaders
        ${PROJECT_SOURCE_DIR}/cpplib/src/Classes
        ${PROJECT_SOURCE_DIR}/cpplib/src/Platform
)

add_custom_command(
        TARGET cpplib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/cpplib/BondLength.ini
                ${CMAKE_CURRENT_BINARY_DIR}/cpplib/Tests/modules/c_modules/BondLength.ini)
add_custom_command(
    TARGET cpplib POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/cpplib/BondLength.ini
            ${CMAKE_SOURCE_DIR}/api_database/django_project/modules/c_modules/BondLength.ini)
if (UNIX)
    add_custom_command(
        TARGET cpplib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_FILE:cpplib>
                ${CMAKE_SOURCE_DIR}/api_database/django_project/modules/c_modules/cpplib.so)
endif (UNIX)

if (WIN32)
    add_custom_command(
        TARGET cpplib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_FILE:cpplib>
                ${CMAKE_SOURCE_DIR}/api_database/django_project/modules/c_modules/cpplib.dll)
endif (WIN32)

endif(CMAKE_BUILD_TYPE MATCHES Release)


 add_subdirectory ("cpplib/Tests")
 add_subdirectory ("VnE/Source/Extensions/ChemPackSource")
