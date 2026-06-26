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

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# Politics
# ============================================================================
if(CMAKE_VERSION VERSION_LESS 3.15)
    message(STATUS "CMake version 3.15+ recomended")
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    cmake_policy(SET CMP0141 NEW)
    set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT 
        "$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>")
endif()

# ============================================================================
# Optimizations (interface)
# ============================================================================
set(CMAKE_SHOW_INCLUDES OFF)
add_library(optimizations INTERFACE)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_compile_options(optimizations INTERFACE 
            $<$<CONFIG:Release>:
                -O3
                -DNDEBUG
                -march=native
                -mtune=native
                -mfma
                -flto
                -funroll-loops
                -fprefetch-loop-arrays
                -fno-exceptions
                -D_HAS_EXCEPTIONS=0
                -fno-rtti
                -ffast-math
            >
        )
        message(STATUS "Linux: max performance with safe math optimizations")
        
    elseif(MSVC)
        foreach(flag_var 
                CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_RELWITHDEBINFO)
            string(REPLACE "/GR" "" ${flag_var} "${${flag_var}}")
            string(REPLACE "/EHс" "" ${flag_var} "${${flag_var}}")
            string(REPLACE "/EHsс" "" ${flag_var} "${${flag_var}}")
        endforeach()


        target_compile_options(optimizations INTERFACE 
            $<$<CONFIG:Release>:
                /O2
                /Ot
                /Oi
                /Gy
                /GL
                /DNDEBUG
                /arch:AVX2
                /fp:fast
                /D_HAS_EXCEPTIONS=0
                /EHс-
                /EHsс-
                /EHs-
                /GR-
            >
        )
        target_link_options(optimizations INTERFACE 
            $<$<CONFIG:Release>:
                /LTCG
                /OPT:REF
                /OPT:ICF
            >
        )
        message(STATUS "Windows: AVX2 + LTCG with fast FP model")
        
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
        target_compile_options(optimizations INTERFACE 
            $<$<CONFIG:Release>:
                -O3
                -DNDEBUG
                -march=native
                -mtune=native
                -flto
                -fno-math-errno
                -fno-signed-zeros
                -fno-trapping-math
            >
        )
        message(STATUS "Clang/Unix: safe optimizations enabled")
    endif()

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
# OBJECT library
# ============================================================================
add_library(cpplib_objects OBJECT ${CPPLIB_SOURCES})

target_include_directories(cpplib_objects 
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/cpplib/src/BaseHeaders
        ${PROJECT_SOURCE_DIR}/cpplib/src/Classes
        ${PROJECT_SOURCE_DIR}/cpplib/src/Platform
)
set_target_properties(cpplib_objects PROPERTIES
    PREFIX ""
    OUTPUT_NAME "cpplib_objects"
    LINKER_LANGUAGE C
    POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(cpplib_objects PRIVATE optimizations)



# ============================================================================
# CPPLIB creates only in Release
# ============================================================================
if(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_CONFIGURATION_TYPES)
    add_library(cpplib SHARED ${CPPLIB_PY_SOURCES})
    
    # Include dirs
    target_include_directories(cpplib 
        PUBLIC 
            ${PROJECT_SOURCE_DIR}/cpplib/src/BaseHeaders
            ${PROJECT_SOURCE_DIR}/cpplib/src/Classes
            ${PROJECT_SOURCE_DIR}/cpplib/src/Platform
        PRIVATE
            ${Python_INCLUDE_DIRS}
    )
    
    set_target_properties(cpplib PROPERTIES
        PREFIX ""
        OUTPUT_NAME "cpplib"
        LINKER_LANGUAGE C
        POSITION_INDEPENDENT_CODE ON
    )
    
    target_link_libraries(cpplib PRIVATE cpplib_objects optimizations)
    # ============================================================================
    # Link with Python (Stable ABI / Limited API)
    # ============================================================================
    
    if(WIN32)
        find_library(PYTHON_SABI_LIBRARY_RELEASE
            NAMES python3
            PATHS ${Python_LIBRARY_DIRS}
            REQUIRED
        )
        
        find_library(PYTHON_SABI_LIBRARY_DEBUG
            NAMES python3_d
            PATHS ${Python_LIBRARY_DIRS}
            REQUIRED
        )
        
        target_link_libraries(cpplib PRIVATE
            optimized ${PYTHON_SABI_LIBRARY_RELEASE}
            debug ${PYTHON_SABI_LIBRARY_DEBUG}
        )
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
    
    add_custom_command(
        TARGET cpplib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:cpplib>
            ${CMAKE_SOURCE_DIR}/module/cpplib.${MODULE_EXTENSION}
        COMMENT "Copying module to ${CMAKE_SOURCE_DIR}/module/"
    )    
endif()

# ============================================================================
# Tests
# ============================================================================
if(EXISTS ${CMAKE_SOURCE_DIR}/cpplib/Tests AND WIN32)
    add_subdirectory(cpplib/Tests)
    message(STATUS "Tests enabled for Windows build")
elseif(EXISTS ${CMAKE_SOURCE_DIR}/cpplib/Tests)
    message(STATUS "Tests directory exists but disabled (requires Windows)")
endif()

# add_compile_definitions(DEBUG_MESS_ON)