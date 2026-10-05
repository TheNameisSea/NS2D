cmake_minimum_required(VERSION 3.20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(NavierStokes2D
        VERSION 0.1.0
        DESCRIPTION "A simple 2D Navier-Stokes solver"
        HOMEPAGE_URL "github.com/TheNameisSea/cppNavierStokesSolver"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

FetchContent_Declare(json 
                    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
                    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
                    SYSTEM)
FetchContent_MakeAvailable(json)

# Add the source files for the library
add_library(ns2d STATIC src/grid.cpp src/field.cpp src/boundary.cpp 
                        src/operators.cpp src/projection.cpp 
                        src/simulation.cpp src/poisson.cpp)
target_include_directories(ns2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_library(ns2d_warnings INTERFACE)
target_compile_options(ns2d_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion)
target_link_libraries(ns2d PRIVATE ns2d_warnings)

add_library(ns2d_io STATIC src/io/config.cpp src/io/vtk_writer.cpp)
target_link_libraries(ns2d_io PUBLIC ns2d)
target_link_libraries(ns2d_io PUBLIC nlohmann_json::nlohmann_json)
target_link_libraries(ns2d_io PRIVATE ns2d_warnings)

add_executable(cavity apps/cavity.cpp)
target_link_libraries(cavity PRIVATE ns2d_io ns2d_warnings)


# Test section
option(NS2D_BUILD_TESTS "Build unit tests" ON)
if(NS2D_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

