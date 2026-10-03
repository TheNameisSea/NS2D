# Include the FetchContent module
include(FetchContent)

# Declare Catch2 dependency
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.16.0 
  GIT_SHALLOW TRUE
  )

# Make Catch2 available to the build graph
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

add_executable(test_grid test_grid.cpp )
target_link_libraries(test_grid PRIVATE ns2d Catch2::Catch2WithMain ns2d_warnings)

add_executable(test_field test_field.cpp )
target_link_libraries(test_field PRIVATE ns2d Catch2::Catch2WithMain ns2d_warnings)

add_executable(test_bc test_bc.cpp )
target_link_libraries(test_bc PRIVATE ns2d Catch2::Catch2WithMain ns2d_warnings)

catch_discover_tests(test_grid)
catch_discover_tests(test_field)
catch_discover_tests(test_bc)