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

add_executable(test_sanity test_sanity.cpp)
target_link_libraries(test_sanity PRIVATE ns2d Catch2::Catch2WithMain ns2d_warnings)

catch_discover_tests(test_sanity)