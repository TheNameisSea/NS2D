# Include the FetchContent module
include(FetchContent)

# Declare Catch2 dependency
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.16.0 # Replace with the latest stable version if desired
  GIT_SHALLOW TRUE
  )

# Make Catch2 available to the build graph
FetchContent_MakeAvailable(Catch2)

message(STATUS "tests/ included")