# Catch2 (test framework)
include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.16.0
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# ns2d_add_test(<name> <library>)
# Builds tests/<name>.cpp into an executable linked to <library>, Catch2 and our
# warnings, and registers every TEST_CASE in it with ctest.
function(ns2d_add_test name library)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ${library} Catch2::Catch2WithMain ns2d_warnings)
  catch_discover_tests(${name})
endfunction()

# Physics library tests
ns2d_add_test(test_grid       ns2d)
ns2d_add_test(test_field      ns2d)
ns2d_add_test(test_bc         ns2d)
ns2d_add_test(test_operators  ns2d)
ns2d_add_test(test_projection ns2d)
ns2d_add_test(test_poisson    ns2d)
ns2d_add_test(test_simulation ns2d)

# I/O library tests
ns2d_add_test(test_config ns2d_io)
target_compile_definitions(test_config PRIVATE NS2D_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
ns2d_add_test(test_vtk     ns2d_io)
ns2d_add_test(test_summary ns2d_io)
