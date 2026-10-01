#include "ns2d/ns2d.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("simple test", "[sanity]"){
    REQUIRE(checkSrc2d() == 1);
}