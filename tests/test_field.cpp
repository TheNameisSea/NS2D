#include "ns2d/grid.h"
#include "ns2d/field.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

TEST_CASE("field test", "[field]"){

    Grid test_grid(4, 3, 2.0, 1.0);
    Field test_field(test_grid, 1.5);

    SECTION("Field size and indexing") {
        CHECK(test_field.size() == 30);
        CHECK(test_field.index(0, 0) == 0);
        CHECK(test_field.index(test_grid.getNi() - 1, test_grid.getNj() - 1) 
                                == static_cast<std::size_t>(test_grid.size() - 1));
        CHECK(test_field.index(1, 2) - test_field.index(0, 2) == 1);
        CHECK(test_field.index(0, 1) - test_field.index(0, 0) == static_cast<std::size_t>(test_grid.getNi()));
    }

    SECTION("Field access") {
        CHECK(test_field(0, 0) == 1.5);
        CHECK(test_field(3, 2) == 1.5);
        CHECK(test_field(test_grid.getNi() - 1, test_grid.getNj() - 1) == 1.5);
    }

    SECTION("Ghost cells and modification") {

        test_field(2, 1) = 7.0;
        CHECK(test_field(2, 1) == 7.0);
        const Field& cf = test_field;
        CHECK(cf(2, 1) == 7.0);

        test_field(0, 2) = 9.0;
        CHECK(test_field(1, 2) == 1.5);
        CHECK(test_field(test_grid.getNi() - 1, 1) == 1.5);
        CHECK(test_field.values()[test_field.index(2, 1)] == 7.0);
        test_field.fill(-2.0);
        CHECK(test_field(0, 0) == -2.0);
        CHECK(test_field(1, 2) == -2.0);
        CHECK(test_field(test_grid.getNi() - 1, test_grid.getNj() - 1) == -2.0);

        CHECK(test_field.index(0, 2) - test_field.index(test_grid.getNi() - 1, 1) == 1);
    }

    SECTION("Field independence"){
        Field another_field(test_grid, 0.0);
        CHECK(another_field(0, 0) == 0.0);
        CHECK(another_field(3, 2) == 0.0);
        CHECK(another_field(test_grid.getNi() - 1, test_grid.getNj() - 1) == 0.0);

        test_field(1, 1) = 5.0;
        CHECK(another_field(1, 1) == 0.0); // Ensure another_field remains unchanged
        CHECK(test_field(1, 1) == 5.0); // Ensure test_field is updated correctly
    }

}