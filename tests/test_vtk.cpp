#include "ns2d/grid.h"
#include "ns2d/field.h"
#include "ns2d/io/vtk_writer.h"
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> readLines(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

// Index of the first line equal to `text`, or -1
int findLine(const std::vector<std::string>& lines, const std::string& text) {
    for (std::size_t k = 0; k < lines.size(); ++k) {
        if (lines[k] == text) {
            return static_cast<int>(k);
        }
    }
    return -1;
}

} // namespace

TEST_CASE("writeVtk writes a valid legacy VTK file", "[vtk]") {
    const Grid g(3, 2, 3.0, 1.0);   // dx = 1, dy = 0.5
    Field u(g, 0.0);
    Field v(g, 0.0);
    Field p(g, 0.0);
    for (int j = 1; j <= g.getNy(); ++j) {
        for (int i = 1; i <= g.getNx(); ++i) {
            p(i, j) = 10.0 * j + i;   // easy to check the ordering
        }
    }

    const std::string path = (std::filesystem::temp_directory_path() / "ns2d_test.vtk").string();
    writeVtk(path, g, u, v, p, 1.5);
    const auto lines = readLines(path);
    REQUIRE(lines.size() > 10);

    SECTION("header") {
        CHECK(lines[0] == "# vtk DataFile Version 3.0");
        CHECK(lines[1] == "ns2d t=1.5");
        CHECK(lines[2] == "ASCII");
        CHECK(lines[3] == "DATASET STRUCTURED_POINTS");
        CHECK(lines[4] == "DIMENSIONS 4 3 1");      // points = cells + 1
        CHECK(lines[5] == "ORIGIN 0 0 0");
        CHECK(lines[6] == "SPACING 1 0.5 1");
        CHECK(lines[7] == "CELL_DATA 6");
    }

    SECTION("pressure values in x-fastest order") {
        const int k = findLine(lines, "SCALARS pressure double 1");
        REQUIRE(k >= 0);
        CHECK(lines[static_cast<std::size_t>(k + 1)] == "LOOKUP_TABLE default");
        const std::vector<std::string> expected{"11", "12", "13", "21", "22", "23"};
        for (std::size_t n = 0; n < expected.size(); ++n) {
            INFO("value " << n);
            CHECK(lines[static_cast<std::size_t>(k + 2) + n] == expected[n]);
        }
    }

    SECTION("all three data sets are present with one entry per cell") {
        const int kw = findLine(lines, "SCALARS vorticity double 1");
        const int kv = findLine(lines, "VECTORS velocity double");
        REQUIRE(kw >= 0);
        REQUIRE(kv >= 0);
        CHECK(kv - kw == 2 + 6);                                   // header lines + 6 values
        CHECK(static_cast<int>(lines.size()) - kv - 1 == 6);       // 6 vectors at the end
        CHECK(lines.back() == "0 0 0");
    }

    std::filesystem::remove(path);
}

TEST_CASE("writeVtk throws when the file cannot be created", "[vtk]") {
    const Grid g(2, 2, 1.0, 1.0);
    const Field f(g, 0.0);
    CHECK_THROWS_AS(writeVtk("/nonexistent/dir/x.vtk", g, f, f, f, 0.0), std::runtime_error);
}
