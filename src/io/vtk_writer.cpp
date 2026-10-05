#include "ns2d/io/vtk_writer.h"
#include "ns2d/operators.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

void writeVtk(const std::string& path, const Grid& g,
              const Field& u, const Field& v, const Field& p, double time) {

    // Move everything to the cell centers
    Field uc(g);
    Field vc(g);
    Field omega(g);
    cellCenteredVelocity(u, v, g, uc, vc);
    vorticity(u, v, g, omega);

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write VTK file: " + path);
    }

    const int nx = g.getNx();
    const int ny = g.getNy();

    out << std::setprecision(10);
    out << "# vtk DataFile Version 3.0\n";
    out << "ns2d t=" << time << '\n';
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << nx + 1 << ' ' << ny + 1 << " 1\n";   // points = cells + 1
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << g.getDx() << ' ' << g.getDy() << " 1\n";
    out << "CELL_DATA " << nx * ny << '\n';

    // VTK wants x varying fastest: j outer, i inner over the interior cells
    auto writeScalar = [&](const char* name, const Field& f) {
        out << "SCALARS " << name << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                out << f(i, j) << '\n';
            }
        }
    };

    writeScalar("pressure", p);
    writeScalar("vorticity", omega);

    out << "VECTORS velocity double\n";
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            out << uc(i, j) << ' ' << vc(i, j) << " 0\n";
        }
    }

    if (!out) {
        throw std::runtime_error("error while writing VTK file: " + path);
    }
}
