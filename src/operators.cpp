#include "ns2d/operators.h"
#include "ns2d/grid.h"
#include <cassert>

IndexRange interiorRange(const Grid& g, Staggering s){

    const int nx = g.getNx();
    const int ny = g.getNy();

    switch (s)
    {
    case Staggering::Cell:
        return IndexRange{1, nx, 1, ny};
    
    case Staggering::UFace:
        return IndexRange{1, nx-1, 1, ny};

    case Staggering::VFace:
        return IndexRange{1, nx, 1, ny-1};
    }

    assert(false && "unknown staggering");

    return IndexRange{};

}

void laplacian(const Field& f, const Grid& g, Staggering s, Field& out){

    assert((out.getNi() == f.getNi()) && (out.getNj() == f.getNj()));

    const auto [i0, i1, j0, j1] = interiorRange(g, s);
    
    const double dx = g.getDx();
    const double dy = g.getDy();

    // Precompute coefficients
    const double idx2 = 1.0 / (dx * dx);
    const double idy2 = 1.0 / (dy * dy);

    for (int j = j0; j <= j1; ++j){

        for (int i = i0; i <= i1; ++i){

            const double c = f(i, j);
            out(i, j) = (f(i+1, j) - 2.0*c + f(i-1, j)) * idx2
                        + (f(i, j+1) - 2.0*c + f(i, j-1)) * idy2;

        }

    }
}


