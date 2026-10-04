#include "ns2d/poisson.h"

void applyPressureNeumann(Field& p, const Grid& grid){

    const int nx = grid.getNx();
    const int ny = grid.getNy();

    for (int j = 1; j <= ny; ++j){
        p(0, j) = p(1, j);
        p(nx+1, j) = p(nx, j);
            
    }

    for (int i = 1; i <= nx; ++i){
        p(i, 0) = p(i, 1);
        p(i, ny+1) = p(i, ny);
    }

}

