#include "ns2d/boundary.h"

void BoundaryConditions::applyNoSlipWall(Field& u, Field& v, const Grid& g) const {

    const int nx = g.getNx();
    const int ny = g.getNy();

    // Phase 1: Normal velocities
    for (int j = 0; j <= ny+1; ++j) {
        u(0, j) = 0; // Left wall
        u(nx, j) = 0; // Right wall
    }
    for (int i = 0; i <= nx+1; ++i) {
        v(i, 0) = 0; // Bottom wall
        v(i, ny) = 0; // Top wall
    }

    // Phase 2: Tangential velocities via ghosts
    for (int i=0; i <= nx; ++i) {
        u(i, 0) = 2 * bottom - u(i, 1); // Bottom wall
        u(i, ny+1) = 2*top - u(i, ny); // Top wall
    }
    for (int j=0; j <= ny; ++j) {
        v(0, j) = 2 * left - v(1, j); // Left wall
        v(nx+1, j) = 2 * right - v(nx, j); // Right wall
    }


}