#include "ns2d/grid.h"
#include <stdexcept>


Grid::Grid(int nx_, int ny_, double Lx_, double Ly_) : nx(nx_), ny(ny_), Lx(Lx_), Ly(Ly_), dx(Lx_ / nx_), dy(Ly_ / ny_), Ni(nx_ + 2 * ghost), Nj(ny_ + 2 * ghost) {

    if (!(nx > 0 && ny > 0 && Lx > 0.0 && Ly > 0.0)) {
        throw std::invalid_argument("Invalid grid parameters");
    }
}

int Grid::getNx() const { return nx; }
int Grid::getNy() const { return ny; }
double Grid::getLx() const { return Lx; }
double Grid::getLy() const { return Ly; }
double Grid::getDx() const { return dx; }
double Grid::getDy() const { return dy; }
int Grid::getNi() const { return Ni; }
int Grid::getNj() const { return Nj; }

// Pressure position
double Grid::xp(int i) const {
    return (i - 0.5) * dx;
}

double Grid::yp(int j) const {
    return (j - 0.5) * dy;
}


// u velocity position
double Grid::xu(int i) const {
    return i * dx;
}


// v velocity position
double Grid::yv(int j) const {
    return j * dy;
}

int Grid::size() const {
    return Ni * Nj;
}



