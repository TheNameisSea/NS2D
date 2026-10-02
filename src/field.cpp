#include "ns2d/field.h"
#include <algorithm>

Field::Field(const Grid& grid, double initial) 
            : Ni(grid.getNi()), Nj(grid.getNj()), 
            data(static_cast<std::size_t>(grid.size()), initial) {

}

std::size_t Field::size() const {
    return data.size();
}

void Field::fill(double value) {
    std::fill(data.begin(), data.end(), value);
}

int Field::getNi() const {
    return Ni;
}

int Field::getNj() const {
    return Nj;
}


