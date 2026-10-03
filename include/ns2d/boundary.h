#pragma once
#include "ns2d/field.h"
#include "ns2d/grid.h"

struct BoundaryConditions {

    // Top is the tanential wall velocity
    // u along the top and bottom walls
    double top{0.0};
    double bottom{0.0};
    // v along the left and right walls
    double left{0.0};
    double right{0.0};

    // Apply no slip wall
    void applyNoSlipWall(Field& u, Field& v, const Grid& g) const;

};