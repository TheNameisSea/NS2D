#pragma once
#include "ns2d/field.h"
#include "ns2d/grid.h"

enum class Staggering {

    Cell,
    UFace,
    VFace

};

struct IndexRange {

    int iBegin{};
    int iEnd{};
    int jBegin{};
    int jEnd{};

};

IndexRange interiorRange(const Grid&, Staggering);

void laplacian(const Field& f, const Grid& g, Staggering s, Field& out);


