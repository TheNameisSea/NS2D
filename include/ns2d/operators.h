#pragma once
#include "ns2d/field.h"
#include "ns2d/grid.h"

enum class Staggering {

    Cell,
    UFace,
    VFace

};

enum class AdvectionScheme {
    Central,
    Upwind2
};

struct IndexRange {

    int iBegin{};
    int iEnd{};
    int jBegin{};
    int jEnd{};

};

IndexRange interiorRange(const Grid&, Staggering);

void laplacian(const Field& f, const Grid& g, Staggering s, Field& out);

inline double uAtV(const Field& u, int i, int j){
    return  0.25 * ( u(i-1, j+1) + u(i, j+1) + u(i, j) + u(i-1, j) );
}
inline double vAtU(const Field& v, int i, int j){
    return  0.25 * ( v(i, j) + v(i+1, j) + v(i, j-1) + v(i+1, j-1) );
}

double derivative(const Field& f, int i, int j, int di, int dj, double h, double a, AdvectionScheme scheme);

// out holds +(u·∇)u
void advectionU(const Field& u, const Field& v, const Grid& g, AdvectionScheme scheme, Field& out);
void advectionV(const Field& u, const Field& v, const Grid& g, AdvectionScheme scheme, Field& out);



