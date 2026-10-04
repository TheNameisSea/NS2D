#include "ns2d/projection.h"
#include "ns2d/field.h"
#include "ns2d/operators.h"
#include <cassert>

    
PredictorWork::PredictorWork(const Grid& g) : advU(g), advV(g), lapU(g), lapV(g){
}

void predictor(const Field& u, const Field& v, const Grid& g,
               double dt, double Re, AdvectionScheme scheme,
               PredictorWork& work, Field& uStar, Field& vStar) {
    
    assert((u.size() == v.size()) && (u.size() == uStar.size()) && (v.size() == vStar.size()));
    assert(&uStar != &u);
    assert(&vStar != &v);

    uStar = u;
    vStar = v;

    advectionU(u, v, g, scheme, work.advU);
    advectionV(u, v, g, scheme, work.advV);
    laplacian(u, g, Staggering::UFace, work.lapU);
    laplacian(v, g, Staggering::VFace, work.lapV);

    const double nu = 1.0 / Re;  
    {                             
        auto [i0, i1, j0, j1] = interiorRange(g, Staggering::UFace);
        for (int j = j0; j <= j1; ++j){
            for (int i = i0;  i<= i1; ++i){
                uStar(i, j) = u(i, j) + dt * ( -work.advU(i, j) + nu * work.lapU(i, j));
            }
        }
    }

    {
        auto [i0, i1, j0, j1] = interiorRange(g, Staggering::VFace);
        for (int j = j0; j <= j1; ++j){
            for (int i = i0;  i<= i1; ++i){
            vStar(i, j) = v(i, j) + dt * ( -work.advV(i, j) + nu * work.lapV(i, j));
            }
        }
    }

                
}

void divergence(const Field& u, const Field& v, const Grid& g, Field& out) {

    const double dx = g.getDx();
    const double dy = g.getDy();

    const auto[i0, i1, j0, j1] = interiorRange(g, Staggering::Cell);
    for (int j = j0; j <= j1; ++j){
        for (int i = i0; i <= i1; ++i){
            out(i,j) = ( u(i,j) - u(i-1,j) ) / dx + ( v(i,j) - v(i,j-1) ) / dy;
        }
    }
}

void corrector(const Field& uStar, const Field& vStar, const Field& p,
               const Grid& g, double dt, Field& u, Field& v){
    
    assert((u.size() == v.size()) && (u.size() == uStar.size()) && (v.size() == vStar.size()));
    assert(&uStar != &u);
    assert(&vStar != &v);
    
    u = uStar;
    v = vStar;
    
    const double dt_dx = dt / g.getDx();
    const double dt_dy = dt / g.getDy();

    {
        const auto[i0, i1, j0, j1] = interiorRange(g, Staggering::UFace);
        for (int j = j0; j <= j1; ++j){
            for (int i = i0; i <= i1; ++i){
                u(i,j) = uStar(i,j) - dt_dx * ( p(i+1,j) - p(i,j) );
            }
        }
    }

    {
        const auto[i0, i1, j0, j1] = interiorRange(g, Staggering::VFace);
        for (int j = j0; j <= j1; ++j){
            for (int i = i0; i <= i1; ++i){
                v(i,j) = vStar(i,j) - dt_dy * ( p(i,j+1) - p(i,j) );
            }
        }
    }

}
