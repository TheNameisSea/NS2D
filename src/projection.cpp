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
