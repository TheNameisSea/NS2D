#pragma once
#include "ns2d/field.h"
#include "ns2d/grid.h"
#include "ns2d/operators.h"

struct PredictorWork {
    Field advU, advV, lapU, lapV;

    explicit PredictorWork(const Grid& g);
};

void predictor(const Field& u, const Field& v, const Grid& g,
               double dt, double Re, AdvectionScheme scheme,
               PredictorWork& work, Field& uStar, Field& vStar);


