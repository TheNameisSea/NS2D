#include "ns2d/io/config.h"
#include "ns2d/operators.h"
#include "ns2d/simulation.h"
#include <fstream> 
#include <stdexcept>

namespace {
    AdvectionScheme schemeFromString(const std::string& s){
        if (s == "central"){
            return AdvectionScheme::Central;
        }
        if (s == "upwind2"){
            return AdvectionScheme::Upwind2;
        }
        throw std::invalid_argument("unknown advection scheme: " + s);
    }
}

CaseConfig parseConfig(const nlohmann::json& j){

    // default
    double steadyTol = 1e-6;   
    int maxSteps = 100000;   
    int printEvery = 100;   

    const auto& jg = j.at("grid");                         // required section
    Grid grid( jg.at("nx").get<int>(), jg.at("ny").get<int>(),
            jg.at("Lx").get<double>(), jg.at("Ly").get<double>() );

    SimulationParams prm{};                                // start from the defaults
    prm.Re = j.at("Re").get<double>();                 // required
    prm.cfl = j.value("cfl", prm.cfl);                  // optional: keep default if missing
    prm.scheme = schemeFromString(j.value("scheme", std::string{"upwind2"}));
    prm.bc = BoundaryConditions{ .top = j.value("lidVelocity", 1.0) };

    // j.contains("poisson")
    if (j.contains("poisson")){                                   // j.contains("poisson")
        const auto& jp = j.at("poisson");
        prm.poissonTol = jp.value("tol", prm.poissonTol);
        prm.poissonMaxIter = jp.value("maxIter", prm.poissonMaxIter);
        if (jp.value("solver", std::string{"jacobi"}) != "jacobi"){
            throw std::invalid_argument("Invalid solver type");
        }
        

    }
    
    if (j.contains("steady")){                                   // j.contains("poisson")
        const auto& js = j. at("steady");
        steadyTol = js.value("tol", steadyTol);
        maxSteps = js.value("maxSteps", maxSteps);

    }

    if (j.contains("output")){
        const auto& jo = j.at("output");
        printEvery = jo.value("printEvery", printEvery);
    }

    return CaseConfig{ grid, prm, steadyTol, maxSteps, printEvery };

}

CaseConfig loadConfig(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open config file: " + path);
    }
    return parseConfig(nlohmann::json::parse(file));

}