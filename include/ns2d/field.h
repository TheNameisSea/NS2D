#pragma once
#include <vector>
#include "ns2d/grid.h"
#include <cassert>

class Field {

    private:
        int Ni{};
        int Nj{};
        std::vector<double> data;

    public:
        // Constructor
        Field(const Grid& grid, double initial=0.0);
        std::size_t size() const;
        std::size_t index(int i, int j) const{
            assert(i >= 0 && i < Ni);
            assert(j >= 0 && j < Nj);
    
            return static_cast<std::size_t>(i + Ni * j);
        }

        // Access operator
        double& operator()(int i, int j) {
            return data[index(i, j)];
        }

        // Const access operator
        double operator()(int i, int j) const {
            return data[index(i, j)];
        }

        // Fill the field with a specific value
        void fill(double value);

        // Getters for Ni and Nj
        int getNi() const;
        int getNj() const;

        // Access to the underlying data
        const std::vector<double>& values() const{
            return data;
        }

        // Non-const access to the underlying data
        std::vector<double>& values(){
            return data;
        }

};