#pragma once

// Use MAC grid
// Pressure at cell centers, 
// u-velocity at vertical cell faces, 
// v-velocity at horizontal cell faces
// u of cell is at the right face of the cell, 
// v of cell is at the top face of the cell

class Grid {

    private:
        int nx{};
        int ny{};
        double Lx{};
        double Ly{};

        double dx{};
        double dy{};

        // Including ghost cells
        int Ni{};
        int Nj{};

        static constexpr int ghost{1};

    public:
        // Constructor
        Grid(int nx, int ny, double Lx, double Ly);

        // Getters
        int getNx() const;
        int getNy() const;
        double getLx() const;
        double getLy() const;
        double getDx() const;
        double getDy() const;
        int getNi() const;
        int getNj() const;

        // Pressure position
        double xp(int i) const;
        double yp(int j) const;

        // u velocity position
        double xu(int i) const;

        // v velocity position
        double yv(int j) const;

        // Size of the grid (including ghost cells)
        int size() const;

};
