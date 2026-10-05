#pragma once
#include "ns2d/field.h"
#include "ns2d/grid.h"
#include <string>

// Write pressure, vorticity and velocity at the cell centers to a legacy ASCII VTK file
// (STRUCTURED_POINTS, CELL_DATA) that ParaView can open. u and v must have their BCs applied.
// Throws std::runtime_error if the file cannot be written.
void writeVtk(const std::string& path, const Grid& g,
              const Field& u, const Field& v, const Field& p, double time);
