// Author: Evgeny Berendeev
// Budker Institute of Nuclear Physics of Siberian Branch Russian Academy of Sciences
// beren@inp.nsk.su
// (c) 2022, for licensing details see the LICENSE file

#pragma once

#ifndef UTIL_H
#define UTIL_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <algorithm>

// Define basic types
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> Operator;
typedef Eigen::VectorXd Field;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;

const double dt = 0.1;
const double dx = 0.2;
const double dy = 0.2;
const double dz = 0.2;

const int cellsX = 2;
const int cellsY = 2;
const int cellsZ = 2;
const Vector3d totalSize(cellsX+3,cellsY+3,cellsZ+3);
//general indexing routine (row major)
const int Nc = 3;
const int Nx = cellsX+3;
const int Ny = cellsY+3;
const int Nz = cellsZ+3;
inline constexpr int vind(int x, int y, int z, int c)
{
	return(z + Nz*(y + Ny*(x + Nx*c)));
}
void stencil_identity(Operator& identity);
void stencil_mat_xx(vector<Trip>& trips);
#endif
