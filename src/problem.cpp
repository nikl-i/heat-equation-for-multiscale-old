#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "errhandle.h"
#include "problem.h"
#include "settings.h"

#ifndef SETTINGS_NUM_OF_ITERATIONS_TO_SKIP
#define SETTINGS_NUM_OF_ITERATIONS_TO_SKIP 1
#endif

#ifndef SETTINGS_BCX0_TYPE
#define SETTINGS_BCX0_TYPE 0
#endif

#ifndef SETTINGS_BCQX0_FUNCTION
#define SETTINGS_BCQX0_FUNCTION(n,j,k) 0.0
#endif

#ifndef SETTINGS_BCPX0_FUNCTION
#define SETTINGS_BCPX0_FUNCTION(n,j,k) 0.0
#endif

#ifndef SETTINGS_BCX1_TYPE
#define SETTINGS_BCX1_TYPE 0
#endif

#ifndef SETTINGS_BCQX1_FUNCTION
#define SETTINGS_BCQX1_FUNCTION(n,j,k) 0.0
#endif

#ifndef SETTINGS_BCPX1_FUNCTION
#define SETTINGS_BCPX1_FUNCTION(n,j,k) 0.0
#endif

#ifndef SETTINGS_HEATSRC1D_FUNCTION
#define SETTINGS_HEATSRC1D_FUNCTION(U,n,i) 0.0
#endif
////////////////////////////////////////////////////////////////////////////////
#ifndef SETTINGS_DY
#define SETTINGS_DY 0.0
#endif

#ifndef SETTINGS_NY
#define SETTINGS_NY 0
#endif

#ifndef SETTINGS_BCY0_TYPE
#define SETTINGS_BCY0_TYPE 0
#endif

#ifndef SETTINGS_BCQY0_FUNCTION
#define SETTINGS_BCQY0_FUNCTION(n,i,k) 0.0
#endif

#ifndef SETTINGS_BCPY0_FUNCTION
#define SETTINGS_BCPY0_FUNCTION(n,i,k) 0.0
#endif

#ifndef SETTINGS_BCY1_TYPE
#define SETTINGS_BCY1_TYPE 0
#endif

#ifndef SETTINGS_BCQY1_FUNCTION
#define SETTINGS_BCQY1_FUNCTION(n,i,k) 0.0
#endif

#ifndef SETTINGS_BCPY1_FUNCTION
#define SETTINGS_BCPY1_FUNCTION(n,i,k) 0.0
#endif

#ifndef SETTINGS_HEATSRC2D_FUNCTION
#define SETTINGS_HEATSRC2D_FUNCTION(U,n,i,j) 0.0
#endif

////////////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_DZ
#define SETTINGS_DZ 0.0
#endif

#ifndef SETTINGS_NZ
#define SETTINGS_NZ 0
#endif

#ifndef SETTINGS_BCZ0_TYPE
#define SETTINGS_BCZ0_TYPE 0
#endif

#ifndef SETTINGS_BCQZ0_FUNCTION
#define SETTINGS_BCQZ0_FUNCTION(n,i,j) 0.0
#endif

#ifndef SETTINGS_BCPZ0_FUNCTION
#define SETTINGS_BCPZ0_FUNCTION(n,i,j) 0.0
#endif

#ifndef SETTINGS_BCZ1_TYPE
#define SETTINGS_BCZ1_TYPE 0
#endif

#ifndef SETTINGS_BCQZ1_FUNCTION
#define SETTINGS_BCQZ1_FUNCTION(n,i,j) 0.0
#endif

#ifndef SETTINGS_BCPZ1_FUNCTION
#define SETTINGS_BCPZ1_FUNCTION(n,i,j) 0.0
#endif

#ifndef SETTINGS_HEATSRC_FUNCTION
#define SETTINGS_HEATSRC_FUNCTION(U,n,i,j,k) 0.0
#endif


#ifndef SETTINGS_ELECTRIC_FILE
#define SETTINGS_ELECTRIC_FILE ""
#endif

////////////////////////////////////////////////////////////////////////////////

using namespace std;

Problem::Problem()
{
	initCondFile = SETTINGS_INITCOND_FILE;
	outputFile = SETTINGS_SOLUTION_FILE;
	electricIntensityFile = SETTINGS_ELECTRIC_FILE;

	numberOfIterationsToSkip = SETTINGS_NUM_OF_ITERATIONS_TO_SKIP;
	maxt = SETTINGS_MAXT;
	tau = SETTINGS_TAU;
	maxstep = maxt/tau;
	maxm = maxstep/numberOfIterationsToSkip;
	print();
}

void Problem::loadElectricIntensity()
{
	if (electricIntensityFile.empty()) return;

	E = new double [totalSize];
	cout << "Loading electric field condition... ";
	ifstream input(electricIntensityFile.c_str(), ios::binary|ios::in);
	input.read((char *) E, totalSize * sizeof(double));
	cout << "Done." << endl;
	return;
}

void Problem::clearElectricIntensity()
{
	if (E != NULL) delete [] E;
	return;
}

void Problem::print() const
{
 	cout << "================================ Files: ========================================" << endl;
 	cout << "Initial condidion file: " << initCondFile << endl;
 	cout << "Output file   : " << outputFile << endl;
 	cout << "=============================== Parameters: ====================================" << endl;
 	cout << "tau  = " << tau << endl;
 	cout << "maxt = " << maxt << endl;
 	cout << "Number of iterations to skip = " << numberOfIterationsToSkip << endl;
}

Problem1d::Problem1d()
{
	boundCondTypeX0 = SETTINGS_BCX0_TYPE;
	boundCondTypeX1 = SETTINGS_BCX1_TYPE;
	Dx = SETTINGS_DX;
	Nx = SETTINGS_NX;
	totalSize = Nx + 1;
	print();
}

void Problem1d::print() const
{
	cout << "Nx   = " << Nx << endl;
	cout << "Dx   = " << Dx << endl;
	cout << "hx    = " << 1.0 / Nx << endl;
	cout << "Boundary condition type at x0 = " << boundCondTypeX0 << endl;
	cout << "Boundary condition type at x1 = " << boundCondTypeX1 << endl;
}

Problem2d::Problem2d()
{
	boundCondTypeY0 = SETTINGS_BCY0_TYPE;
	boundCondTypeY1 = SETTINGS_BCY1_TYPE;
	Dy = SETTINGS_DY;
	Ny = SETTINGS_NY;
	totalSize = (Nx+1)*(Ny+1);
	print();

}

void Problem2d::print() const
{
	cout << "Ny   = " << Ny << endl;
	cout << "Dy   = " << Dy << endl;
	cout << "hy    = " << 1.0 / Ny << endl;
	cout << "Boundary condition type at y0 = " << boundCondTypeY0 << endl;
	cout << "Boundary condition type at y1 = " << boundCondTypeY1 << endl;
}

Problem3d::Problem3d()
{
	boundCondTypeZ0 = SETTINGS_BCZ0_TYPE;
	boundCondTypeZ1 = SETTINGS_BCZ1_TYPE;
	Dz = SETTINGS_DZ;
	Nz = SETTINGS_NZ;
	totalSize = (Nx+1)*(Ny+1)*(Nz+1);
	print();

#ifdef SETTINGS_ELECTRIC_FILE
#endif
}

void Problem3d::print() const
{
	cout << "Nz   = " << Nz << endl;
	cout << "Dz   = " << Dz << endl;
	cout << "hz    = " << 1.0 / Nz << endl;
	cout << "Boundary condition type at z0 = " << boundCondTypeZ0 << endl;
	cout << "Boundary condition type at z1 = " << boundCondTypeZ1 << endl;
	cout << endl;
}

///////////////////////////////////////////////////////////
double Problem1d::heatSrc(const double *U, const int n, const int i)
{
	return SETTINGS_HEATSRC1D_FUNCTION(U,n,i);
}

double Problem2d::heatSrc(const double *U, const int n, const int i,
						const int j)
{
	return SETTINGS_HEATSRC2D_FUNCTION(U,n,i,j);
}

double Problem3d::heatSrc(const double *U, const int n, const int i,
						const int j, const int k)
{
	return SETTINGS_HEATSRC_FUNCTION(U,n,i,j,k);
}

///////////////////////////////////////////////////////////

double Problem1d::qx0(const int n, const int j, const int k)
{
	return SETTINGS_BCQX0_FUNCTION(n,j,k);
}

double Problem1d::px0(const int n, const int j, const int k)
{
	return SETTINGS_BCPX0_FUNCTION(n,j,k);
}

double Problem1d::qx1(const int n, const int j, const int k)
{
	return SETTINGS_BCQX1_FUNCTION(n,j,k);
}

double Problem1d::px1(const int n, const int j, const int k)
{
	return SETTINGS_BCPX1_FUNCTION(n,j,k);
}

///////////////////////////////////////////////////////////

double Problem2d::qy0(const int n, const int i, const int k)
{
	return SETTINGS_BCQY0_FUNCTION(n,j,k);
}

double Problem2d::py0(const int n, const int i, const int k)
{
	return SETTINGS_BCPY0_FUNCTION(n,j,k);
}

double Problem2d::qy1(const int n, const int i, const int k)
{
	return SETTINGS_BCQY1_FUNCTION(n,j,k);
}

double Problem2d::py1(const int n, const int i, const int k)
{
	return SETTINGS_BCPY1_FUNCTION(n,j,k);
}

///////////////////////////////////////////////////////////

double Problem3d::qz0(const int n, const int i, const int j)
{
	return SETTINGS_BCQZ0_FUNCTION(n,j,k);
}

double Problem3d::pz0(const int n, const int i, const int j)
{
	return SETTINGS_BCPZ0_FUNCTION(n,j,k);
}

double Problem3d::qz1(const int n, const int i, const int j)
{
	return SETTINGS_BCQZ1_FUNCTION(n,j,k);
}

double Problem3d::pz1(const int n, const int i, const int j)
{
	return SETTINGS_BCPZ1_FUNCTION(n,j,k);
}
