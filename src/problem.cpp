#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "errhandle.h"
#include "problem.h"
#include "settings.h"

using namespace std;

Problem::Problem()
{
	initCondFile = "/tmp/heat/u0.bin";
	outputFile = "/tmp/heat/u.bin";
	numberOfIterationsToSkip = SETTINGS_NUM_OF_ITERATIONS_TO_SKIP;
	maxt = SETTINGS_MAXT;
	tau = SETTINGS_TAU;
	maxstep = maxt/tau;
	maxm = maxstep/numberOfIterationsToSkip;
	print();
}

void Problem::print()
{
 	cout << "==================== Files: ====================" << endl;
 	cout << "init cond file: " << initCondFile << endl;
 	cout << "output file   : " << outputFile << endl;
 	cout << "================================================" << endl;
 	cout << "================= Parameters: ==================" << endl;
 	cout << "tau  = " << tau << endl;
 	cout << "maxt = " << maxt << endl;
 	cout << "iters2skip = " << numberOfIterationsToSkip << endl;
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

void Problem1d::print()
{
	cout << "Nx   = " << Nx << endl;
	cout << "Dx   = " << Dx << endl;
	cout << "hx    = " << 1.0 / Nx << endl;
	cout << "boundary cond type at x0 = " << boundCondTypeX0 << endl;
	cout << "boundary cond type at x1 = " << boundCondTypeX1 << endl;
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

void Problem2d::print()
{
	cout << "Ny   = " << Ny << endl;
	cout << "Dy   = " << Dy << endl;
	cout << "hy    = " << 1.0 / Ny << endl;
	cout << "boundary cond type at y0 = " << boundCondTypeY0 << endl;
	cout << "boundary cond type at y1 = " << boundCondTypeY1 << endl;
}

Problem3d::Problem3d()
{
	boundCondTypeZ0 = SETTINGS_BCZ0_TYPE;
	boundCondTypeZ1 = SETTINGS_BCZ1_TYPE;
	Dz = SETTINGS_DZ;
	Nz = SETTINGS_NZ;
	totalSize = (Nx+1)*(Ny+1)*(Nz+1);
	print();
}

void Problem3d::print()
{
	cout << "Nz   = " << Nz << endl;
	cout << "Dz   = " << Dz << endl;
	cout << "hz    = " << 1.0 / Nz << endl;
	cout << "boundary cond type at z0 = " << boundCondTypeZ0 << endl;
	cout << "boundary cond type at z1 = " << boundCondTypeZ1 << endl;
	cout << endl;
}

///////////////////////////////////////////////////////////

double Problem1d::qx0(const int n, const int j, const int k)
{
	return SETTINGS_BCQX0_FUNCTION(n,j,k);
}

double Problem1d::qx1(const int n, const int j, const int k)
{
	return SETTINGS_BCQX1_FUNCTION(n,j,k);
}

double Problem1d::px0(const int n, const int j, const int k)
{
	return SETTINGS_BCPX0_FUNCTION(n,j,k);
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

double Problem2d::qy1(const int n, const int i, const int k)
{
	return SETTINGS_BCQY1_FUNCTION(n,j,k);
}

double Problem2d::py0(const int n, const int i, const int k)
{
	return SETTINGS_BCPY0_FUNCTION(n,j,k);
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

double Problem3d::qz1(const int n, const int i, const int j)
{
	return SETTINGS_BCQZ1_FUNCTION(n,j,k);
}

double Problem3d::pz0(const int n, const int i, const int j)
{
	return SETTINGS_BCPZ0_FUNCTION(n,j,k);
}

double Problem3d::pz1(const int n, const int i, const int j)
{
	return SETTINGS_BCPZ1_FUNCTION(n,j,k);
}
