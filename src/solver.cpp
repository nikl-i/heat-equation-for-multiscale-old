#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "errhandle.h"
#include "lapacke.h"
#include "problem.h"
#include "solution.h"
#include "solver.h"
using namespace std;

Solver::Solver(Problem1d problem)
{
	ipiv1 = new lapack_int [problem.Nx+1];
	Un = new double [problem.totalSize];  // current time step
	B = new double [problem.totalSize];
	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);

	cout << "Initializiation of the solver... ";
	hx = 1.0 / problem.Nx;
	double coef = 0.5 * problem.tau / (hx*hx);
	coef_ax = coef * problem.Dx;
	coef_cx = coef * problem.Dx;
	coef_bx = coef_ax + coef_cx;

	nrhs = 1;
	ldb = nrhs;
	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
			   get_alphaX(problem), get_betaX(problem), get_gammaX(problem), get_deltaX(problem));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	cout << "Done."  << endl;
	cout << "================================================================================" << endl;

	return;
}

void Solver::solve(Problem1d problem, Solution *U)
{
	cout << "Solving the problem... "  << endl;
	memcpy(Un, U->get_data(0), problem.totalSize * sizeof(double));

	for (int n = 0; n < problem.maxstep; n++)
	{
		if ((problem.boundCondTypeX0 == 3) || (problem.boundCondTypeX1 == 3))
		{
			cout << "Matrix update... ";
			initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
					   get_alphaX(problem, n), get_betaX(problem), get_gammaX(problem), get_deltaX(problem, n));
			check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
			cout << "Done." << endl;
		}

		for (int i = 1; i < problem.Nx; i++)
			B[i] = coef_ax*Un[i-1] + coef_cx*Un[i+1] + (1.0 - coef_bx)*Un[i] + problem.heatSrc(Un,n,i);
		B[0] = get_b0X(problem, Un, n);
		B[problem.Nx] = get_bNX(problem, Un, n);

		check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
								   nrhs, dl1, d1, du1, du21, ipiv1, B, ldb));
		switchPointers(&Un, &B);

		if ((n+1) % problem.numberOfIterationsToSkip == 0)
		{
			cout << "Iteration n = " << n << " done." << endl;
			memcpy(U->get_data((n+1) / problem.numberOfIterationsToSkip), Un, problem.totalSize * sizeof(double));
		}
	}
	cout << "Done."  << endl;
	cout << "================================================================================" << endl;
	return;
}


Solver::Solver(Problem2d problem)
{
	cout << "Initializiation of the 2d cpu solver... ";
	hx = 1.0 / problem.Nx;
	double coef = 0.5 * problem.tau / (hx*hx);
	coef_ax = coef * problem.Dx;
	coef_cx = coef * problem.Dx;
	coef_bx = coef_ax + coef_cx;

	hy = 1.0 / problem.Ny;
	coef = 0.5 * problem.tau / (hy*hy);
	coef_ay = coef * problem.Dy;
	coef_cy = coef * problem.Dy;
	coef_by = coef_ay + coef_cy;

	nrhs = 1;
	ldb = nrhs;
	ipiv1 = new lapack_int [problem.Nx+1];
	ipiv2 = new lapack_int [problem.Ny+1];

	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);
	allocateMatrix(&dl2, &d2, &du2, &du22, problem.Ny+1);

	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
			   get_alphaX(problem), get_betaX(problem), get_gammaX(problem), get_deltaX(problem));
	initMatrix(dl2, d2, du2, problem.Ny+1, -coef_ay, 1.0 + coef_by, -coef_cy,
			   get_alphaY(problem), get_betaY(problem), get_gammaY(problem), get_deltaY(problem));

	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));
	Un = new double [problem.totalSize];  // current time step. Will be initialized later
	B = new double [problem.totalSize];
	memset(B, 0, problem.totalSize*sizeof(double));
	temp = new double [problem.totalSize];

	cout << "Done."  << endl;
	cout << "================================================================================" << endl;

	return;
}

void Solver::solve(Problem2d problem, Solution *U)
{
	cout << "Solving the problem... "  << endl;
	memcpy(Un, U->get_data(0), problem.totalSize * sizeof(double));

	for (int n = 0; n < problem.maxstep; n++)
	{
 	/**** First direction ****/
#pragma omp parallel for
   		for (int j = 1; j < problem.Ny; j++)
   		{
  			double *Bj = B + j*(problem.Nx+1);

			if ((problem.boundCondTypeX0 == 3) || (problem.boundCondTypeX1 == 3))
			{
				cout << "Matrix update... ";
				initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
						   get_alphaX(problem, n), get_betaX(problem), get_gammaX(problem), get_deltaX(problem, n));
				check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
				cout << "Done." << endl;
			}

 			for (int i = 1; i < problem.Nx; i++) // don't forget to handle boundaries
 			{
				Bj[i] = coef_ax*Un[j*(problem.Nx+1)+i+1] +
					coef_cx*Un[j*(problem.Nx+1)+i-1] +
					2.0 * coef_ay*Un[(j+1)*(problem.Nx+1)+i] +
					2.0 * coef_cy*Un[(j-1)*(problem.Nx+1)+i] +
					(1.0 - coef_bx - 2.0*coef_by) * Un[j*(problem.Nx+1)+i] +
					problem.heatSrc(Un,n,i,j);
			}
			Bj[0] = get_b0X(problem, Un, n);
			Bj[problem.Nx] = get_bNX(problem, Un, n);

  			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
									   nrhs, dl1, d1, du1, du21, ipiv1, Bj, ldb));
  		}

 		copyTransposed(temp, B, problem.Nx+1, problem.Ny+1);
		switchPointers(&B, &temp);

		// Second direction
#pragma omp parallel for
  		for (int i = 1; i < problem.Nx; i++)
  		{
 			double *Bi = B + i*(problem.Nx+1);
			if ((problem.boundCondTypeY0 == 3) || (problem.boundCondTypeY1 == 3))
			{
				cout << "Matrix update... ";
				initMatrix(dl2, d2, du2, problem.Ny+1, -coef_ay, 1.0 + coef_by, -coef_cy,
						   get_alphaY(problem, n), get_betaY(problem), get_gammaY(problem),
						   get_deltaY(problem, n));
				check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));
				cout << "Done." << endl;
			}

			for (int j = 1; j < problem.Ny; j++)
			{
			 	Bi[j] += - coef_ay * Un[(j+1)*(problem.Nx + 1) + i] -
					coef_cy * Un[(j-1)*(problem.Nx + 1) + i] +
					coef_by * Un[ j   *(problem.Nx + 1) + i] +
					problem.heatSrc(Un,n,i,j);
			}

			Bi[0] = get_b0Y(problem, Un, n);
			Bi[problem.Ny] = get_bNY(problem, Un, n);

 			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Ny+1,
									   nrhs, dl2, d2, du2, du22, ipiv2, Bi, ldb));
 		}

 		copyTransposed(temp, B, problem.Nx+1, problem.Ny+1);
		switchPointers(&Un, &temp);

		if ((n+1) % problem.numberOfIterationsToSkip == 0)
		{
		 	cout << "Iteration n = " << n << " done." << endl;
		 	memcpy(U->get_data((n+1) / problem.numberOfIterationsToSkip), Un, problem.totalSize * sizeof(double));
		}
	}
	cout << "Done."  << endl;
	cout << "================================================================================" << endl;
	return;
}

Solver::Solver(Problem3d problem)
{
	cout << "Initializiation of the solver... ";

	hx = 1.0 / problem.Nx;
	double coef = 0.5 * problem.tau / (hx*hx);
	coef_ax = coef * problem.Dx;
	coef_cx = coef * problem.Dx;
	coef_bx = coef_ax + coef_cx;

	hy = 1.0 / problem.Ny;
	coef = 0.5 * problem.tau / (hy*hy);
	coef_ay = coef * problem.Dy;
	coef_cy = coef * problem.Dy;
	coef_by = coef_ay + coef_cy;

	hz = 1.0 / problem.Nz;
	coef = 0.5 * problem.tau / (hz*hz);
	coef_az = coef * problem.Dz;
	coef_cz = coef * problem.Dz;
	coef_bz = coef_az + coef_cz;

	ipiv1 = new lapack_int [problem.Nx+1];
	ipiv2 = new lapack_int [problem.Ny+1];
	ipiv3 = new lapack_int [problem.Nz+1];

	nrhs = 1;
	ldb = nrhs;

	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);
	allocateMatrix(&dl2, &d2, &du2, &du22, problem.Ny+1);
	allocateMatrix(&dl3, &d3, &du3, &du23, problem.Nz+1);

	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
			   get_alphaX(problem), get_betaX(problem), get_gammaX(problem), get_deltaX(problem));
	initMatrix(dl2, d2, du2, problem.Ny+1, -coef_ay, 1.0 + coef_by, -coef_cy,
			   get_alphaY(problem), get_betaY(problem), get_gammaY(problem), get_deltaY(problem));
	initMatrix(dl3, d3, du3, problem.Nz+1, -coef_az, 1.0 + coef_bz, -coef_cz,
			   get_alphaZ(problem), get_betaZ(problem), get_gammaZ(problem), get_deltaZ(problem));

	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nz+1, dl3, d3, du3, du23, ipiv3));

	Un = new double [problem.totalSize];  // current time step
	B = new double [problem.totalSize];
	temp = new double [problem.totalSize];
	memset(B, 0, problem.totalSize*sizeof(double));

	cout << "Done."  << endl;
	cout << "================================================================================" << endl;
	return;
}


void Solver::solve(Problem3d problem, Solution *U)
{
	cout << "Solving the problem... "  << endl;
	memcpy(Un, U->get_data(0), problem.totalSize * sizeof(double));

	for (int n = 0; n < problem.maxstep; n++)
	{
		/**** First direction ****/
#pragma omp parallel for
		for (int k = 1; k < problem.Nz; k++)
		{
			for (int j = 1; j < problem.Ny; j++)
			{
				double *Bkj = B + (k*(problem.Ny+1)+ j)*(problem.Nx+1);

				if ((problem.boundCondTypeX0 == 3) || (problem.boundCondTypeX1 == 3))
				{
					cout << "Matrix update... ";
					initMatrix(dl1, d1, du1, problem.Nx+1, -coef_ax, 1.0 + coef_bx, -coef_cx,
							   get_alphaX(problem, n), get_betaX(problem), get_gammaX(problem), get_deltaX(problem, n));
					check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
					cout << "Done." << endl;
				}

				for (int i = 1; i < problem.Nx; i++) // don't forget to handle boundaries
				{
					Bkj[i] = coef_ax*Un[(k    *(problem.Ny+1)+ j  )*(problem.Nx+1)+ i+1] +
						coef_cx*Un[(k    *(problem.Ny+1)+ j  )*(problem.Nx+1)+ i-1] +
						2.0 * coef_ay*Un[(k    *(problem.Ny+1)+ j+1)*(problem.Nx+1)+ i  ] +
						2.0 * coef_cy*Un[(k    *(problem.Ny+1)+ j-1)*(problem.Nx+1)+ i  ] +
						2.0 * coef_az*Un[((k+1)*(problem.Ny+1)+ j  )*(problem.Nx+1)+ i  ] +
						2.0 * coef_cz*Un[((k-1)*(problem.Ny+1)+ j  )*(problem.Nx+1)+ i  ] +
            			(1.0 - coef_ax - coef_cx - 2.0*coef_ay - 2.0*coef_cy -
					 	 2.0*coef_az - 2.0*coef_cz) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i]
						+ problem.heatSrc(Un,n,i,j,k);
				}

				Bkj[0] = get_b0X(problem, Un, n);
				Bkj[problem.Nx] = get_bNX(problem, Un, n);

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
										   nrhs, dl1, d1, du1, du21, ipiv1, Bkj, ldb));
			}
		}

 		copyTransposed3d(temp, B, problem.Nx+1, problem.Ny+1, problem.Nz+1);
		switchPointers(&B, &temp);

		// Second direction
#pragma omp parallel for
		for (int i = 1; i < problem.Nx; i++)
		{
		  	for (int k = 1; k < problem.Nz; k++)
		  	{
		 		double *Bik = B + (i*(problem.Nz+1)+ k)*(problem.Ny+1);
				if ((problem.boundCondTypeY0 == 3) || (problem.boundCondTypeY1 == 3))
				{
					cout << "Matrix update... ";
					initMatrix(dl2, d2, du2, problem.Ny+1, -coef_ay, 1.0 + coef_ay + coef_cy, -coef_cy,
							   get_alphaY(problem, n), get_betaY(problem), get_gammaY(problem),
							   get_deltaY(problem, n));
					check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));
					cout << "Done." << endl;
				}

		 		for (int j = 1; j < problem.Ny; j++)
		 		{
					Bik[j] += - coef_ay*Un[(k    *(problem.Ny+1)+ j+1)*(problem.Nx+1)+ i] -
						coef_cy*Un[(k    *(problem.Ny+1)+ j-1)*(problem.Nx+1)+ i] +
						(coef_ay + coef_cy) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i]
						+ problem.heatSrc(Un,n,i,j,k);
		 		}
				Bik[0] = get_b0Y(problem, Un, n);
				Bik[problem.Ny] = get_bNY(problem, Un, n);

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Ny+1,
										   nrhs, dl2, d2, du2, du22, ipiv2, Bik, ldb));

		 	}
		}
 		copyTransposed3d(temp, B, problem.Ny+1, problem.Nz+1, problem.Nx+1);
		switchPointers(&B, &temp);

        // Third direction
#pragma omp parallel for
		for (int j = 1; j < problem.Ny; j++)
		{
		 	for (int i = 1; i < problem.Nx; i++)
		 	{
		 		double *Bji = B + (j*(problem.Nx+1)+ i)*(problem.Nz+1);
				if ((problem.boundCondTypeZ0 == 3) || (problem.boundCondTypeZ1 == 3))
				{
					cout << "Matrix update... ";
					initMatrix(dl3, d3, du3, problem.Nz+1, -coef_az, 1.0 + coef_az + coef_cz, -coef_cz,
							   get_alphaZ(problem, n), get_betaZ(problem), get_gammaZ(problem),
							   get_deltaZ(problem, n));
					check_error(LAPACKE_dgttrf( (lapack_int) problem.Nz+1, dl3, d3, du3, du23, ipiv3));
					cout << "Done." << endl;
				}

				for (int k = 1; k < problem.Nz; k++)
				{
					Bji[k] += - coef_az*Un[((k+1)*(problem.Ny+1)+ j)*(problem.Nx+1)+ i] -
						coef_cz*Un[((k-1)*(problem.Ny+1)+ j)*(problem.Nx+1)+ i] +
						(coef_az + coef_cz) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i] +
						problem.heatSrc(Un,n,i,j,k);
				}

				Bji[0] = get_b0Z(problem, Un, n);
				Bji[problem.Nz] = get_bNZ(problem, Un, n);

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nz+1,
										   nrhs, dl3, d3, du3, du23, ipiv3, Bji, ldb));

		 	}
		}

 		copyTransposed3d(temp, B, problem.Ny+1, problem.Nz+1, problem.Nx+1);
		switchPointers(&Un, &temp);

		if ((n+1) % problem.numberOfIterationsToSkip == 0)
		{
		 	cout << "Iteration n = " << n << " done." << endl;
		 	memcpy(U->get_data((n+1) / problem.numberOfIterationsToSkip), Un, problem.totalSize * sizeof(double));
		}

	}
	cout << "Done."  << endl;
	cout << "================================================================================" << endl;
	return;
}


Solver::~Solver()
{
	if (ipiv1 != NULL) delete [] ipiv1;
	if (B != NULL) delete [] B;
	if (temp != NULL) delete [] temp;
	if (Un != NULL) delete [] Un;
	if (dl1 != NULL) delete [] dl1;
	if (d1 != NULL) delete [] d1;
	if (du1 != NULL) delete [] du1;
	if (du21 != NULL) delete [] du21;

	if (ipiv2 != NULL) delete [] ipiv2;
	if (dl2 != NULL) delete [] dl2;
	if (d2 != NULL) delete [] d2;
	if (du2 != NULL) delete [] du2;
	if (du22 != NULL) delete [] du22;

	if (ipiv3 != NULL) delete [] ipiv3;
	if (dl3 != NULL) delete [] dl3;
	if (d3 != NULL) delete [] d3;
	if (du3 != NULL) delete [] du3;
	if (du23 != NULL) delete [] du23;


	return;
}

double Solver::get_alphaX(Problem1d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeX0 == 2)
		coef = 1.0 + coef_bx;
	else if (problem.boundCondTypeX0 == 3)
		coef = 1.0 + coef_bx + 2.0*coef_cx*hx*problem.px0(n);
	return coef;
}

double Solver::get_alphaY(Problem2d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeY0 == 2)
		coef = 1.0 + coef_by;
	else if (problem.boundCondTypeY0 == 3)
		coef = 1.0 + coef_by + 2.0*coef_cy*hy*problem.py0(n);
	return coef;
}

double Solver::get_alphaZ(Problem3d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeZ0 == 2)
		coef = 1.0 + coef_bz;
	else if (problem.boundCondTypeZ0 == 3)
		coef = 1.0 + coef_bz + 2.0*coef_cz*hz*problem.pz0(n);
	return coef;
}

double Solver::get_betaX(Problem1d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeX0 > 1)
		coef = -coef_bx;
	return coef;
}

double Solver::get_betaY(Problem2d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeY0 > 1)
		coef = -coef_by;
	return coef;
}

double Solver::get_betaZ(Problem3d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeZ0 > 1)
		coef = -coef_bz;
	return coef;
}

double Solver::get_gammaX(Problem1d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeX1 > 1)
		coef = -coef_bx;
	return coef;
}

double Solver::get_gammaY(Problem2d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeY1 > 1)
		coef = -coef_by;
	return coef;
}

double Solver::get_gammaZ(Problem3d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeZ1 > 1)
		coef = -coef_bz;
	return coef;
}

double Solver::get_deltaX(Problem1d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeX1 == 2)
		coef = 1.0 + coef_bx;
	else if (problem.boundCondTypeX1 == 3)
		coef = 1.0 + coef_bx + 2.0*coef_ax*hx*problem.px1(n);
	return coef;
}

double Solver::get_deltaY(Problem2d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeY1 == 2)
		coef = 1.0 + coef_by;
	else if (problem.boundCondTypeY1 == 3)
		coef = 1.0 + coef_by + 2.0*coef_ay*hy*problem.py1(n);
	return coef;
}

double Solver::get_deltaZ(Problem3d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeZ1 == 2)
		coef = 1.0 + coef_bz;
	else if (problem.boundCondTypeZ1 == 3)
		coef = 1.0 + coef_bz + 2.0*coef_az*hz*problem.pz1(n);
	return coef;
}

double Solver::get_b0X(Problem1d problem, double *Un, int n)
{
	double res = 0.0;
	if (problem.boundCondTypeX0 == 1)
		res = problem.qx0(n);
	else if (problem.boundCondTypeX0 == 2)
		res = coef_bx * Un[1] + (1.0 - coef_bx)* Un[0] -
			2.0*coef_cx*hx*(problem.qx0(n+1) + problem.qx0(n));
	else if (problem.boundCondTypeX0 == 3)
		res = coef_bx * Un[1] + (1.0 - coef_bx - 2.0 * coef_cx* hx * problem.px0(n))* Un[0] -
			2.0 * coef_cx * hx * (problem.qx0(n+1) + problem.qx0(n));
	return res;
}

double Solver::get_b0Y(Problem2d problem, double *Un, int n)
{
	double res = 0.0;
	if (problem.boundCondTypeY0 == 1)
		res = problem.qy0(n);
	else if (problem.boundCondTypeY0 == 2)
		res = coef_by * Un[1] + (1.0 - coef_by)* Un[0] -
			2.0*coef_cy*hy*(problem.qy0(n+1) + problem.qy0(n));
	else if (problem.boundCondTypeY0 == 3)
		res = coef_by * Un[1] + (1.0 - coef_by - 2.0 * coef_cy* hy * problem.py0(n))* Un[0] -
			2.0 * coef_cy * hy * (problem.qy0(n+1) + problem.qy0(n));
	return res;
}

double Solver::get_b0Z(Problem3d problem, double *Un, int n)
{
	double res = 0.0;
	if (problem.boundCondTypeZ0 == 1)
		res = problem.qz0(n);
	else if (problem.boundCondTypeZ0 == 2)
		res = coef_bz * Un[1] + (1.0 - coef_bz)* Un[0] -
			2.0*coef_cz*hz*(problem.qz0(n+1) + problem.qz0(n));
	else if (problem.boundCondTypeZ0 == 3)
		res = coef_bz * Un[1] + (1.0 - coef_bz - 2.0 * coef_cz* hz * problem.pz0(n))* Un[0] -
			2.0 * coef_cz * hz * (problem.qz0(n+1) + problem.qz0(n));
	return res;
}

double Solver::get_bNX(Problem1d problem, double *Un, int n)
{
	double res = 0.0;

	if (problem.boundCondTypeX1 == 1)
		res = problem.qx1(n);
	else if (problem.boundCondTypeX1 == 2)
		res = coef_bx * Un[problem.Nx-1] +
			(1.0 - coef_bx)* Un[problem.Nx] +
			2.0*coef_ax*hx*(problem.qx1(n+1) + problem.qx1(n));
	else if (problem.boundCondTypeX1 == 3)
		res = coef_bx * Un[problem.Nx-1] +
			(1.0 - coef_bx - 2.0 * coef_ax* hx * problem.px1(n))*Un[problem.Nx-1] -
			2.0 * coef_ax * hx * (problem.qx1(n+1) + problem.qx1(n));

	return res;
}

double Solver::get_bNY(Problem2d problem, double *Un, int n)
{
	double res = 0.0;

	if (problem.boundCondTypeY1 == 1)
		res = problem.qy1(n);
	else if (problem.boundCondTypeY1 == 2)
		res = coef_by * Un[problem.Ny-1] +
			(1.0 - coef_by)* Un[problem.Ny] +
			2.0*coef_ay*hy*(problem.qy1(n+1) + problem.qy1(n));
	else if (problem.boundCondTypeY1 == 3)
		res = coef_by * Un[problem.Ny-1] +
			(1.0 - coef_by - 2.0 * coef_ay* hy * problem.py1(n))*Un[problem.Ny-1] -
			2.0 * coef_ay * hy * (problem.qy1(n+1) + problem.qy1(n));

	return res;
}

double Solver::get_bNZ(Problem3d problem, double *Un, int n)
{
	double res = 0.0;

	if (problem.boundCondTypeZ1 == 1)
		res = problem.qz1(n);
	else if (problem.boundCondTypeZ1 == 2)
		res = coef_bz * Un[problem.Nz-1] +
			(1.0 - coef_bz)* Un[problem.Nz] +
			2.0*coef_az*hz*(problem.qz1(n+1) + problem.qz1(n));
	else if (problem.boundCondTypeZ1 == 3)
		res = coef_bz * Un[problem.Nz-1] +
			(1.0 - coef_bz - 2.0 * coef_az* hz * problem.pz1(n))*Un[problem.Ny-1] -
			2.0 * coef_ay * hy * (problem.qy1(n+1) + problem.qy1(n));

	return res;
}

void Solver::allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n)
{
	*dl = new double [n-1];
	*d = new double [n];
	*du = new double [n-1];
	*du2 = new double [n-2];

	return;
}

void Solver::initMatrix(double *dl, double *d, double *du, const int N,
				const double coef_a, const double coef_b, const double coef_c,
				const double coef_b0, const double coef_c0,
				const double coef_aN, const double coef_bN)
{
	for (int i = 1; i < N-1; i++)
		du[i] = coef_c;
	for (int i = 1; i < N-1; i++)
		d[i] = coef_b;
	for (int i = 0; i < N-2; i++)
		dl[i] = coef_a;

	d[0] = coef_b0; // border cond 1.0
	du[0] = coef_c0; // 0.0
	dl[N-2] = coef_aN; // 0.0
	d[N-1] = coef_bN; // border cond 1.0

	return;
}

void Solver::switchPointers(double **A, double **B)
{
	double *tmpPtr = *A;
	*A = *B;
	*B = tmpPtr;
	return;
}

void Solver::copyTransposed(double *B, const double *A, const int N1, const int N2)
{
	for (int j = 0; j < N2; j++)
		for (int i = 0; i < N1; i++)
		   B[i*N2 + j] = A[j*N1 + i];
	return;
}


void Solver::copyTransposed3d(double *B, const double *A, const int N1,
							  const int N2, const int N3)
{
	for (int k = 0; k < N3; k++)
		for (int j = 0; j < N2; j++)
			for (int i = 0; i < N1; i++)
				B[(i*N3+k)*N2 + j] = A[(k*N2 +j)*N1 + i];
	return;
}
