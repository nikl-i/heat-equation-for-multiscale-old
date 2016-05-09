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
	ipiv1 = (lapack_int *) malloc((problem.Nx+1) * sizeof(lapack_int));
	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);
	Un = (double *) malloc(problem.totalSize * sizeof(double));  // current time step
	B = (double *) malloc(problem.totalSize * sizeof(double));

	cout << "Initializiation of the solver...";
	h1 = 1.0 / problem.Nx;
	coef1 = 0.5 * problem.tau / (h1*h1);
	coef_a1 = coef1 * problem.Dx;
	coef_c1 = coef1 * problem.Dx;
	coef_b1 = coef_a1 + coef_c1;
	nrhs = 1;
	ldb = nrhs;
	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_a1, 1.0 + coef_b1, -coef_c1,
			   get_d1(problem), get_d2(problem), get_d3(problem), get_d4(problem));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	cout << "Done."  << endl;

	return;
}

Solver::Solver(Problem2d problem)
{
	cout << "Initializiation of the solver...";
	h1 = 1.0 / problem.Nx;
	coef1 = 0.5 * problem.tau / (h1*h1);
	coef_a1 = coef1 * problem.Dx;
	coef_c1 = coef1 * problem.Dx;
	coef_b1 = 1.0 + coef_a1 + coef_c1;

	h2 = 1.0 / problem.Ny;
	coef2 = 0.5 * problem.tau / (h2*h2);
	coef_a2 = coef2 * problem.Dy;
	coef_c2 = coef2 * problem.Dy;
	coef_b2 = 1.0 + coef_a2 + coef_c2;

	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);
	allocateMatrix(&dl2, &d2, &du2, &du22, problem.Ny+1);

	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_a1, 1.0 + coef_a1 + coef_c1,
			   -coef_c1, get_d1(problem), get_d2(problem), get_d3(problem), get_d4(problem));
	initMatrix(dl2, d2, du2, problem.Ny+1, -coef_a2, 1.0 + coef_a2 + coef_c2,
			   -coef_c2, get_d1(problem), get_d2(problem), get_d3(problem), get_d4(problem));

	ipiv1 = (lapack_int *) malloc((problem.Nx+1) * sizeof(lapack_int));
	ipiv2 = (lapack_int *) malloc((problem.Ny+1) * sizeof(lapack_int));
	nrhs = 1;
	ldb = nrhs;

	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));

	B = (double *) malloc(problem.totalSize * sizeof(double));
	Un = (double *) malloc(problem.totalSize * sizeof(double));  // current time step
	cout << "Done."  << endl;

	return;
}

Solver::Solver(Problem3d problem)
{
	cout << "Initializiation of the solver...";

	h1 = 1.0 / problem.Nx;
	coef1 = 0.5 * problem.tau / (h1*h1);
	coef_a1 = coef1 * problem.Dx;
	coef_c1 = coef1 * problem.Dx;
	coef_b1 = 1.0 + coef_a1 + coef_c1;

	h2 = 1.0 / problem.Ny;
	coef2 = 0.5 * problem.tau / (h2*h2);
	coef_a2 = coef2 * problem.Dy;
	coef_c2 = coef2 * problem.Dy;
	coef_b2 = 1.0 + coef_a2 + coef_c2;

	h3 = 1.0 / problem.Nz;
	coef3 = 0.5 * problem.tau / (h3*h3);
	coef_a3 = coef3 * problem.Dz;
	coef_c3 = coef3 * problem.Dz;
	coef_b3 = 1.0 + coef_a3 + coef_c3;

	B = (double *) malloc(problem.totalSize * sizeof(double));

	allocateMatrix(&dl1, &d1, &du1, &du21, problem.Nx+1);
	allocateMatrix(&dl2, &d2, &du2, &du22, problem.Ny+1);
	allocateMatrix(&dl3, &d3, &du3, &du23, problem.Nz+1);

	initMatrix(dl1, d1, du1, problem.Nx+1, -coef_a1, 1.0 + coef_a1 + coef_c1, -coef_c1, 1.0, 0.0, 0.0, 1.0);
	initMatrix(dl2, d2, du2, problem.Ny+1, -coef_a2, 1.0 + coef_a2 + coef_c2, -coef_c2, 1.0, 0.0, 0.0, 1.0);
	initMatrix(dl3, d3, du3, problem.Nz+1, -coef_a3, 1.0 + coef_a3 + coef_c3, -coef_c3, 1.0, 0.0, 0.0, 1.0);

	ipiv1 = (lapack_int *) malloc((problem.Nx+1) * sizeof(lapack_int));
	ipiv2 = (lapack_int *) malloc((problem.Ny+1) * sizeof(lapack_int));
	ipiv3 = (lapack_int *) malloc((problem.Nz+1) * sizeof(lapack_int));

	nrhs = 1;
	ldb = nrhs;

	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Ny+1, dl2, d2, du2, du22, ipiv2));
	check_error(LAPACKE_dgttrf( (lapack_int) problem.Nz+1, dl3, d3, du3, du23, ipiv3));

	cout << "Done."  << endl;

	return;
}


Solver::~Solver()
{
	if (ipiv1 != NULL) free(ipiv1);
	if (B != NULL) free(B);
	if (Un != NULL) free(Un);
	if (dl1 != NULL) free(dl1);
	if (d1 != NULL) free(d1);
	if (du1 != NULL) free(du1);
	if (du21 != NULL) free(du21);

	if (ipiv2 != NULL) free(ipiv2);
	if (dl2 != NULL) free(dl2);
	if (d2 != NULL) free(d2);
	if (du2 != NULL) free(du2);
	if (du22 != NULL) free(du22);

	if (ipiv3 != NULL) free(ipiv3);
	if (dl3 != NULL) free(dl3);
	if (d3 != NULL) free(d3);
	if (du3 != NULL) free(du3);
	if (du23 != NULL) free(du23);

	return;
}

double Solver::get_d1(Problem1d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeX0 == 2)
		coef = 1.0 + coef_b1;
	else if (problem.boundCondTypeX0 == 3)
		coef = 1.0 + coef_b1 + 2.0*coef_c1*h1*problem.px0(n);
	return coef;
}

double Solver::get_d2(Problem1d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeX0 > 1)
		coef = -coef_b1;
	return coef;
}

double Solver::get_d3(Problem1d problem)
{
	double coef = 0.0;
	if (problem.boundCondTypeX1 > 1)
		coef = -coef_b1;
	return coef;
}

double Solver::get_d4(Problem1d problem, int n)
{
	double coef = 1.0;
	if (problem.boundCondTypeX1 == 2)
		coef = 1.0 + coef_b1;
	else if (problem.boundCondTypeX1 == 3)
		coef = 1.0 + coef_b1 + 2.0*coef_a1*h1*problem.px1(n);
	return coef;
}

double Solver::get_b0(Problem1d problem, double *Un, int n)
{
	double res = 0.0;
	if (problem.boundCondTypeX0 == 1)
		res = problem.qx0(n);
	else if (problem.boundCondTypeX0 == 2)
		res = coef_b1 * Un[1] + (1.0 - coef_b1)* Un[0] -
			2.0*coef_c1*h1*(problem.qx0(n+1) + problem.qx0(n));
	else if (problem.boundCondTypeX0 == 3)
		res = coef_b1 * Un[1] + (1.0 - coef_b1 - 2.0 * coef_c1* h1 * problem.px0(n))* Un[0] -
			2.0 * coef_c1 * h1 * (problem.qx0(n+1) + problem.qx0(n));
	return res;
}

double Solver::get_bN(Problem1d problem, double *Un, int n)
{
	double res = 0.0;
	if (problem.boundCondTypeX1 == 1)
		res = problem.qx1(n);
	else if (problem.boundCondTypeX1 == 2)
		res = coef_b1 * Un[problem.Nx-1] +
			(1.0 - coef_b1)* Un[problem.Nx] +
			2.0*coef_a1*h1*(problem.qx1(n+1) + problem.qx1(n));
	else if (problem.boundCondTypeX1 == 3)
		res = coef_b1 * Un[problem.Nx-1] +
			(1.0 - coef_b1 - 2.0 * coef_a1* h1 * problem.px1(n))*Un[problem.Nx-1] -
			2.0 * coef_a1 * h1 * (problem.qx1(n+1) + problem.qx1(n));
	return res;
}

void Solver::solve(Problem1d problem, Solution *U)
{
	cout << "================================================" << endl;
	cout << "a_x = " << coef_a1 << " | b_x = " << coef_b1 << " | c_x = " << coef_c1  << endl;
	cout << "================================================" << endl;
    /*************************************************************/
    /******************** Main loop ******************************/
    /*************************************************************/
	memcpy(Un, U->U[0], problem.totalSize * sizeof(double));

	for (int n = 0; n < problem.maxstep; n++)
	{
		if ((problem.boundCondTypeX0 == 3) || (problem.boundCondTypeX1 == 3))
		{
			cout << "Matrix update... ";
			initMatrix(dl1, d1, du1, problem.Nx+1, -coef_a1, 1.0 + coef_b1, -coef_c1,
					   get_d1(problem, n), coef_d2, coef_d3, get_d4(problem, n));
			check_error(LAPACKE_dgttrf( (lapack_int) problem.Nx+1, dl1, d1, du1, du21, ipiv1));
			cout << "Done." << endl;
		}

		for (int i = 1; i < problem.Nx; i++)
			B[i] = coef_a1*Un[i-1] + coef_c1*Un[i+1] + (1.0 - coef_b1)*Un[i];
		B[0] = get_b0(problem, Un, n);
		B[problem.Nx] = get_bN(problem, Un, n);

		check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
								   nrhs, dl1, d1, du1, du21, ipiv1, B, ldb));
		temp = Un;
		Un = B;
		B = temp;

		if ((n+1) % problem.numberOfIterationsToSkip == 0)
		{
			cout << "Iteration t=" << n << "|" << (n+1) / problem.numberOfIterationsToSkip <<" done." << endl;
			memcpy(U->U[(n+1) / problem.numberOfIterationsToSkip ], Un, problem.totalSize * sizeof(double *));
		}
	}
    /*************************************************************/
    /******************** End of main loop ***********************/
    /*************************************************************/
	return;
}

void Solver::solve(Problem2d problem, Solution *U)
{

	cout << "================================================" << endl;
	cout << "a_x = " << coef_a1 << " | b_x = " << coef_b1 << " | c_x = " << coef_c1  << endl;
	cout << "a_y = " << coef_a2 << " | b_y = " << coef_b2 << " | c_y = " << coef_c2  << endl;
	cout << "================================================" << endl;

	/*************************************************************/
	/******************** Main loop ******************************/
	/*************************************************************/
	for (int n = 0; n < problem.maxstep; n++)
	{
		double *Un = U->U[n];
		/**** First direction ****/
		for (int j = 1; j < problem.Ny; j++)
		{
			double *Bj = B + j*(problem.Nx+1);
			for (int i = 1; i < problem.Nx; i++) // don't forget to handle boundaries
			{
				Bj[i] = coef_a1*Un[j*(problem.Nx+1)+i+1] +
					coef_c1*Un[j*(problem.Nx+1)+i-1] +
					2.0 * coef_a2*Un[(j+1)*(problem.Nx+1)+i] +
					2.0 * coef_c2*Un[(j-1)*(problem.Nx+1)+i] +
					(1.0 - coef_a1 - coef_c1 - 2.0*coef_a2 - 2.0*coef_c2) * Un[j*(problem.Nx+1)+i];
			}

			Bj[0] = 0.0;
			Bj[problem.Nx] = 0.0;

			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
		 							   nrhs, dl1, d1, du1, du21, ipiv1, Bj, ldb));
		}

 		assignTransposed(&B, B, problem.Nx+1, problem.Ny+1);

        // Second direction
		for (int i = 1; i < problem.Nx; i++)
		{
			double *Bi = B + i*(problem.Nx+1);
			for (int j = 1; j < problem.Ny; j++)
			{
			 	Bi[j] += - coef_a2 * Un[(j+1)*(problem.Nx + 1) + i] -
					coef_c2 * Un[(j-1)*(problem.Nx + 1) + i] +
					(coef_a2 + coef_c2)* Un[ j   *(problem.Nx + 1) + i];
			}

			Bi[0] = 0.0;
			Bi[problem.Ny] = 0.0;

			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Ny+1,
									   nrhs, dl2, d2, du2, du22, ipiv2, Bi, ldb));

			assignTransposed(&(U->U[n+1]), B, problem.Ny+1, problem.Nx+1); // look for not allocate U[n] at the begining
		}

		if (n % problem.numberOfIterationsToSkip == 0)
			cout << "Iteration t=" << n*problem.tau << " done." << endl;
	}
	/*************************************************************/
	/******************** End of main loop ***********************/
	/*************************************************************/
	return;
}

void Solver::solve(Problem3d problem, Solution *U)
{
	cout << "================================================" << endl;
	cout << "a_x = " << coef_a1 << " | b_x = " << coef_b1 << " | c_x = " << coef_c1  << endl;
	cout << "a_y = " << coef_a2 << " | b_y = " << coef_b2 << " | c_y = " << coef_c2  << endl;
	cout << "a_z = " << coef_a3 << " | b_z = " << coef_b3 << " | c_z = " << coef_c3  << endl;
	cout << "================================================" << endl;
	/*************************************************************/
	/******************** Main loop ******************************/
	/*************************************************************/
	for (int n = 0; n < problem.maxstep; n++)
	{
		double *Un = U->U[n];
		/**** First direction ****/
		for (int k = 1; k < problem.Nz; k++)
		{
			for (int j = 1; j < problem.Ny; j++)
			{
				double *Bkj = B + (k*(problem.Ny+1)+ j)*(problem.Nx+1);
				for (int i = 1; i < problem.Nx; i++) // don't forget to handle boundaries
				{
					Bkj[i] = coef_a1*Un[(k    *(problem.Ny+1)+ j  )*(problem.Nx+1)+ i+1] +
						coef_c1*Un[(k    *(problem.Ny+1)+ j  )*(problem.Nx+1)+ i-1] +
						2.0 * coef_a2*Un[(k    *(problem.Ny+1)+ j+1)*(problem.Nx+1)+ i  ] +
						2.0 * coef_c2*Un[(k    *(problem.Ny+1)+ j-1)*(problem.Nx+1)+ i  ] +
						2.0 * coef_a3*Un[((k+1)*(problem.Ny+1)+ j  )*(problem.Nx+1)+ i  ] +
						2.0 * coef_c3*Un[((k-1)*(problem.Ny+1)+ j  )*(problem.Nx+1)+ i  ] +
            			(1.0 - coef_a1 - coef_c1 - 2.0*coef_a2 - 2.0*coef_c2 -
					 	 2.0*coef_a3 - 2.0*coef_c3) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i];
				}

				Bkj[0] = 0.0;
				Bkj[problem.Nx] = 0.0;

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nx+1,
										   nrhs, dl1, d1, du1, du21, ipiv1, Bkj, ldb));

			}
		}
	   	assignTransposed3d(&B, B, problem.Nx+1, problem.Ny+1, problem.Nz+1);

		// Second direction
		for (int i = 1; i < problem.Nx; i++)
		{
		  	for (int k = 1; k < problem.Nz; k++)
		  	{
		 		double *Bik = B + (i*(problem.Nz+1)+ k)*(problem.Ny+1);
		 		for (int j = 1; j < problem.Ny; j++)
		 		{
					Bik[j] += - coef_a2*Un[(k    *(problem.Ny+1)+ j+1)*(problem.Nx+1)+ i] -
						coef_c2*Un[(k    *(problem.Ny+1)+ j-1)*(problem.Nx+1)+ i] +
						(coef_a2 + coef_c2) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i];
		 		}
				Bik[0] = 0.0;
				Bik[problem.Ny] = 0.0;

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Ny+1,
										   nrhs, dl2, d2, du2, du22, ipiv2, Bik, ldb));

		 	}
		}
		assignTransposed3d(&B, B, problem.Ny+1, problem.Nz+1, problem.Nx+1);

        // Third direction
		for (int j = 1; j < problem.Ny; j++)
		{
		 	for (int i = 1; i < problem.Nx; i++)
		 	{
		 		double *Bji = B + (j*(problem.Nx+1)+ i)*(problem.Nz+1);
				for (int k = 1; k < problem.Nz; k++)
				{
					Bji[k] += - coef_a3*Un[((k+1)*(problem.Ny+1)+ j)*(problem.Nx+1)+ i] -
						coef_c3*Un[((k-1)*(problem.Ny+1)+ j)*(problem.Nx+1)+ i] +
						(coef_a3 + coef_c3) * Un[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i];
				}

				Bji[0] = 0.0;
				Bji[problem.Nz] = 0.0;

				check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) problem.Nz+1,
										   nrhs, dl3, d3, du3, du23, ipiv3, Bji, ldb));

		 	}
		}
		assignTransposed3d(&(U->U[n+1]), B, problem.Nz+1, problem.Nx+1, problem.Ny+1);
		// for (int k = 0; k < problem.Nz+1; k++)
		// 	for (int j = 0; j < problem.Ny+1; j++)
		// 		for (int i = 0; i < problem.Nx+1; i++)
		// 			U[n+1][(k*(problem.Ny+1)+j)*(problem.Nx+1)+i] = B[(k*(problem.Ny+1)+j)*(problem.Nx+1)+i];

		// cout << U[1][50*(problem.Nx+1)+50] << endl;

		if (n % problem.numberOfIterationsToSkip == 0)
			cout << "Iteration t=" << n*problem.tau << " done." << endl;

	}

	/*************************************************************/
	/******************** End of main loop ***********************/
	/*************************************************************/

	return;
}

void Solver::allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n)
{
	*dl = (double *) malloc((n-1) * sizeof(double));
	*d = (double *) malloc(n * sizeof(double));
	*du = (double *) malloc((n-1) * sizeof(double));
	*du2 = (double *) malloc((n-2) * sizeof(double));

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

void Solver::assignTransposed(double **B, double *A, int N1, int N2)
{
	double *temp;
	temp = (double *) malloc(N1*N2*sizeof(double));
	for (int j = 0; j < N2; j++)
		for (int i = 0; i < N1; i++)
			temp[i*N2 + j] = A[j*N1 + i];
	if (*B != NULL) free(*B);
	*B = temp;
	return;
}

void Solver::assignTransposed3d(double **B, double *A, int N1, int N2, int N3)
{
	double *temp;
	temp = (double *) malloc(N1*N2*N3*sizeof(double));

	for (int k = 0; k < N3; k++)
		for (int j = 0; j < N2; j++)
			for (int i = 0; i < N1; i++)
				temp[(i*N3+k)*N2 + j] = A[(k*N2 +j)*N1 + i];
	if (*B != NULL) free(*B);
	*B = temp;
	return;
}
