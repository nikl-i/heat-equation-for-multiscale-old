#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <sstream>

#include "lapacke.h"

const double pi = 3.141592653589793238462643383279502884197;

void check_error(const int info);
void print_matrix(char desc, int m, int n, double* a, int lda);

struct parameters_t
{
 	double Dx, q, maxt, tau;
	int Nx;
};

void read_args(const int argc, char **argv, std::string *filename_settings,
			   std::string *filename_init_cond)
{
	if (argc == 3)
	{
		*filename_settings = argv[1];
		*filename_init_cond = argv[2];
	}
	return;
}

void read_parameters(const std::string filename_settings, parameters_t *pars)
{

	std::string line;
	std::ifstream file(filename_settings.c_str());

	while(std::getline(file, line))
	{
		std::string token;
		std::string keyname;
		std::istringstream tokens(line) ;
		while(tokens >> token)
		{
			if (token == "#") break;
			if (token == "=") continue;

			if (keyname.empty())
				keyname = token;
			else
			{
				double value = atof(token.c_str());

				if (keyname == "tau")	{
					(*pars).tau = value;
				} else if (keyname == "maxt") {
					(*pars).maxt = value;
				} else if (keyname ==  "Dx") {
					(*pars).Dx = value;
				} else if (keyname == "q") {
					(*pars).q = value;
				} else if (keyname == "Nx") {
					(*pars).Nx = (int) value;
				}
				keyname.clear();
			}
		}
	}
	std::cout << "================ Parameters: =============" << std::endl;
	std::cout << "tau  = " << (*pars).tau << std::endl;
	std::cout << "maxt = " << (*pars).maxt << std::endl;
	std::cout << "Nx   = " << (*pars).Nx << std::endl;
	std::cout << "Dx   = " << (*pars).Dx << std::endl;
	std::cout << "q    = " << (*pars).q << std::endl;
	std::cout << "==========================================" << std::endl;

	return;
}


void read_init(const std::string filename_init_cond, const size_t size, double *U)
{
	std::ifstream input( filename_init_cond.c_str(), std::ios::binary|std::ios::in );

	input.read((char *) &U, size * sizeof(double));


	for (unsigned int i = 0; i < size; i++)
		std::cout << U[i];

	return;
}


int main(int argc, char **argv)
{
	double *U, *B1, *d, *du, *du2, *dl;
	double h, h2, coef, coef_a, coef_b, coef_c;
	int maxstep, n, i;

	lapack_int *ipiv; // for Lapack
	lapack_int ldb, nrhs, info; // for Lapack

//	double start, finish; // timer
	parameters_t pars;
	std::string filename_settings, filename_init_cond;

	read_args(argc, argv, &filename_settings, &filename_init_cond);


	std::cout << "init cond file: " << filename_init_cond << std::endl;
	std::cout << "settings  file: " << filename_settings << std::endl;

	read_parameters(filename_settings, &pars);


	maxstep = (int) pars.maxt/pars.tau;
	h = 1.0 / pars.Nx;
	h2 = h*h;

	U = (double *) malloc((pars.Nx+1) * sizeof(double));
	read_init(filename_init_cond, pars.Nx, U);




	return 0;

	B1 = (double *) malloc((pars.Nx+1) * sizeof(double));
	ipiv = (lapack_int *) malloc((pars.Nx+1) * sizeof(lapack_int));

	dl = (double *) malloc(pars.Nx * sizeof(double));
	d = (double *) malloc((pars.Nx+1) * sizeof(double));
	du = (double *) malloc(pars.Nx * sizeof(double));
	du2 = (double *) malloc((pars.Nx-1) * sizeof(double));

	coef = 0.5 * pars.tau / h2;
	coef_a = - coef * pars.Dx;
	coef_c = - coef * pars.Dx;
	coef_b = 1.0 + coef_a + coef_c + 0.5 * pars.tau * pars.q;

	du[0] = 0.0;
	for (int i = 1; i < pars.Nx; i++)
		du[i] = coef_c;

	d[0] = 1.0; // border cond
	for (int i = 1; i < pars.Nx; i++)
		d[i] = coef_b;
	d[pars.Nx] = 1.0; // border cond

	for (int i = 0; i < pars.Nx-1; i++)
		dl[i] = coef_a;
	dl[pars.Nx-1] = 0.0;


	nrhs = 1;
	ldb = nrhs;

	info = LAPACKE_dgttrf( (lapack_int) pars.Nx+1, dl, d, du, du2, ipiv);
	check_error(info);

	for (n = 1; n< maxstep; n++)
	{

		B1[0] = 0.0;
		for (i = 1; i < pars.Nx; i++)
			B1[i] = coef_a*U[i-1] + coef_c*U[i+1] + (1.0 - coef_b)*U[i];
		B1[pars.Nx] = 0.0;

		info = LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) pars.Nx+1,
						  nrhs, dl, d, du, du2, ipiv, B1, ldb);
		check_error(info);

		for (i = 0; i < pars.Nx+1; i++)
			U[i] = B1[i];

		print_matrix('B', pars.Nx+1, nrhs, B1, ldb);
	}

	return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix(char desc, int m, int n, double* a, int lda)
{
	int i, j;
	printf("\n----%c:----\n", desc);
	for( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %6.2f", a[j+i*lda] );
		printf( "\n" );
	}
	return;
}

void check_error(const int info)
{
	if (info == 0) return;
	if (info < 0)
		printf("Error: the %d-th argument had an illegal value", -info);
	else
		printf("Error: U(%d,%d) is exactly zero.", info, info);
    exit(1);
}
