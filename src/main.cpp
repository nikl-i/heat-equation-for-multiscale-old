/*****************************************************************
\todo Testing script should generate solution and compare
+ files by elements and calculate norm
\todo Write good info file
\todo interfaces
\todo Time measurement
\todo Clean makefile
\todo Comment the code
\todo Check w\o OpenBLAS

****************************************************************/

#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "lapacke.h"

const double pi = 3.141592653589793238462643383279502884197;

class Exception: public std::exception
{
public:
    explicit Exception(const char* message):
		msg_(message) {}
    explicit Exception(const std::string& message):
		msg_(message) {}
    virtual ~Exception() throw (){}
    virtual const char* msg() const throw (){
		return msg_.c_str();
    }
protected:
    std::string msg_;
};

struct parameters_t
{
	double maxt, tau;
 	double Dx, q;
	int Nx;
// 	double Dx, Dy, Dz;
//  double qx, qy, qz;
//	int Nx, Ny, Nz;
};


/* Interfaces */
void check_error(const int info);
void print_matrix(char desc, int m, int n, double* a, int lda);
void read_init(const std::string filename_init_cond, const size_t size, double *U);
void read_parameters(const std::string filename_settings, parameters_t *pars);
void read_args(const int argc, char **argv, std::string *filename_settings,
			   std::string *filename_init_cond,	std::string *filename_output);

inline void exists (const std::string& name)
{
	std::ifstream f(name.c_str());
    if (!f.good())
		throw Exception("File doesn't exists:" + name);
}

inline void check_error(const int info)
{
	if (info != 0)
	{
		std::stringstream ss;
		if (info < 0)
			ss << "Lapack: the " << -info << "-th argument had an illegal value";
		else
			ss << "Lapack: U(" << info << "," << info << ") is exactly zero.";
		throw Exception(ss.str());
	}
}

int main(int argc, char **argv)
{
	double *U = NULL, *B1 = NULL, *d = NULL, *du = NULL, *du2 = NULL, *dl = NULL;
	lapack_int *ipiv = NULL; // for Lapack
	try
	{
		double h, h2, coef, coef_a, coef_b, coef_c;
		int maxstep, n, i;

		lapack_int ldb, nrhs; // for Lapack

//	double start, finish; // timer
		parameters_t pars;
		std::string filename_settings, filename_init_cond, filename_output;

		/* Get parameters */
		read_args(argc, argv, &filename_settings, &filename_init_cond, &filename_output);
		read_parameters(filename_settings, &pars);

		//maxstep = (int) pars.maxt/pars.tau;
		maxstep = 1;
		h = 1.0 / pars.Nx;
		h2 = h*h;

		/* Allocation */
		int totalSize = pars.Nx+1;

		U = (double *) malloc(totalSize * sizeof(double));
		B1 = (double *) malloc(totalSize * sizeof(double));

		ipiv = (lapack_int *) malloc(pars.Nx * sizeof(lapack_int));
		dl = (double *) malloc(pars.Nx * sizeof(double));
		d = (double *) malloc((pars.Nx+1) * sizeof(double));
		du = (double *) malloc(pars.Nx * sizeof(double));
		du2 = (double *) malloc((pars.Nx-1) * sizeof(double));

		/* init */
		read_init(filename_init_cond, pars.Nx, U);
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

		/* Solution */
		nrhs = 1;
		ldb = nrhs;

		check_error(LAPACKE_dgttrf( (lapack_int) pars.Nx+1, dl, d, du, du2, ipiv));

		for (n = 1; n < maxstep; n++)
		{
			B1[0] = 0.0;
			for (i = 1; i < pars.Nx; i++)
				B1[i] = coef_a*U[i-1] + coef_c*U[i+1] + (1.0 - coef_b)*U[i];
			B1[pars.Nx] = 0.0;

			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) pars.Nx+1,
									   nrhs, dl, d, du, du2, ipiv, B1, ldb));

			for (i = 0; i < pars.Nx+1; i++)
				U[i] = B1[i];

			print_matrix('B', pars.Nx+1, nrhs, B1, ldb);
		}
		std::cout << "Done." << std::endl;
	}
	catch(Exception error)
	{
		std::cerr << error.msg() << std::endl;
		std::cerr << "Error. Program was terminated." << std::endl;
	}

	/* Deallocation */
	if (U != NULL) free(U);
	if (B1 != NULL) free(B1);

//	if (ipiv != NULL) free(ipiv);*/
	if (dl != NULL) free(dl);
	if (d != NULL) free(d);
	if (du != NULL) free(du);
	if (du2 != NULL) free(du2);

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

void read_args(const int argc, char **argv, std::string *filename_settings,
			   std::string *filename_init_cond,	std::string *filename_output)
{
	for (int i = 1; i < argc; i += 2)
	{
		if (!strcmp(argv[i], "-s"))
			*filename_settings = argv[i+1];
		else if (!strcmp(argv[i], "-i"))
			*filename_init_cond = argv[i+1];
		else if (!strcmp(argv[i], "-o"))
			*filename_output = argv[i+1];
		else
			throw Exception("Unknown option!");
	}

	if ((*filename_settings).empty())
		*filename_settings = "settings";
	if ((*filename_init_cond).empty())
		*filename_init_cond = "u0.bin";
	if ((*filename_output).empty())
		*filename_output = "u.bin";

	exists(*filename_settings);
	exists(*filename_init_cond);

	std::cout << "================ Files: ==================" << std::endl;
	std::cout << "init cond file: " << *filename_init_cond << std::endl;
	std::cout << "settings file : " << *filename_settings << std::endl;
	std::cout << "output file   : " << *filename_output << std::endl;
	std::cout << "==========================================" << std::endl;

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
			{
				std::transform(token.begin(), token.end(), token.begin(), ::tolower);
				keyname = token;
			}
			else
			{
				double value = atof(token.c_str());

				if (keyname == "tau")
					(*pars).tau = value;
				else if (keyname == "maxt")
					(*pars).maxt = value;
				else if (keyname ==  "dx")
					(*pars).Dx = value;
				else if (keyname == "q")
					(*pars).q = value;
				else if (keyname == "nx")
					(*pars).Nx = (int) value;
				else
					throw(Exception("Unknown parameter name: " + token));

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
	std::ifstream input(filename_init_cond.c_str(), std::ios::binary|std::ios::in);
	input.read((char *) U, size * sizeof(double));

	for (unsigned int i = 0; i < size+1; i++)
		std::cout << i << ":" << U[i] << std::endl;
	return;
}
