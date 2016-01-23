/*****************************************************************
\todo Add __FUNCNAME__ to Exception
\todo 1d 2d and 3d versions
\todo cpu and gpu versions
\todo batched execultion for 2d and 3d

\todo Time measurement
\todo Check w\o OpenBLAS
\todo Rethink plotting script
\todo Comment the code
****************************************************************/
#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "lapacke.h"
#define ITERATION_TO_OUTPUT 100

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
 	double Dx, Qx;
//  int numOfDims =
	int Nx;
// 	double Dx, Dy, Dz;
//  double Qx, Qy, Qz;
//	int Nx, Ny, Nz;
// 	double Lx, Ly, Lz;
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




void allocateMatrix(double **dl, double **d, double **du, double **du2,
				const int n)
{
	*dl = (double *) malloc(pars.Nx * sizeof(double));
	*d = (double *) malloc((pars.Nx+1) * sizeof(double));
	*du = (double *) malloc(pars.Nx * sizeof(double));
	*du2 = (double *) malloc((pars.Nx-1) * sizeof(double));



	return;
}

void initMatrix(double *dl, double *d, double *du, const int N,
			const double tau, const double D)
{
	double h, coef, coef_a, coef_b, coef_c;

	h = 1.0 / N;
	coef = - 0.5 * tau / (h*h);

	coef_a = coef * D;
	coef_c = coef * D;
	coef_b = 1.0 + coef_a + coef_c;

	du[0] = 0.0;
	for (int i = 1; i < N; i++)
		du[i] = coef_c;

	d[0] = 1.0; // border cond
	for (int i = 1; i < N; i++)
		d[i] = coef_b;
	d[N] = 1.0; // border cond

	for (int i = 0; i < N-1; i++)
		dl[i] = coef_a;
	dl[N-1] = 0.0;

	return;
}

void Output(string filename_output, int maxstep, double **U)
{
		for (int n = 0; n <= maxstep; n++)
		{
			std::ofstream fileOut (filename_output.c_str(),
								   std::ios::binary | std::ios::out | std::ios::trunc);
			fileOut.write((char*) U[n], totalSize*sizeof(double));
		}
}



int main(int argc, char **argv)
{
	double *U = NULL, *B = NULL;
	double *d = NULL, *du = NULL, *du2 = NULL, *dl = NULL;
//	double *d1 = NULL, *du1 = NULL, *du21 = NULL, *dl1 = NULL;
//	double *d2 = NULL, *du2 = NULL, *du22 = NULL, *dl2 = NULL;
//	double *d3 = NULL, *du3 = NULL, *du23 = NULL, *dl3 = NULL;

	try
	{
		double h, coef, coef_a, coef_b, coef_c;
//		double h1, coef1, coef_a1, coef_b1, coef_c1;
//		double h2, coef2, coef_a2, coef_b2, coef_c2;
//		double h3, coef3, coef_a3, coef_b3, coef_c3;

		int maxstep;
		parameters_t pars;
		std::string filename_settings, filename_init_cond, filename_output;

		/* Get parameters */
		read_args(argc, argv, &filename_settings, &filename_init_cond, &filename_output);
		read_parameters(filename_settings, &pars);

		maxstep = pars.maxt/pars.tau;

		/* Allocate1d */
		int totalSize = pars.Nx + 1;
		if (DIM > 1)
			totalSize *= pars.Ny + 1;
		if (DIM > 2)
			totalSize *= pars.Nz + 1;

		allocateData(U, B, totalSize, maxstep);
		B = (double *) malloc(totalSize * sizeof(double));
		U = (double **) malloc((maxstep +1) * sizeof(double *));
		for (int n = 0; i <= maxstep; i++)
			U[n] = (double *) malloc(totalSize * sizeof(double));

		allocateMatrix(&dl, &d, &du, &du2, pars.Nx);
		// allocateMatrix(dl1, d1, du1, du21, pars.Nx);
		// if (DIM > 1)
		// 	allocate1d(dl2, d2, du2, du22, pars.Ny);
		// if (DIM > 2)
		// 	allocate1d(dl3, d3, du3, du23, pars.Nz);


		/* initialize1d(filename_init_cond, pars, dl, d, du, du2, U) */
		read_init(filename_init_cond, pars.Nx, U);

		/****init*/
		initMatrix(dl, d, du, pars.Nx, pars.tau, pars.Dx);

		// initMatrix(dl1, d1, du1, du21, pars.Nx);
		// if (DIM > 1)
		// 	initMatrix(dl2, d2, du2, du22, pars.Ny);
		// if (DIM > 2)
		// 	initMatrix(dl3, d3, du3, du23, pars.Nz);


		/* Solver1d_init_cpu(U, pars, dl, d, du, du2) */
		lapack_int ldb, nrhs; // for Lapack
		lapack_int *ipiv = NULL; // for Lapack
//  	lapack_int *ipiv1 = NULL; // for Lapack
//	    lapack_int *ipiv2 = NULL; // for Lapack
//	    lapack_int *ipiv3 = NULL; // for Lapack
		ipiv = (lapack_int *) malloc(pars.Nx * sizeof(lapack_int));

		nrhs = 1;
		ldb = nrhs;
		check_error(LAPACKE_dgttrf( (lapack_int) pars.Nx+1, dl, d, du, du2, ipiv));

		/*************************************************************/
		/******************** Main loop ******************************/
		/*************************************************************/
		for (int n = 0; n < maxstep; n++)
		{
			B[0] = 0.0;
			double *Un = U[n]; // current time step
			for (int i = 1; i < pars.Nx; i++)
				B[i] = - coef_a*Un[i-1] - coef_c*Un[i+1] + (1.0 + coef_a + coef_c)*Un[i];
			B[pars.Nx] = 0.0;

			check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) pars.Nx+1,
									   nrhs, dl, d, du, du2, ipiv, B1, ldb));

			for (int i = 0; i < pars.Nx+1; i++)
				U[n+1][i] = B1[i];

			if (n % ITERATION_TO_OUTPUT == 0)
				std::cout << "Iteration t=" << n*pars.tau << " done." << std::endl;
		}
		/*************************************************************/
		/******************** End of main loop ***********************/
		/*************************************************************/

//	if (ipiv != NULL) free(ipiv);*/

		Output(filename_output, maxstep, U);
		std::cout << "Done." << std::endl;

	}
	catch(Exception error)
	{
		std::cerr << error.msg() << std::endl;
		std::cerr << "Error. Program was terminated." << std::endl;
	}

	/* Deallocation */
	if (U != NULL) free(U);
	if (B != NULL) free(B);

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
					pars->tau = value;
				else if (keyname == "maxt")
					pars->maxt = value;
				else if (keyname ==  "dx")
					pars->Dx = value;
				else if (keyname == "q")
					pars->Qx = value;
				else if (keyname == "nx")
					pars->Nx = (int) value;
				else
					throw(Exception("Unknown parameter name: " + token));

				keyname.clear();
			}
		}
	}

	std::cout << "================ Parameters: =============" << std::endl;
	std::cout << "tau  = " << pars->tau << std::endl;
	std::cout << "maxt = " << pars->maxt << std::endl;
	std::cout << "Nx   = " << pars->Nx << std::endl;
	std::cout << "Dx   = " << pars->Dx << std::endl;
	std::cout << "Qx    = " << pars->Qx << std::endl;

	std::cout << "h    = " << 1.0 / pars->Nx << std::endl;
	std::cout << "==========================================" << std::endl;

	return;
}

void read_init(const std::string filename_init_cond, const size_t size, double *U)
{
	std::ifstream input(filename_init_cond.c_str(), std::ios::binary|std::ios::in);
	input.read((char *) U, size * sizeof(double));

	// for (unsigned int i = 0; i < size+1; i++)
	// 	std::cout << i << ":" << U[i] << std::endl;
	return;
}
