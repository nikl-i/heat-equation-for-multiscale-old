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

enum processor_t {CPU, GPU};

struct parameters_t
{
	double maxt, tau;
 	double Dx, Dy, Dz;
	double Qx, Qy, Qz;
	int Nx, Ny, Nz;
 	double Lx, Ly, Lz;
};

/* Interfaces */

void Solver1d_cpu(double **U, parameters_t pars);
void check_error(const int info);
void print_matrix(char desc, int m, int n, double* a, int lda);
void read_init(const std::string filename_init_cond, const size_t size, double *U);

void read_args(const int argc, char **argv, std::string *filename_settings,
			   std::string *filename_init_cond,	std::string *filename_output,
			   processor_t *proc, int *dim);
void read_parameters(const std::string filename_settings, parameters_t *pars);

void allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n);
void initMatrix(double *dl, double *d, double *du, const int N,
				const double tau, const double D);
void Output(const std::string filename_output, const int maxstep, const int totalSize, double **U);

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
	int dim = 1;
	processor_t proc = CPU;
	double **U = NULL;

	try
	{
		int maxstep;
		parameters_t pars;
		std::string filename_settings, filename_init_cond, filename_output;

		/* Get parameters */
		read_args(argc, argv, &filename_settings, &filename_init_cond, &filename_output, &proc, &dim);
		if (proc == CPU)
			std::cout << "proc = CPU" << std::endl;
		else if (proc == GPU)
			std::cout << "proc = GPU" << std::endl;
		std::cout << "dim  = " << dim << std::endl;
		read_parameters(filename_settings, &pars);

		return 0;
		maxstep = pars.maxt/pars.tau;

		int totalSize = pars.Nx + 1;
		if (dim > 1)
			totalSize *= pars.Ny + 1;
		if (dim > 2)
		 	totalSize *= pars.Nz + 1;

		U = (double **) malloc((maxstep +1) * sizeof(double *));
		for (int n = 0; n <= maxstep; n++)
			U[n] = (double *) malloc(totalSize * sizeof(double));

		read_init(filename_init_cond, totalSize, U[0]);

		if (proc == CPU)
		{
			if (dim == 1)
				Solver1d_cpu(U, pars);
			else if (dim == 2)
				Solver1d_cpu(U, pars);
			else if (dim == 3)
				Solver1d_cpu(U, pars);
			else
				throw Exception("Invalid dimensions!");
		}
		else if (proc == GPU)
		{
			if (dim == 1)
				Solver1d_cpu(U, pars);
			else if (dim == 2)
				Solver1d_cpu(U, pars);
			else if (dim == 3)
				Solver1d_cpu(U, pars);
			else
				throw Exception("Invalid dimensions!");
		}
		else
			throw Exception("Invalid processor type!");

		Output(filename_output, maxstep, totalSize, U);
		std::cout << "Done." << std::endl;

	}
	catch(Exception error)
	{
		std::cerr << error.msg() << std::endl;
		std::cerr << "Error. Program was terminated." << std::endl;
	}

	/* Deallocation */
	if (U != NULL) free(U);


	return 0;
}


void Solver1d_cpu(double **U, parameters_t pars)
{
	double *B = NULL;
//	double *d = NULL, *du = NULL, *du2 = NULL, *dl = NULL;
	double *d1 = NULL, *du1 = NULL, *du21 = NULL, *dl1 = NULL;
	// double *d2 = NULL, *du2 = NULL, *du22 = NULL, *dl2 = NULL;
	// double *d3 = NULL, *du3 = NULL, *du23 = NULL, *dl3 = NULL;
    // double h, coef, coef_a, coef_b, coef_c;
	double h1, coef1, coef_a1, coef_b1, coef_c1;
	// double h2, coef2, coef_a2, coef_b2, coef_c2;
	// double h3, coef3, coef_a3, coef_b3, coef_c3;

	lapack_int ldb, nrhs; // for Lapack
	lapack_int *ipiv1 = NULL; // for Lapack
	// lapack_int *ipiv2 = NULL; // for Lapack
	// lapack_int *ipiv3 = NULL; // for Lapack
	//allocateData(U, B, totalSize, maxstep);

	int totalSize = pars.Nx + 1;
		// if (dim > 1)
		// 	totalSize *= pars.Ny + 1;
		// if (dim > 2)
		//  	totalSize *= pars.Nz + 1;

	int maxstep = pars.maxt/pars.tau;

	B = (double *) malloc(totalSize * sizeof(double));
	allocateMatrix(&dl1, &d1, &du1, &du21, pars.Nx);
	initMatrix(dl1, d1, du1, pars.Nx, pars.tau, pars.Dx);


	ipiv1 = (lapack_int *) malloc(pars.Nx * sizeof(lapack_int));
 	// if (dim > 1)
	// 	ipiv1 = (lapack_int *) malloc(pars.Ny * sizeof(lapack_int));
	// if (dim > 2)
	// 	ipiv2 = (lapack_int *) malloc(pars.Nz * sizeof(lapack_int));

	nrhs = 1;
	ldb = nrhs;
	check_error(LAPACKE_dgttrf( (lapack_int) pars.Nx+1, dl1, d1, du1, du21, ipiv1));

	/*************************************************************/
	/******************** Main loop ******************************/
	/*************************************************************/
	for (int n = 0; n < maxstep; n++)
	{
		B[0] = 0.0;
		double *Un = U[n]; // current time step
		for (int i = 1; i < pars.Nx; i++)
			B[i] = - coef_a1*Un[i-1] - coef_c1*Un[i+1] + (1.0 + coef_a1 + coef_c1)*Un[i];
		B[pars.Nx] = 0.0;

		check_error(LAPACKE_dgttrs(LAPACK_ROW_MAJOR, 'N', (lapack_int) pars.Nx+1,
								   nrhs, dl1, d1, du1, du21, ipiv1, B, ldb));

		for (int i = 0; i < pars.Nx+1; i++)
			U[n+1][i] = B[i];

		if (n % ITERATION_TO_OUTPUT == 0)
			std::cout << "Iteration t=" << n*pars.tau << " done." << std::endl;
	}
	/*************************************************************/
	/******************** End of main loop ***********************/
	/*************************************************************/

    // allocateMatrix(&dl, &d, &du, &du2, pars.Nx);
	// allocateMatrix(&dl1, &d1, &du1, &du21, pars.Nx);
	// if (dim > 1)
	//  	allocateMatrix(&dl2, &d2, &du2, &du22, pars.Ny);
	// if (dim > 2)
	//  	allocateMatrix(&dl3, &d3, &du3, &du23, pars.Nz);


	/* initialize1d(filename_init_cond, pars, dl, d, du, du2, U) */

	/****init*/
	// initMatrix(dl1, d1, du1, pars.Nx, pars.tau, pars.Dx);

	// initMatrix(dl1, d1, du1, pars.Nx, pars.tau, pars.Dx);
	// if (dim > 1)
  	// 	initMatrix(dl2, d2, du2, pars.Ny, pars.tau, pars.Dy);
	// if (dim > 2)
	// 	initMatrix(dl3, d3, du3, pars.Nz, pars.tau, pars.Dz);



//	if (ipiv != NULL) free(ipiv);*/

	// if (B != NULL) free(B);

	// if (dl1 != NULL) free(dl1);
	// if (d1 != NULL) free(d1);
	// if (du1 != NULL) free(du1);
	// if (du21 != NULL) free(du21);

	// if (dim > 1)
	// {
	// 	if (dl2 != NULL) free(dl2);
	// 	if (d2 != NULL) free(d2);
	// 	if (du2 != NULL) free(du2);
	// 	if (du22 != NULL) free(du22);
	// }

	// if (dim > 2)
	// {
	// 	if (dl3 != NULL) free(dl3);
	// 	if (d3 != NULL) free(d3);
	// 	if (du3 != NULL) free(du3);
	// 	if (du23 != NULL) free(du23);
	// }



	return;
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
			   std::string *filename_init_cond,	std::string *filename_output,
			   processor_t *proc, int *dim)
{
	for (int i = 1; i < argc; i += 2)
	{
		if (!strcmp(argv[i], "-s"))
			*filename_settings = argv[i+1];
		else if (!strcmp(argv[i], "-i"))
			*filename_init_cond = argv[i+1];
		else if (!strcmp(argv[i], "-o"))
			*filename_output = argv[i+1];
		else if (!strcmp(argv[i], "-d"))
			*dim = (int) strtol (argv[i+1], NULL, 10);
		else if (!strcmp(argv[i], "-p"))
		{
			char *p = argv[i+1];
			for ( ; *p; ++p) *p = tolower(*p);
			if (!strcmp(argv[i+1], "cpu"))
				*proc = CPU;
			else if (!strcmp(argv[i+1], "gpu"))
				*proc = GPU;
		}
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


void allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n)
{
	*dl = (double *) malloc(n * sizeof(double));
	*d = (double *) malloc((n+1) * sizeof(double));
	*du = (double *) malloc(n * sizeof(double));
	*du2 = (double *) malloc((n-1) * sizeof(double));

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

void Output(const std::string filename_output, const int maxstep, const int totalSize, double **U)
{
		for (int n = 0; n <= maxstep; n++)
		{
			std::ofstream fileOut (filename_output.c_str(),
								   std::ios::binary | std::ios::out | std::ios::trunc);
			fileOut.write((char*) U[n], totalSize*sizeof(double));
		}
}
