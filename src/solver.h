class Solver
{
	double h1, coef1, coef_a1, coef_c1, coef_b1;
	double h2, coef2, coef_a2, coef_b2, coef_c2;
	double h3, coef3, coef_a3, coef_b3, coef_c3;
	double coef_d1, coef_d2, coef_d3, coef_d4;
	int nrhs, ldb;

	int *ipiv1 = NULL, *ipiv2 = NULL, *ipiv3 = NULL; // for Lapack
	double *B = NULL, *Un = NULL, *temp = NULL;
	double *d1 = NULL, *du1 = NULL, *du21 = NULL, *dl1 = NULL;
	double *d2 = NULL, *du2 = NULL, *du22 = NULL, *dl2 = NULL;
	double *d3 = NULL, *du3 = NULL, *du23 = NULL, *dl3 = NULL;

	double get_d1(Problem1d problem, int n = 0);
	double get_d2(Problem1d problem);
	double get_d3(Problem1d problem);
	double get_d4(Problem1d problem, int n = 0);
	double get_b0(Problem1d problem, double *Un, int n = 0);
	double get_bN(Problem1d problem, double *Un, int n = 0);

	void allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n);
	void initMatrix(double *dl, double *d, double *du, const int N,
				const double coef_a, const double coef_b, const double coef_c,
				const double coef_b0, const double coef_c0,
				const double coef_aN, const double coef_bN);
	void assignTransposed3d(double **B, double *A,
						const int N1, const int N2, const int N3);
	void assignTransposed(double **B, double *A, const int N1, const int N2);

public:
	Solver(Problem1d problem);
	Solver(Problem2d problem);
	Solver(Problem3d problem);
	~Solver();
	void solve(Problem1d problem, Solution *U);
	void solve(Problem2d problem, Solution *U);
	void solve(Problem3d problem, Solution *U);
};
