enum Direction {X, Y, Z};

class Solver
{
	double hx, coefx, coef_ax, coef_cx, coef_bx;
	double hy, coefy, coef_ay, coef_by, coef_cy;
	double hz, coefz, coef_az, coef_bz, coef_cz;
	int nrhs, ldb;

	int *ipiv1 = NULL, *ipiv2 = NULL, *ipiv3 = NULL; // for Lapack
	double *B = NULL, *Un = NULL, *temp = NULL;
	double *d1 = NULL, *du1 = NULL, *du21 = NULL, *dl1 = NULL;
	double *d2 = NULL, *du2 = NULL, *du22 = NULL, *dl2 = NULL;
	double *d3 = NULL, *du3 = NULL, *du23 = NULL, *dl3 = NULL;

	double get_alphaX(Problem1d problem, int n = 0);
	double get_alphaY(Problem2d problem, int n = 0);
	double get_alphaZ(Problem3d problem, int n = 0);
	double get_betaX(Problem1d problem);
	double get_betaY(Problem2d problem);
	double get_betaZ(Problem3d problem);
	double get_gammaX(Problem1d problem);
	double get_gammaY(Problem2d problem);
	double get_gammaZ(Problem3d problem);
	double get_deltaX(Problem1d problem, int n = 0);
	double get_deltaY(Problem2d problem, int n = 0);
	double get_deltaZ(Problem3d problem, int n = 0);
	double get_b0X(Problem1d problem, double *Un, int n = 0);
	double get_bNX(Problem1d problem, double *Un, int n = 0);
	double get_b0Y(Problem2d problem, double *Un, int n = 0);
	double get_bNY(Problem2d problem, double *Un, int n = 0);
	double get_b0Z(Problem3d problem, double *Un, int n = 0);
	double get_bNZ(Problem3d problem, double *Un, int n = 0);

	void allocateMatrix(double **dl, double **d, double **du, double **du2,
					const int n);
	void initMatrix(double *dl, double *d, double *du, const int N,
				const double coef_a, const double coef_b, const double coef_c,
				const double coef_b0, const double coef_c0,
				const double coef_aN, const double coef_bN);
	void copyTransposed3d(double *B, const double *A, const int N1,
						  const int N2, const int N3);
	void copyTransposed(double *B, const double *A, const int N1, const int N2);
	void switchPointers(double **A, double **B);

public:
	Solver(Problem1d problem);
	Solver(Problem2d problem);
	Solver(Problem3d problem);
	~Solver();
	void solve(Problem1d problem, Solution *U);
	void solve(Problem2d problem, Solution *U);
	void solve(Problem3d problem, Solution *U);
};
