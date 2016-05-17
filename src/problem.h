class Problem
{
	friend class Solver;
	std::string initCondFile, outputFile;
	int numberOfIterationsToSkip;
	double maxt, tau;
	int maxstep, maxm;
	double heatSrc(const double *U, const int n, const int i = 0,
				   const int j = 0, const int k = 0);
protected:
	int Nx, Ny, Nz;
	int totalSize;
	double *E = NULL;
	void loadElectricIntensity(const std::string filename);
public:
	Problem();
	void clearElectricIntensity();
	void print() const;
	int get_maxm() const {return maxm;};
	int get_totalSize() const {return totalSize;};
	std::string get_initCondFile() const {return initCondFile;};
	std::string get_outputFile() const {return outputFile;};
};

class Problem1d : public Problem
{
	friend class Solver;
protected:
	double Dx;
	int boundCondTypeX0, boundCondTypeX1;
	double qx0(const int n, const int j = 0, const int k = 0);
	double px0(const int n, const int j = 0, const int k = 0);
	double qx1(const int n, const int j = 0, const int k = 0);
	double px1(const int n, const int j = 0, const int k = 0);
public:
	Problem1d();
	void print() const;
};

class Problem2d : public Problem1d
{
	friend class Solver;
protected:
	double Dy;
	int boundCondTypeY0, boundCondTypeY1;
	double qy0(const int n, const int i = 0, const int k = 0);
	double py0(const int n, const int i = 0, const int k = 0);
	double qy1(const int n, const int i = 0, const int k = 0);
	double py1(const int n, const int i = 0, const int k = 0);
public:
	Problem2d();
	void print() const;
};

class Problem3d : public Problem2d
{
	friend class Solver;
	double Dz;
	int boundCondTypeZ0, boundCondTypeZ1;
	double qz0(const int n, const int i = 0, const int j = 0);
	double pz0(const int n, const int i = 0, const int j = 0);
	double qz1(const int n, const int i = 0, const int j = 0);
	double pz1(const int n, const int i = 0, const int j = 0);
public:
	Problem3d();
	void print() const;
};
