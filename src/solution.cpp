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
#include "solution.h"

using namespace std;

void Solution::save(const string outputFile)
{
	cout << "Saving solution... ";
	ofstream fileOut (outputFile.c_str(),
	 				  ios::binary | ios::out | ios::trunc);
	for (int n = 0; n <= maxm; n++)
	  	fileOut.write((char*) U[n],totalSize*sizeof(double));
	cout << "Done." << endl;
}

Solution::Solution(Problem problem)
{
	maxm = problem.get_maxm();
	totalSize = problem.get_totalSize();
	U = new double * [maxm+1];
	for (int n = 0; n <= maxm; n++)
	{
		U[n] = new double [totalSize];
		memset(U[n], 0, totalSize*sizeof(double));
	}
	loadInitCond(problem.get_initCondFile());
}

Solution::~Solution()
{
	for (int n = 0; n <= maxm; n++)
		if (U[n] != NULL) delete [] U[n];
	if (U != NULL) delete [] U;
	return;
}

void Solution::loadInitCond(const string filename)
{
	cout << "Loading initial condition... ";
	ifstream input(filename.c_str(), ios::binary|ios::in);
	input.read((char *) U[0], totalSize * sizeof(double));
	cout << "Done." << endl;
	return;
}

// void readArgs(const int argc, char **argv, string *settingsFile,
// 			   string *initCondFile,	string *outputFile,
// 			   processor_t *proc);

// void readParameters(const string settingsFile, parameters_t *pars);

// /** Get boundary conditions */
// if (problem.boundCondTypeX0 > 0)
// {
// 	q = (double **) malloc(2*problem.dim  * sizeof(double*));
// 	qx0 = (double *) malloc((maxstep + 1) * sizeof(double));
// 	//readFile(filename_bcqX0, maxstep + 1, qx0);
// 	for (int n = 0; n < maxstep +1; n++)
// 		qx0[n] = 0.0;

// }
// if (problem.boundCondTypeX1 > 0)
// {
// 	if (q == NULL) q = (double **) malloc(2*problem.dim  * sizeof(double*));
// 	qx1 = (double *) malloc((maxstep + 1) * sizeof(double));
// 	for (int n = 0; n < maxstep +1; n++)
// 		qx1[n] = 0.0;
// 	//readFile(filename_bcqX1, maxstep + 1, qx0);
// }

// if (problem.boundCondTypeX0 == 3)
// {
// 	p = (double **) malloc(2*problem.dim  * sizeof(double*));
// 	px0 = (double *) malloc((maxstep + 1) * sizeof(double));
// 	//readFile(filename_bcpX0, maxstep + 1, px0);
// }
// if (problem.boundCondTypeX1 == 3)
// {
// 	if (p == NULL) p = (double **) malloc(2*problem.dim  * sizeof(double*));
// 	px1 = (double *) malloc((maxstep + 1) * sizeof(double));
// 	//readFile(filename_bcpX1, maxstep + 1, px0);
// }

/* Get parameters */
// readArgs(argc, argv, &settingsFile, &initCondFile, &outputFile, &proc);
// if (proc == CPU)
// 	cout << "proc = CPU" << endl;
// else if (proc == GPU)
// 	cout << "proc = GPU" << endl;
// readParameters(settingsFile, &pars);

// void readArgs(const int argc, char **argv, string *settingsFile,
// 			   string *initCondFile,	string *outputFile,
// 			   processor_t *proc)
// {
// 	for (int i = 1; i < argc; i += 2)
// 	{
// 		if (!strcmp(argv[i], "-s"))
// 			*settingsFile = argv[i+1];
// 		else if (!strcmp(argv[i], "-i"))
// 			*initCondFile = argv[i+1];
// 		else if (!strcmp(argv[i], "-o"))
// 			*outputFile = argv[i+1];
// 		else if (!strcmp(argv[i], "-p"))
// 		{
// 			char *p = argv[i+1];
// 			for ( ; *p; ++p) *p = tolower(*p);
// 			if (!strcmp(argv[i+1], "cpu"))
// 				*proc = CPU;
// 			else if (!strcmp(argv[i+1], "gpu"))
// 				*proc = GPU;
// 		}
// 		else
// 			throw Exception("Unknown option!");
// 	}

// 	if ((*settingsFile).empty())
// 		*settingsFile = "settings";
// 	if ((*initCondFile).empty())
// 		*initCondFile = "u0.bin";
// 	if ((*outputFile).empty())
// 		*outputFile = "u.bin";

// 	exists(*settingsFile);
// 	exists(*initCondFile);


// 	return;
// }

// void readParameters(const string settingsFile, parameters_t *pars)
// {
// 	string line;
// 	ifstream file(settingsFile.c_str());

// 	while(getline(file, line))
// 	{
// 		string token;
// 		string keyname;
// 		istringstream tokens(line) ;
// 		while(tokens >> token)
// 		{
// 			if (token == "#") break;
// 			if (token == "=") continue;

// 			if (keyname.empty())
// 			{
// 				transform(token.begin(), token.end(), token.begin(), ::tolower);
// 				keyname = token;
// 			}
// 			else
// 			{
// 				double value = atof(token.c_str());

// 				if (keyname == "dim")
// 					pars->dim = (int) value;
// 				else if (keyname == "tau")
// 					pars->tau = value;
// 				else if (keyname == "maxt")
// 					pars->maxt = value;
// 				else if (keyname == "nx")
// 					pars->Nx = (int) value;
// 				else if (keyname ==  "dx")
// 					pars->Dx = value;
// 				else if (keyname == "qx")
// 					pars->Qx = value;
// 				else if (keyname == "ny")
// 					pars->Ny = (int) value;
// 				else if (keyname ==  "dy")
// 					pars->Dy = value;
// 				else if (keyname == "qy")
// 					pars->Qy = value;
// 				else if (keyname == "nz")
// 					pars->Nz = (int) value;
// 				else if (keyname ==  "dz")
// 					pars->Dz = value;
// 				else if (keyname == "qz")
// 					pars->Qz = value;
// 				else if (keyname == "iter2skip")
// 					pars->numberOfIterationsToSkip = value;
// 				else if (keyname == "bcx0")
// 					pars->boundCondTypeX0 = value;
// 				else if (keyname == "bcx1")
// 					pars->boundCondTypeX1 = value;
// 				else
// 					throw(Exception("Unknown parameter name: " + keyname + " = " + token));

// 				keyname.clear();
// 			}
// 		}
// 	}

// 	cout << "================================================" << endl;
// 	return;
// }
