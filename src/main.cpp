/*****************************************************************
\todo Add __FUNCNAME__ to Exception
\todo Time measurement
\todo Check w\o OpenBLAS
\todo Comment the code
****************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>

#include "errhandle.h"
#include "problem.h"
#include "solution.h"
#include "solver.h"
#include "settings.h"

using namespace std;
int main(int argc, char **argv)
{
#ifdef _OPENMP
	cout << "OpenMP is enabled" << endl;
#endif
	try
	{
		SETTINGS_PROBLEM_TYPE problem;
		Solution solution(problem);
		Solver solver(problem);
		problem.loadElectricIntensity();
		solver.solve(problem, &solution);
		solution.save(problem.get_outputFile());
		problem.clearElectricIntensity();

	}
	catch(Exception error)
	{
		cerr << error.msg() << endl;
		cerr << "Error. Program was terminated." << endl;
	}
	return 0;
}
