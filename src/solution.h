class Solution
{
	int maxm, totalSize;
	double **U;
public:
	Solution(const Problem& problem);
	~Solution();
	double **get_data() {return U;};
	double *get_data(const int i) {return U[i];};
	void save(const std::string& outputFile);
	void loadInitCond(const std::string& inputFile);
};
