class Solution
{
	int maxm, totalSize;
public:
	double **U;
	Solution(Problem problem);
	~Solution();
	void save(const std::string outputFile);
	void loadInitCond(const std::string inputFile);
};
