#ifndef SETUP_DP_H
#define SETUP_DP_H

#include "Common.h"
#include "DataPrep.h"

typedef vector<vector<double>> Matrix;

typedef set<int> ProductSet;
typedef tuple<int, ProductSet, int> CacheKey;
typedef map<CacheKey, int> CacheMP;

class setupDP
{
	double CPU = 0;

	bool Solved = false;

	IloNumVarArray C;
	BoolVarArray2 e;
	BoolVarArray2 z;
	NumVarArray2 eta;

	Matrix e_val;

	IloEnv TSPenv;

public:
	setupDP();
	~setupDP();
	
	int MinSetupMultiPeriod(ProductPeriods& PPs, CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last);
	int MinSetupMiddle(CacheMP& cache, const Matrix& Setups, int First, ProductSet& Middle, int Last);
	int SetupBetween(const Matrix& Setups, int i, int j);
	void generateCombinations(const vector<int>& nums, set<int>& current, int index, vector<set<int>>& all);

	double SetupTSPModel(ProductPeriods& PPs, const Matrix& Setups, int First, ProductSet& Middle, int Last);
	//double TSP_Solve(double timeLimit);

};
#endif

