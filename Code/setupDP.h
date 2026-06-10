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

	BoolVarArray2 e;
	IloNumVarArray u;

	IloEnv TSPenv;
	IloModel TSPmodel;
	IloCplex TSPcplex;
	IloObjective TSPobj;
	
	IloRangeArray OutDeg;          // per-node out-degree rows
	IloRangeArray InDeg;           // per-node in-degree rows
	IloRangeArray Bal;             // per-node flow-balance rows (out - in)
	IloRange TotalEdges;           // global edge-count row (used by the all-open case)
	bool ModelBuilt = false;
	vector<int> PrevActiveNodes;   // nodes activated last call, for cheap reset

public:
	setupDP();
	~setupDP();
	
	int MinSetupMultiPeriod(ProductPeriods& PPs, CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last);
	int MinSetupMiddle(CacheMP& cache, const Matrix& Setups, int First, ProductSet& Middle, int Last);
	int SetupBetween(const Matrix& Setups, int i, int j);
	void generateCombinations(const vector<int>& nums, set<int>& current, int index, vector<set<int>>& all);

	double SetupTSPModel(CacheMP& cache, ProductPeriods& PPs, const Matrix& Setups, int First, ProductSet& Middle, int Last);
	void BuildTSPModel(ProductPeriods& PPs, const Matrix& Setups);

};
#endif

