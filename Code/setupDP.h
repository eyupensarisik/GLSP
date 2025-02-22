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

	

public:
	setupDP();
	
	int MinSetupMultiPeriod(CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last);
	int MinSetupMiddle(CacheMP& cache, const Matrix& Setups, int First, ProductSet& Middle, int Last);
	int SetupBetween(const Matrix& Setups, int i, int j);
	
};
#endif

