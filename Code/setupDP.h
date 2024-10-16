#ifndef SETUP_DP_H
#define SETUP_DP_H

#include "Common.h"
#include "DataPrep.h"

typedef vector<vector<double>> Matrix;

class setupDP
{
	double CPU = 0;

	

public:
	setupDP();
	
	int MinSetup(map<pair<int, set<int>>, int>& Cache, const Matrix& Setups, int Prev, set<int>& Current);
	
};
#endif

