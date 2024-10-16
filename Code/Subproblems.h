
#ifndef Subproblems_H
#define Subproblems_H

#include "Common.h"
#include "DataPrep.h"
#include <unordered_map>
#include <list>

struct Subproblems;

typedef pair<IntegerVector, double> Solution;
typedef vector<vector<double>> Matrix;

struct Subproblems
{
	double BestBoundThreshold;
	double CPU = 0;

	IloEnv SPenv;
	IloModel SPmodel;
	IloObjective SPobj;
	IloCplex SPcplex;
	bool Solved = false;

	IloNumVarArray C;
	BoolVarArray2 e;

	ProductPeriods& PP;
	ParameterMap& Parameters;

	IntegerVector InitialSolution;

	int IntegerSolutionLimit = INT_MAX;

	Matrix SP_setup;

	Matrix e_val;

	Subproblems(ProductPeriods& PPIn, ParameterMap& PM);
	~Subproblems();

	void SetupBSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot);
	bool BSP_Solve();
	double GetBSP_CPUTime() { return CPU; }
	double GetBSP_LB();
	double GetBSP_UB();
	Matrix GetBSP_Solutions(int W);
};
#endif
