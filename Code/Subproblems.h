
#ifndef Subproblems_H
#define Subproblems_H

#include "Common.h"
#include "DataPrep.h"
#include "ilcp/cp.h"
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
	BoolVarArray2 z;
	NumVarArray2 eta;

	Matrix SP_setup;
	Matrix e_val;

	IloEnv TSPenv;
	IloModel TSPmodel;
	IloObjective TSPobj;
	IloCplex TSPcplex;

	IloEnv CPenv;
	IloModel CPmodel;
	IloObjective CPobj;
	IloCP CPcplex;

	IloEnv TSP_CPenv;
	IloModel TSP_CPmodel;
	IloObjective TSP_CPobj;
	IloCP TSP_CPcplex;

	IloExprArray nexts;
	IloIntervalSequenceVar V;

	ProductPeriods& PP;
	ParameterMap& Parameters;

	IntegerVector InitialSolution;

	int IntegerSolutionLimit = INT_MAX;

	

	Subproblems(ProductPeriods& PPIn, ParameterMap& PM);
	~Subproblems();

	void SetupBSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot);
	bool BSP_Solve(double timeLimit);
	double GetBSP_CPUTime() { return CPU; }
	double GetBSP_LB();
	double GetBSP_UB();
	Matrix GetBSP_Solutions(int W);

	void SetupTSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot);
	bool TSP_Solve(double timeLimit);
	double GetTSP_CPUTime() { return CPU; }
	double GetTSP_LB();
	double GetTSP_UB();
	Matrix GetTSP_Solutions(int W);

	void SetupCPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot, Matrix x_val);
	bool CP_Solve(double timeLimit);
	double GetCP_CPUTime() { return CPU; }
	double GetCP_Obj();
	double GetCP_Bound();
	double GetCP_Gap();

	void SetupTSP_CPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot, Matrix x_val);
	bool TSP_CP_Solve(double timeLimit);
	double GetTSP_CP_CPUTime() { return CPU; }
	double GetTSP_CP_Obj();
	double GetTSP_CP_Bound();
	double GetTSP_CP_Gap();
};
#endif
