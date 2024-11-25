
#ifndef TwoPhase_H
#define TwoPhase_H

#include "Common.h"
#include "DataPrep.h"
#include <unordered_map>
#include <list>

struct TwoPhase;

typedef pair<IntegerVector, double> Solution;

struct TwoPhase
{
	double BestBoundThreshold;
	double CPU = 0;

	IloEnv env;
	IloModel model;
	IloObjective obj;
	IloCplex cplex;
	bool Solved = false;

	NumVarArray2 I; //Inventory level of product p at the end of macro-period t
	NumVarArray2 q; //Amount of product p produced in micro-period s

	BoolVarArray2 x;
	IloNumVar theta;
	IloNumVar gamma;
	IloNumVarArray theta_t;

	ProductPeriods& PP;
	ParameterMap& Parameters;

	IntegerVector InitialSolution;

	int IntegerSolutionLimit = INT_MAX;

	int LB_theta;

	TwoPhase(ProductPeriods& PPIn, ParameterMap& PM);
	~TwoPhase();

	void SetupModel();
	void SetIntegerSolutionLimit(int i) { IntegerSolutionLimit = i; }
	bool Solve(double timeLimit, double* TwoPhase_Iter, double* TwoPhase_Cut, double* TwoPhase_CPU, double* TwoPhase_UB, double* TwoPhase_LB, double* SP_Cons_CPU, double* SP_Solve_CPU, double* MP_CPU, int SPtype);
	double GetCPUTime() { return CPU; }
	double GetLB();
	double GetUB();
};
#endif
