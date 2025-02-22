#ifndef GLSP_H
#define GLSP_H

#include "Common.h"
#include "DataPrep.h"
#include "setupDP.h"

class GLSP
{
	ProductPeriods& PP;
	ParameterMap& Parameters;

	double CPU = 0;

	IloEnv env;
	IloModel model;
	IloObjective obj;
	IloCplex cplex;
	bool Solved = false;

	NumVarArray2 I; //Inventory level of product p at the end of macro-period t
	NumVarArray2 q; //Amount of product p produced in micro-period s

	BoolVarArray2 y; //1, if the machine is prepared for the product p in micro-period s; 0, otherwise
	NumVarArray3 z; //1 if the transition from product p to r occurs at the beginning of micro - period s; 0 otherwise

	BoolVarArray3 zz; //1 if the transition from product p to r occurs at the beginning of micro - period s; 0 otherwise
	NumVarArray2 yy;
	NumVarArray2 rr;
	BoolVarArray2 w;

	setupDP seq_DP;
	CacheMP SOcache;
	
public:
	GLSP(ProductPeriods& PPIn, ParameterMap& PM);
	~GLSP();

	void SetupModel_S(double timeLimit);
	void SetupModel_SO(double timeLimit);
	void SetupModel_NF(double timeLimit);
	void SetupModel_CC(double timeLimit);
	void SetupModel_TF(double timeLimit);
	bool Solve(double timeLimit);
	double GetCPUTime() { return CPU; }
	double GetLB();
	double GetUB();
	double GetGap();
	int GetConsts();
	int GetVars();
	void GetSolutions(int SolP, int SolT, int S);
	void GetSolutions_CC(int SolP, int SolT, int S);
	void GetSolutions_TF(int SolP, int SolT, int S);
	void GetSolutions_SO(int SolP, int SolT, int S);
};
#endif
