#ifndef TwoPhase_H
#define TwoPhase_H

#include "Common.h"
#include "DataPrep.h"
#include <unordered_map>
#include <list>
#include <mutex>

struct TwoPhase;

struct Worker
{
	TwoPhase* pTwoPhase;
	size_t CallCount;
	size_t CutCount;
	double CPU;

	Worker(TwoPhase* pTwoPhase);

	virtual bool separate(ProductPeriods& PP, const IloNum thetaVal, const IloNumArray& zSol, double& OptimalCost, IloExpr& cutLhs, bool applyLifting) = 0;
};

struct WorkerWW : public Worker
{
	WorkerWW(TwoPhase* pTwoPhase);

	bool separate(ProductPeriods& PP, const IloNum thetaVal, const IloNumArray& zSol, double& OptimalCost, IloExpr& cutLhs, bool applyLifting);
};

struct TwoPhaseCallback : public IloCplex::Callback::Function
{
	TwoPhase* pTwoPhase;
	vector<Worker*> workers;

	TwoPhaseCallback(TwoPhase* pTwoPhase);
	~TwoPhaseCallback();

	void invoke(const IloCplex::Callback::Context& context);
};

typedef pair<IntegerVector, double> Solution;

struct TwoPhase
{
	mutex Mutex;
	unique_ptr<TwoPhaseCallback> pCallback;
	double BestBoundThreshold;

	double CPU = 0;

	IloEnv env;
	IloModel model;
	IloObjective obj;
	IloCplex cplex;
	bool Solved = false;

	NumVarArray2 I; //Inventory level of product p at the end of macro-period t
	NumVarArray2 q; //Amount of product p produced in micro-period s

	BoolVarArray2 y; //1, if the machine is prepared for the product p in micro-period s; 0, otherwise
	BoolVarArray3 z; //1 if the transition from product p to r occurs at the beginning of micro - period s; 0 otherwise

	ProductPeriods& PP;

	IntegerVector InitialSolution;

	TwoPhase(ProductPeriods& PPIn);
	~TwoPhase();

	void SetupModel();

	bool Solve(double timeLimit);
	double GetCPUTime() { return CPU; }
	double GetLB();
	double GetUB();
	double GetGap();
	void GetSolutions(int SolP, int SolT, int S);
	double GetCallbackCPU();
	size_t GetCallCount();
};
#endif
