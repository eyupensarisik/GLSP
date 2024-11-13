
#ifndef TwoPhaseCallback_H
#define TwoPhaseCallback_H

#include "Common.h"
#include "DataPrep.h"
#include "setupDP.h"
#include <unordered_map>
#include <list>
#include <mutex>

typedef vector<vector<double>> Matrix;

struct TwoPhaseC;

struct Worker
{
	TwoPhaseC* pTwoPhaseC;
	size_t CallCount;
	size_t CutCount;
	double CPU;
	double TotalSPtime;
	double TotalSPcons_time;

	Worker(TwoPhaseC* pTwoPhaseC);

	virtual bool separate(ProductPeriods& PP, ParameterMap& Parameters, int& LB_theta, const IloNum thetaVal, const NumArray2& xSol, double& OptimalCost,
		IloExpr& cutLhs, IloExpr& cutGMs, IloExpr& Sum1, IloExpr& Sum2, IloExpr& GMCut, IloExprArray& DPCut, IloExprArray& cutDPs, 
		Matrix setup_pr) = 0;
};

struct WorkerWW : public Worker
{
	WorkerWW(TwoPhaseC* pTwoPhaseC);

	bool separate(ProductPeriods& PP, ParameterMap& Parameters,  int& LB_theta, const IloNum thetaVal, const NumArray2& xSol, double& OptimalCost,
		IloExpr& cutLhs, IloExpr& cutGMs, IloExpr& Sum1, IloExpr& Sum2, IloExpr& GMCut, IloExprArray& DPCut, IloExprArray& cutDPs, 
		Matrix setup_pr);
};

struct TwoPhaseCallback : public IloCplex::Callback::Function
{
	TwoPhaseC* pTwoPhaseC;
	vector<Worker*> workers;

	TwoPhaseCallback(TwoPhaseC* pTwoPhaseC);
	~TwoPhaseCallback();

	void invoke(const IloCplex::Callback::Context& context);
};

typedef pair<IntegerVector, double> Solution;

struct TwoPhaseC
{
	mutex Mutex;
	unique_ptr<TwoPhaseCallback> pCallback;
	double BestBoundThreshold;

	const double Epsilon = 0.001;
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
	double W_SPtime = 0;
	double W_SPcons_time = 0;

	setupDP sDP;
	map <pair<int, set<int>>, int> Cache;
	Matrix setup_pr;

	list<IloConstraint> GeneratedCuts;

	TwoPhaseC(ProductPeriods& PPIn, ParameterMap& PM);
	~TwoPhaseC();

	void SetupModel();
	void SetIntegerSolutionLimit(int i) { IntegerSolutionLimit = i; }
	bool Solve(double timeLimit);
	double GetCPUTime() { return CPU; }
	double GetLB();
	double GetUB();
	double GetCallbackCPU();
	size_t GetCallCount();
	size_t GetCutCount();
	double GetSPconsCPU();
	double GetSPsolveCPU();
};
#endif