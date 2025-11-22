
#include "Subproblems.h"
#include "DataPrep.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;

Subproblems::Subproblems(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), SPmodel(SPenv), SPcplex(SPenv), CPmodel(CPenv), CPcplex(CPenv)
{
	BestBoundThreshold = DBL_MAX;
}

Subproblems::~Subproblems()
{
	SPenv.end();
	CPenv.end();
}

void Subproblems::SetupBSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot)
{
	// Setup Integer Subroblem (IP)
	SPcplex.setOut(SPenv.getNullStream());
	//Find due-dates and release dates
	vector<int> dd;
	dd.resize(PP.T);
	vector<int> rd;
	rd.resize(PP.T);
	int d_sum = 0;
	for (int t = 0; t < PP.T; ++t) {
		rd[t] = d_sum;
		d_sum += PP.K[t];
		dd[t] = d_sum;
	}

	//Calculate the Big M value
	int TotalCapacity = 0;
	for (int t = 0; t < PP.T; ++t)
		TotalCapacity += PP.K[t];

	int M = TotalCapacity;

	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_a[j] = PP.Products[wp[j]].a;
	SP_a[0] = 0;
	SP_a[W + 1] = 0;

	vector<int> SP_b;
	SP_b.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_b[j] = dd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = d_sum;

	vector<int> SP_r;
	SP_r.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_r[j] = rd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = 0;

	SP_setup.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		SP_setup[j].resize(W + 1);

	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SP_setup[j][l] = PP.Products[wp[j]].s_pr[wp[l]];
		}

	C = CreateNumVarArray(SPenv, W + 2, "C", 0, IloInfinity); //The completion time of job j
	e = CreateBoolVarArray2(SPenv, W + 2, W + 2, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise

	IloExpr SPobj1(SPenv);
	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (j != l)
				SPobj1 += SP_setup[j][l] * e[j][l];

	SPobj = IloMinimize(SPenv, SPobj1);
	SPmodel.add(SPobj);

	SPcplex = IloCplex(SPmodel);

	// Constraint (2)
	for (int j = 0; j < W + 1; ++j) {
		IloExpr Each1(SPenv);
		for (int l = 1; l < W + 2; ++l)
			if (j != l)
				Each1 += e[j][l];
		SPmodel.add(Each1 == 1);
	}

	// Constraint (3)
	for (int j = 1; j < W + 2; ++j) {
		IloExpr Each2(SPenv);
		for (int l = 0; l < W + 1; ++l)
			if (j != l)
				Each2 += e[l][j];
		SPmodel.add(Each2 == 1);
	}

	// Constraint (4)
	for (int j = 0; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SPmodel.add(C[j] + SP_a[l] <= C[l] + M * (1 - e[j][l]));
		}

	// Constraint (5)
	for (int j = 0; j < W + 1; ++j) {
		SPmodel.add(C[j] <= SP_b[j]);
	}
	
	// Constraint (6)
	for (int j = 0; j < W + 1; ++j) {
		SPmodel.add(C[j] >= SP_a[j] + SP_r[j]);
	}

}

bool Subproblems::BSP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		SPcplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		SPcplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		SPcplex.exportModel("IP.lp");
		Solved = SPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetBSP_UB()
{
	return Solved ? SPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetBSP_LB()
{
	return Solved ? SPcplex.getBestObjValue() : DBL_MAX;
}

Matrix Subproblems::GetBSP_Solutions(int W)
{
	e_val.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		e_val[j].resize(W + 1);

	for (int j = 0; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (SPcplex.getValue(e[j][l]) > 0.5)
				e_val[j][l] = SPcplex.getValue(e[j][l]);

	return e_val;
}

void Subproblems::SetupCPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot, Matrix x_val)
{
	// Setup Constraint Programming Subroblem (CP)

	//Find due-dates and release dates
	vector<int> dd;
	dd.resize(PP.T);
	vector<int> rd;
	rd.resize(PP.T);
	int d_sum = 0;
	for (int t = 0; t < PP.T; ++t){
		rd[t] = d_sum;
		for (int p = 0; p < PP.P; ++p) {
			d_sum += x_val[p][t];
		}
		dd[t] = d_sum;
	}

	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W);
	for (int j = 0; j < W; ++j)
		SP_a[j] = PP.Products[wp[j+1]].a;

	vector<int> SP_b;
	SP_b.resize(W);
	for (int j = 0; j < W; ++j)
		SP_b[j] = dd[wt[j+1]];

	vector<int> SP_r;
	SP_r.resize(W);
	for (int j = 0; j < W; ++j)
		SP_r[j] = rd[wt[j+1]];

	SP_setup.resize(W);
	for (int j = 0; j < W; ++j)
		SP_setup[j].resize(W);

	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			SP_setup[j][l] = PP.Products[wp[j+1]].s_pr[wp[l+1]];
		}

	IntArray2 iloSetup = CreateIntArray2(CPenv, W, W);
	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			iloSetup[j][l] = SP_setup[j][l];
		}

	IloIntArray type(CPenv, W);
	for (int j = 0; j < W; ++j)
		type[j] = j;

	IloIntervalVarArray Y(CPenv, W);
	for (int j = 0; j < W; j++)
	{
		Y[j] = IloIntervalVar(CPenv, SP_a[j], "Y");
		Y[j].setPresent();
	}

	IloIntervalSequenceVar V (CPenv, Y, type, "V");

	IloExpr CPobj1(CPenv);
	for (int j = 0; j < W; j++){
		IloIntExpr l = IloTypeOfNext(V, Y[j], j);
		CPobj1 += IloElement(iloSetup[j], l);
	}

	CPobj = IloMinimize(CPenv, CPobj1);
	CPmodel.add(CPobj);

	CPcplex = IloCP(CPmodel);

	CPmodel.add(IloNoOverlap(CPenv, V));

	for (int j = 0; j < W; j++)
		CPmodel.add(IloStartOf(Y[j], 0) >= SP_r[j]);

	for (int j = 0; j < W; j++)
		CPmodel.add(IloEndOf(Y[j], 0) <= SP_b[j]);
}

bool Subproblems::CP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		CPcplex.setParameter(IloCP::Workers, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		CPcplex.setParameter(IloCP::TimeLimit, timeLimit);

	try
	{
		CPcplex.exportModel("CP_GLSP.cpo");
		Solved = CPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetCP_Obj()
{
	return Solved ? CPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetCP_Bound()
{
	return Solved ? CPcplex.getObjBound() : 0;
}

double Subproblems::GetCP_Gap()
{
	return Solved ? CPcplex.getObjGap() : 100;
}