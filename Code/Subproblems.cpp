
#include "Subproblems.h"
#include "DataPrep.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;

Subproblems::Subproblems(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), SPmodel(SPenv), SPcplex(SPenv)
{
	BestBoundThreshold = DBL_MAX;
}

Subproblems::~Subproblems()
{
	SPenv.end();
}

void Subproblems::SetupBSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot)
{
	// Setup BSP Subroblem
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

	vector<int> SP_h;
	SP_h.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_h[j] = PP.Products[wp[j]].h;
	SP_h[0] = 0;
	SP_h[W + 1] = 0;

	SP_setup.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		SP_setup[j].resize(W + 1);

	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SP_setup[j][l] = PP.Products[wp[j]].s_pr[wp[l]];
		}

	C = CreateNumVarArray(SPenv, W + 2, "C", 0, IloInfinity);
	e = CreateBoolVarArray2(SPenv, W + 2, W + 2, "e");

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

	//SPcplex.exportModel("SP_GLSP.lp");

}

bool Subproblems::BSP_Solve()
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		SPcplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));

	try
	{
		//SPcplex.exportModel("BSP.lp");
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