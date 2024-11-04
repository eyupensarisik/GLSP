#include "GLSP_gen.h"
#include "DataPrep.h"
#include <chrono>

GLSP::GLSP(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), model(env), cplex(env)
{
}

GLSP::~GLSP()
{
	env.end();
}


void GLSP::SetupModel_S(double timeLimit)
{
	int P = PP.P;
	int T = PP.T;
	int S = PP.S;
	int M = PP.M;

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity);
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity);

	y = CreateBoolVarArray2(env, P, S, "y");
	z = CreateBoolVarArray3(env, P, P, S, "z");

	// Constraint (2)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal += q[p][s];
			model.add(I[p][t] == I[p][t - 1] + InvBal - Demand);
		}

	// Constraint (2.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal2(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal2 += q[p][s];
			model.add(I[p][t] == InvBal2 - Demand);
		}

	// Constraint (3)
	for (int t = 0; t < T; ++t) {
		IloExpr Cap(env);
		for (int p = 0; p < P; ++p)
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				Cap += PP.Products[p].a * q[p][s];
		model.add(Cap <= PP.K[t]);
	}

	// Constraint (4)
	for (int t = 0; t < T; ++t)
		for (int p = 0; p < P; ++p)
			for (int s = 0; s < S; ++s)
				model.add(q[p][s] * PP.Products[p].a <= PP.K[t] * y[p][s]);

	// Constraint (5)
	for (int s = 0; s < S; ++s) {
		IloExpr Each(env);
		for (int p = 0; p < P; ++p)
			Each += y[p][s];
		model.add(Each == 1);
	}

	// Constraint (6.1)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (s == 0)
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
			else
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1])); 

	// Constraint (6.2)
	/*for (int p = 0; p < P; ++p)
		for (int t = 0; t < T-1; ++t)
			for (int s = 0; s < S; ++s)
				if (s == 0)
					model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
				else if (s == PP.L[t])
					model.add(q[p][s] + q[p][s + 1] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));
				else
					model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));*/


	// Constraint (7)
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 1; s < S; ++s)
				model.add(z[p][r][s] >= y[p][s-1] + y[r][s] - 1);

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			model.add(z[p][r][0] == 0);
	
	// Objective function terms
	// Setup cost
	IloExpr obj1(env);
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			obj1 += (PP.Products[p].f * y[p][s] + PP.Products[p].c * q[p][s]);

	// Production cost
	IloExpr obj2(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				obj2 += PP.Products[p].s_pr[r] * z[p][r][s];

	// Holding cost
	IloExpr obj3(env);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			obj3 += PP.Products[p].h * I[p][t];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	cplex = IloCplex(model);
	
}

void GLSP::SetupModel_NF(double timeLimit)
{
	int P = PP.P;
	int T = PP.T;
	int S = PP.S;
	int M = PP.M;

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity);
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity);

	y = CreateBoolVarArray2(env, P, S, "y");
	z = CreateBoolVarArray3(env, P, P, S, "z");

	// Constraint (2)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal += q[p][s];
			model.add(I[p][t] == I[p][t - 1] + InvBal - Demand);
		}

	// Constraint (2.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal2(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal2 += q[p][s];
			model.add(I[p][t] == InvBal2 - Demand);
		}

	// Constraint (3)
	for (int t = 0; t < T; ++t) {
		IloExpr Cap(env);
		for (int p = 0; p < P; ++p)
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				Cap += PP.Products[p].a * q[p][s];
		model.add(Cap <= PP.K[t]);
	}

	// Constraint (4)
	for (int t = 0; t < T; ++t)
		for (int p = 0; p < P; ++p)
			for (int s = 0; s < S; ++s)
				model.add(q[p][s] * PP.Products[p].a <= PP.K[t] * y[p][s]);

	// Constraint (5)
	for (int s = 0; s < S; ++s) {
		IloExpr Each(env);
		for (int p = 0; p < P; ++p)
			Each += y[p][s];
		model.add(Each == 1);
	}

	// Constraint (6.1)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (s == 0)
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
			else
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));

	// Constraint (6.2)
	/*for (int p = 0; p < P; ++p)
		for (int t = 0; t < T-1; ++t)
			for (int s = 0; s < S; ++s)
				if (s == 0)
					model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
				else if (s == PP.L[t])
					model.add(q[p][s] + q[p][s + 1] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));
				else
					model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));*/
	
	// Constraint (7.1)
	for (int p = 0; p < P; ++p)
		for (int s = 1; s < S; ++s) {
			IloExpr NF1(env);
			for (int r = 0; r < P; ++r)
				NF1 += z[p][r][s];
			model.add(NF1 == y[p][s - 1]);
		}

	// Constraint (7.2)
	for (int r = 0; r < P; ++r)
		for (int s = 0; s < S; ++s) {
			IloExpr NF2(env);
			for (int p = 0; p < P; ++p)
				NF2 += z[p][r][s];
			model.add(NF2 == y[r][s]);
		}
		
		// Objective function terms
		// Setup cost
	IloExpr obj1(env);
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			obj1 += (PP.Products[p].f * y[p][s] + PP.Products[p].c * q[p][s]);

	// Production cost
	IloExpr obj2(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				obj2 += PP.Products[p].s_pr[r] * z[p][r][s];

	// Holding cost
	IloExpr obj3(env);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			obj3 += PP.Products[p].h * I[p][t];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	cplex = IloCplex(model);

}

bool GLSP::Solve(double timeLimit)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		cplex.exportModel("GLSP.lp");
		Solved = cplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return Solved;
}

double GLSP::GetUB()
{
	return Solved ? cplex.getObjValue() : DBL_MAX;
}

double GLSP::GetLB()
{
	return Solved ? cplex.getBestObjValue() : DBL_MAX;
}

double GLSP::GetGap()
{
	return Solved ? cplex.getMIPRelativeGap() : DBL_MAX;
}

int GLSP::GetConsts()
{
	return Solved ? cplex.getNrows() : 0;
}

int GLSP::GetVars()
{
	return Solved ? cplex.getNcols() : 0;
}

void GLSP::GetSolutions(int P, int T, int S)
{
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			if (cplex.getValue(I[p][t]) > 0)
				cout << "Value of I[" << p << "][" << t << "] is " << cplex.getValue(I[p][t]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(q[p][s]) > 0)
				cout << "Value of q[" << p << "][" << s << "] is " << cplex.getValue(q[p][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(y[p][s]) > 0)
				cout << "Value of y[" << p << "][" << s << "] is " << cplex.getValue(y[p][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				if (cplex.getValue(z[p][r][s]) > 0)
					cout << "Value of z[" << p << "][" << r << "][" << s << "] is " << cplex.getValue(z[p][r][s]) << endl;
}