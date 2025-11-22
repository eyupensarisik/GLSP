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

void GLSP::SetupModel_S(double timeLimit) //Setup the standard model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods
	int S = PP.S; //The set of micro - periods

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	y = CreateBoolVarArray2(env, P, S, "y"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise
	z = CreateNumVarArray3(env, P, P, S, "z"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise

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

	// Constraint (6)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (s == 0)
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
			else
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1])); 

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
	total_obj = obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	cplex = IloCplex(model);
	
}

void GLSP::SetupModel_NF(double timeLimit) //Setup the network-flow model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods
	int S = PP.S; //The set of micro - periods

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	y = CreateBoolVarArray2(env, P, S, "y"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise
	z = CreateNumVarArray3(env, P, P, S, "z"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise

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

	// Constraint (6)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (s == 0)
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s]));
			else
				model.add(q[p][s] >= PP.Products[p].m * (y[p][s] - y[p][s - 1]));

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
	total_obj = obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	cplex = IloCplex(model);

}

void GLSP::SetupModel_CC(double timeLimit) //Setup the Clark and Clark (2000) model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods
	int S = PP.S; //The set of micro - periods

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	zz = CreateBoolVarArray3(env, P, P, S, "T"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise

	// Objective function terms
	// Setup cost
	IloExpr obj1(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				obj1 += PP.Products[p].s_pr[r] * zz[p][r][s];

	// Holding cost
	IloExpr obj2(env);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			obj2 += PP.Products[p].h * I[p][t];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

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

	// Constraint (11)
	IloExpr zzExpr(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			zzExpr += zz[r][p][0];
	model.add(zzExpr == 1);

	// Constraint (12)
	for (int t = 0; t < T; ++t)
		for (int p = 0; p < P; ++p)
			for (int s = 0; s < S; ++s){
				IloExpr Cap2(env);
				for (int r = 0; r < P; ++r)
					Cap2 += zz[r][p][s];
				model.add(q[p][s] * PP.Products[p].a <= PP.K[t] * Cap2);
			}
	// Constraint (13)
	for (int p = 0; p < P; ++p)
		for (int s = 1; s < S; ++s) {
			IloExpr Each1(env);
			IloExpr Each2(env);
			for (int r = 0; r < P; ++r){
				Each1 += zz[r][p][s-1];
				Each2 += zz[p][r][s];
			}
			model.add(Each1 == Each2);
		}

	// Constraint (14)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s){
			IloExpr MOQ(env);
			for (int r = 0; r < P; ++r)
				if (p != r)
					MOQ += zz[r][p][s];
			model.add(q[p][s] >= PP.Products[p].m * MOQ);
		}
	
	cplex = IloCplex(model);

}

void GLSP::SetupModel_TF(double timeLimit) //Setup the Time Flow model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods
	int S = PP.S; //The set of micro - periods
	S = P * T;

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	yy = CreateNumVarArray2(env, P, S, "y"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise
	rr = CreateNumVarArray2(env, P, S, "r"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise

	zz = CreateBoolVarArray3(env, P, P, S, "z"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise
	ww = CreateNumVarArray3(env, P, P, S, "w");
	
	vector<vector<int>> st;
	st.resize(P);
	for (int p = 0; p < P; ++p)
		st[p].resize(P);

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			st[p][r] = 0;

	vector<vector<vector<int>>> L;
	L.resize(P);
	for (int p = 0; p < P; ++p){
		L[p].resize(P);
		for (int r = 0; r < P; ++r)
			L[p][r].resize(S);
	}

	vector<vector<vector<int>>> U;
	U.resize(P);
	for (int p = 0; p < P; ++p) {
		U[p].resize(P);
		for (int r = 0; r < P; ++r)
			U[p][r].resize(S);
	}

	vector<int> startK; 
	vector<int> endK;
	startK.resize(T); 
	endK.resize(T);

	int KK = 0;
	for (int t = 0; t < T; ++t){
		startK[t] += KK;

		KK += PP.K[t];

		endK[t] += KK;
	}

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int t = 0; t < T; ++t)
				for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
					if (s == PP.S_b[t])
						L[p][r][s] = startK[t] - st[p][r];
					else
						L[p][r][s] = startK[t];

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int t = 0; t < T; ++t)
				for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
					if (s == PP.S_b[t])
						U[p][r][s] = startK[t];
					else
						U[p][r][s] = endK[t] - st[p][r];


	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s){
				model.add(zz[p][r][s]);
				model.add(ww[p][r][s]);
			}

	// Objective function terms
	// Setup cost
	IloExpr obj1(env);
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			obj1 += (PP.Products[p].f * yy[p][s] + PP.Products[p].c * q[p][s]);

	// Production cost
	IloExpr obj2(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				obj2 += PP.Products[p].s_pr[r] * zz[p][r][s];

	// Holding cost
	IloExpr obj3(env);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			obj3 += PP.Products[p].h * I[p][t];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	// Constraint (2b)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal += q[p][s];
			model.add(I[p][t] == I[p][t - 1] + InvBal - Demand);
		}

	// Constraint (2b.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			double Demand = PP.Products[p].d[t];
			IloExpr InvBal2(env);
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				InvBal2 += q[p][s];
			model.add(I[p][t] == InvBal2 - Demand);
		}

	// Constraint (2c)
	for (int t = 0; t < T; ++t) {
		IloExpr expr1(env);
		for (int p = 0; p < P; ++p)
			for (int s = PP.S_b[t]; s < PP.S_f[t] - 1; ++s) {
				IloExpr exprW(env);
				for (int r = 0; r < P; ++r)
					exprW += ww[r][p][s] + st[r][p] * zz[r][p][s];

				IloExpr exprW2(env);
				for (int k = 0; k < P; ++k)
					exprW2 += ww[p][k][s+1];

			model.add(exprW + PP.Products[p].a * q[p][s] + rr[p][s] == exprW2);
			}
	}


	// Constraint (2d)
	for (int t = 0; t < T; ++t)
		for (int p = 0; p < P; ++p)
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				model.add(PP.Products[p].a * q[p][s] + rr[p][s] <= PP.K[t] * yy[p][s]);

	// Constraint (2e)
	IloExpr yyExpr(env);
	for (int p = 0; p < P; ++p)
		yyExpr += yy[p][0];
	model.add(yyExpr == 1);

	// Constraint (2f)
	for (int p = 0; p < P; ++p)
		for (int s = 1; s < S; ++s) {
			IloExpr expr3(env);
			for (int r = 0; r < P; ++r)
				expr3 += zz[r][p][s];
			model.add(yy[p][s] == expr3);
		}

	// Constraint (2g)
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S - 1; ++s) {
			IloExpr expr4(env);
			for (int r = 0; r < P; ++r)
				expr4 += zz[p][r][s + 1];
			model.add(yy[p][s] == expr4);
		}

	// Constraint (2h.1)
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 1; s < S; ++s)
				model.add(L[p][r][s] * zz[p][r][s] <= ww[p][r][s]);

	// Constraint (2h.2)
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 1; s < S; ++s)
				model.add(ww[p][r][s] <= U[p][r][s] * zz[p][r][s]);

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

void GLSP::GetSolutions_CC(int P, int T, int S)
{
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			if (cplex.getValue(I[p][t]) > 0)
				cout << "Value of I[" << p << "][" << t << "] is " << cplex.getValue(I[p][t]) << endl;

	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			if (cplex.getValue(I[p][t]) > 0)
				cout << "Value of I[" << p << "][" << t << "] is " << cplex.getValue(I[p][t]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(q[p][s]) > 0)
				cout << "Value of q[" << p << "][" << s << "] is " << cplex.getValue(q[p][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				if (cplex.getValue(zz[p][r][s]) > 0)
					cout << "Value of T[" << p << "][" << r << "][" << s << "] is " << cplex.getValue(zz[p][r][s]) << endl;
}

void GLSP::GetSolutions_TF(int P, int T, int S)
{
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				if (cplex.getValue(zz[p][r][s]) > 0)
					cout << "Value of z[" << p << "][" << r << "][" << s << "] is " << cplex.getValue(zz[p][r][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(yy[p][s]) > 0)
				cout << "Value of y[" << p << "][" << s << "] is " << cplex.getValue(yy[p][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(q[p][s]) > 0)
				cout << "Value of q[" << p << "][" << s << "] is " << cplex.getValue(q[p][s]) << endl;

	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			if (cplex.getValue(rr[p][s]) > 0)
				cout << "Value of r[" << p << "][" << s << "] is " << cplex.getValue(rr[p][s]) << endl;
}