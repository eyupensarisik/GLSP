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
	//z = CreateBoolVarArray3(env, P, P, S, "z"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise

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

void GLSP::SetupModel_TF(double timeLimit) //Setup the Clark and Clark (2000) model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods
	int S = PP.S; //The set of micro - periods

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	yy = CreateNumVarArray2(env, P, S, "y"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise
	rr = CreateNumVarArray2(env, P, S, "r"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise

	zz = CreateBoolVarArray3(env, P, P, S, "z"); //1 if the transition from product p to r occurs at the beginning of micro - period s, 0 otherwise

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				model.add(zz[p][r][s]);

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
	for (int t = 0; t < T; ++t){
		IloExpr expr1(env);
		for (int p = 0; p < P; ++p)
			for (int s = PP.S_b[t]; s < PP.S_f[t]; ++s)
				expr1 += PP.Products[p].a * q[p][s] + rr[p][s];
		model.add(expr1 <= PP.K[t]);
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
		for (int s = 0; s < S-1; ++s) {
			IloExpr expr4(env);
			for (int r = 0; r < P; ++r)
				expr4 += zz[p][r][s+1];
			model.add(yy[p][s] == expr4);
		}
	
	cplex = IloCplex(model);

}

void GLSP::SetupModel_SO(double timeLimit) //Setup the standard model
{
	int P = PP.P; //The set of products, denoted by p
	int T = PP.T; //The set of macro - periods

	Matrix setup_pr;
	setup_pr.resize(PP.P);
	for (int p = 0; p < PP.P; ++p)
		setup_pr[p].resize(PP.P);

	for (int p = 0; p < PP.P; ++p)
		for (int r = 0; r < PP.P; ++r)
			setup_pr[p][r] = PP.Products[p].s_pr[r];

	ProductSet SeqSet;
	for (int p = 0; p < PP.P; ++p)
		SeqSet.insert(p);

	seq_DP.MinSetupMiddle(SOcache, setup_pr, -1, SeqSet, -1);
	map<CacheKey, int>::iterator it = SOcache.begin();

	while (it != SOcache.end()) {
		if (get<0>(it->first) == -1 || get<2>(it->first) == -1)
			SOcache.erase(it);

		++it;
	}

	for (int p = 0; p < PP.P; ++p){
		auto [it, inserted] = SOcache.try_emplace(CacheKey(p, { }, p));
		it->second = 0;
		for (int r = 0; r < PP.P; ++r) 
			if(p != r)	{
			auto [it, inserted] = SOcache.try_emplace(CacheKey(p, { }, r));
			it->second = setup_pr[p][r];
			}
	}

	int SO = SOcache.size(); //The set of sequences

	IntegerMatrix g;
	g.resize(P);
	for (int p = 0; p < P; ++p)
		g[p].resize(SO);

	IntegerMatrix f;
	f.resize(P);
	for (int p = 0; p < P; ++p)
		f[p].resize(SO);

	IntegerMatrix l;
	l.resize(P);
	for (int p = 0; p < P; ++p)
		l[p].resize(SO);

	vector<int> sc;
	sc.resize(SO);

	it = SOcache.begin();
	int s = 0;
	while (it != SOcache.end()) {
		sc[s] = it->second;
		f[get<0>(it->first)][s] = 1;
		l[get<2>(it->first)][s] = 1;
		g[get<0>(it->first)][s] = 1;
		g[get<2>(it->first)][s] = 1;
		for (auto i : get<1>(it->first))
			g[i][s] = 1;
		++it;
		++s;
	}

	vector<int> SO_t;
	SO_t.resize(T);

	for (int t = 0; t < T; ++t)
		SO_t[t] = SO;

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity); //Inventory level of product p at the end of macro-period t
	q = CreateNumVarArray2(env, P, T, "q", 0, IloInfinity); //Amount of product p produced in micro - period s

	w = CreateBoolVarArray2(env, SO, T, "w"); //1, if the machine is prepared for the product p in micro - period s; 0, otherwise
	
	// Constraint (50)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t) {
			double Demand = PP.Products[p].d[t];
			model.add(I[p][t] == I[p][t - 1] + q[p][t] - Demand);
		}

	// Constraint (50.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			double Demand = PP.Products[p].d[t];
			model.add(I[p][t] == q[p][t] - Demand);
		}

	// Constraint (58)
	for (int t = 0; t < T; ++t) {
		IloExpr Cap(env);
		for (int p = 0; p < P; ++p)
			Cap += PP.Products[p].a * q[p][t];
		model.add(Cap <= PP.K[t]);
	}

	// Constraint (59)
	for (int t = 0; t < T; ++t) {
		IloExpr Each(env);
		for (int s = 0; s < SO_t[t]; ++s)
			Each += w[s][t];
		model.add(Each == 1);
	}
	
	// Constraint (60)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t){
			IloExpr first(env);
			IloExpr last(env);
			for (int s = 0; s < SO_t[t]; ++s)
				first += f[p][s] * w[s][t];
			for (int s = 0; s < SO_t[t-1]; ++s)
				last += l[p][s] * w[s][t-1];
			model.add(first == last);
		}
	
	// Constraint (61)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t) {
			IloExpr Each2(env);
			for (int s = 0; s < SO_t[t]; ++s) {
				Each2 += g[p][s] * w[s][t];
			}
			model.add(q[p][t] * PP.Products[p].a <= PP.K[t] * Each2);
		}

	// Objective function terms
	// Setup cost
	IloExpr obj2(env);
	for (int t = 0; t < T; ++t)
		for (int s = 0; s < SO; ++s)
			obj2 += sc[s] * w[s][t];

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

bool GLSP::Solve(double timeLimit)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		//cplex.exportModel("GLSP.lp");
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

void GLSP::GetSolutions_SO(int P, int T, int S)
{
	for (int t = 0; t < T; ++t)
		for (int s = 0; s < SOcache.size(); ++s)
			if (cplex.getValue(w[s][t]) > 0)
				cout << "Value of w[" << s << "][" << t << "] is " << cplex.getValue(w[s][t]) << endl;

}