#include "Common.h"

typedef vector<vector<double>> Matrix;
#if 0

int main()
{
	IloEnv env;
	IloModel model(env);

	int P = 2;
	int T = 2;
	int S = 4;
	int M = 99;

	vector<int> S_b;
	S_b.resize(T);
	//for (int t = 0; t < T; ++t)
		//S_b[t] = 1;
	S_b[0] = 0;
	S_b[1] = 2;

	vector<int> S_f;
	S_f.resize(T);
	//for (int t = 0; t < T; ++t)
		//S_f[t] = 1;
	S_f[0] = 2;
	S_f[1] = 4;

	vector<int> L;
	L.resize(T);
	//for (int t = 0; t < T; ++t)
		//L[t] = 1;
	L[0] = 1;
	L[1] = 3;

	vector<int> K;
	K.resize(T);
	for (int t = 0; t < T; ++t)
		K[t] = 2;

	vector<int> a;
	a.resize(P);
	for (int p = 0; p < P; ++p)
		a[p] = 1;

	vector<int> m;
	m.resize(P);
	//for (int p = 0; p < P; ++p)
		//m[p] = 1;
	m[0] = 1;
	m[1] = 2;

	vector<int> h;
	h.resize(P);
	for (int p = 0; p < P; ++p)
		h[p] = 1;

	Matrix s_pr;
	s_pr.resize(P);
	for (int p = 0; p < P; ++p)
		s_pr[p].resize(P);

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r) {
			s_pr[p][r] = 0;
		}

	Matrix d;
	d.resize(P);
	for (int p = 0; p < P; ++p)
		d[p].resize(T);

	//for (int p = 0; p < P; ++p)
		//for (int t = 0; t < T; ++t) {
			//d[p][t] = rand() % 11;
		//}
	d[0][0] = 1;
	d[0][1] = 0;
	d[1][0] = 0;
	d[1][1] = 3;

	vector<int> f;
	f.resize(P);
	for (int p = 0; p < P; ++p)
		f[p] = 0;

	vector<int> c;
	c.resize(P);
	for (int p = 0; p < P; ++p)
		c[p] = 0;

	NumVarArray2 I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity);
	NumVarArray2 q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity);

	BoolVarArray2 y = CreateBoolVarArray2(env, P, S, "y");
	BoolVarArray3 z = CreateBoolVarArray3(env, P, P, S, "z");

	IloExpr obj1(env);
	for (int p = 0; p < P; ++p)
		for (int s = 0; s < S; ++s)
			obj1 += (f[p] * y[p][s] + c[p] * q[p][s]);

	IloExpr obj2(env);
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 0; s < S; ++s)
				obj2 += s_pr[p][r] * z[p][r][s];

	IloExpr obj3(env);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			obj3 += h[p] * I[p][t];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2 + obj3;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	// Constraint (2)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t){
			IloExpr InvBal(env);
			for (int s = S_b[t]; s < S_f[t]; ++s)
				InvBal += q[p][s];
			model.add(I[p][t] == I[p][t-1] + InvBal - d[p][t]);
		}

	// Constraint (2.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			IloExpr InvBal2(env);
			for (int s = S_b[t]; s < S_f[t]; ++s)
				InvBal2 += q[p][s];
			model.add(I[p][t] == InvBal2 - d[p][t]);
		}
	
	// Constraint (3)
	for (int t = 0; t < T; ++t) {
		IloExpr Cap(env);
		for (int p = 0; p < P; ++p)
			for (int s = S_b[t]; s < S_f[t]; ++s)
				Cap += a[p] * q[p][s];
			model.add(Cap <= K[t]);
		}
	
	// Constraint (4)
	for (int t = 0; t < T; ++t) 
		for (int p = 0; p < P; ++p)
			for (int s = 0; s < S; ++s)
				model.add(q[p][s] * a[p] <= K[t] * y[p][s]);

	// Constraint (5)
	for (int s = 0; s < S; ++s) {
		IloExpr Each(env);
		for (int p = 0; p < P; ++p)
			Each += y[p][s];
		model.add(Each == 1);
	}

	// Constraint (6.1)
	/*for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			for (int s = 0; s < S; ++s)
				if (s == 0)
					model.add(q[p][s] >= m[p] * (y[p][s]));
				else
					model.add(q[p][s] >= m[p] * (y[p][s] - y[p][s - 1])); */

	// Constraint (6.2)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			for (int s = 0; s < S; ++s)
				if (s == 0)
					model.add(q[p][s] >= m[p] * (y[p][s]));
				else if (s == 1)
					model.add(q[p][s] + q[p][s + 1] >= m[p] * (y[p][s] - y[p][s - 1]));
				else
					model.add(q[p][s] >= m[p] * (y[p][s] - y[p][s - 1]));
					

	// Constraint (7)
	/*for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			for (int s = 1; s < S; ++s)
				model.add(z[p][r][s] >= y[p][s-1] + y[r][s] - 1);*/

	// Constraint (7.1)
	for (int p = 0; p < P; ++p)
		for (int s = 1; s < S; ++s) {
			IloExpr NF1(env);
			for (int r = 0; r < P; ++r)
				NF1 += z[p][r][s];
			model.add(NF1 == y[p][s-1]);
		}

	// Constraint (7.2)
	for (int r = 0; r < P; ++r)
		for (int s = 0; s < S; ++s) {
			IloExpr NF2(env);
			for (int p = 0; p < P; ++p)
				NF2 += z[p][r][s];
			model.add(NF2 == y[r][s]);
		}
		

	IloCplex cplex(model);
	cplex.exportModel("GLSP.lp");

	bool result = cplex.solve();

	if (result && cplex.isPrimalFeasible())
	{
		IloNum UB = cplex.getObjValue();
		IloNum LB = cplex.getBestObjValue();

		cout << "UB: " << UB << " LB: " << LB << endl;

		for (int p = 0; p < P; ++p)
			for (int t = 0; t < T; ++t) {
				IloNum I_val = cplex.getValue(I[p][t]);
				cout << "I_ " << p << t << " : " << I_val << endl;
			}

		for (int p = 0; p < P; ++p)
			for (int s = 0; s < S; ++s) {
				IloNum q_val = cplex.getValue(q[p][s]);
				cout << "q_ " << p << s << " : " << q_val << endl;
			}
	}

	return 0;
}
#endif // 0