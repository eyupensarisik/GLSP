#include "Common.h"

typedef vector<vector<double>> Matrix;
#if 0



int main()
{
	IloEnv env;
	IloModel model(env);

	int J = 4;
	int M = 99;

	vector<int> a;
	a.resize(J+1);
	for (int j = 1; j < J; ++j)
		a[j] = 1;
	a[0] = 0;
	a[4] = 0;

	vector<int> b;
	b.resize(J+1);
	//for (int j = 1; j < J; ++j)
		//b[j] = 1;
	b[0] = 0;
	b[1] = 1;
	b[2] = 2;
	b[3] = 4;
	b[4] = 99;

	vector<int> h;
	h.resize(J+1);
	for (int j = 1; j < J; ++j)
		h[j] = 1;
	h[0] = 0;
	h[4] = 0;

	Matrix setup;
	setup.resize(J+1);
	for (int j = 0; j < J+1; ++j)
		setup[j].resize(J+1);

	for (int j = 0; j < J+1; ++j)
		for (int l = 0; l < J+1; ++l) {
			setup[j][l] = 0;
		}

	IloNumVarArray C = CreateNumVarArray(env, J + 1, "C", 0, IloInfinity);

	BoolVarArray2 e = CreateBoolVarArray2(env, J + 1, J + 1, "e");

	IloExpr obj1(env);
	for (int j = 1; j < J; ++j)
		obj1 += (h[j] * (a[j] * (b[j] - C[j])));

	IloExpr obj2(env);
	for (int j = 0; j < J; ++j)
		for (int l = 1; l < J; ++l)
			if (j != l)
				obj2 += setup[j][l] * e[j][l];

	IloExpr total_obj(env);
	total_obj = obj1 + obj2;

	IloObjective obj = IloMinimize(env, total_obj);
	model.add(obj);

	// Constraint (2)
	for (int j = 0; j < J; ++j) {
			IloExpr Each1(env);
			for (int l = 1; l < J+1; ++l)
				if (j != l)
					Each1 += e[j][l];
			model.add(Each1 == 1);
		}

	// Constraint (3)
	for (int j = 1; j < J+1; ++j) {
		IloExpr Each2(env);
		for (int l = 1; l < J + 1; ++l)
			if (j != l)
				Each2 += e[l][j];
		model.add(Each2 == 1);
	}

	// Constraint (4)
	for (int j = 0; j < J; ++j)
		for (int l = 1; l < J + 1; ++l)
			model.add(C[j] + a[l] <= C[l] + M * (1 - e[j][l]));

	// Constraint (5)
	for (int j = 1; j < J; ++j)
		model.add(C[j] <= b[j]);

	model.add(C[0] == 0);

	
	IloCplex cplex(model);
	cplex.exportModel("BSP.lp");

	bool result = cplex.solve();

	if (result && cplex.isPrimalFeasible())
	{
		IloNum UB = cplex.getObjValue();
		IloNum LB = cplex.getBestObjValue();

		cout << "UB: " << UB << " LB: " << LB << endl;

		for (int j = 1; j < J; ++j)
			for (int l = 1; l < J; ++l) {
				IloNum e_val = cplex.getValue(e[j][l]);
				cout << "e_ " << j << l << " : " << e_val << endl;
			}

		for (int j = 0; j < J + 1; ++j) {
				IloNum C_val = cplex.getValue(C[j]);
				cout << "C_ " << j << " : " << C_val << endl;
			}
	}

	return 0;
}
#endif // 1