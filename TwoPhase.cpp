#include "Common.h"

typedef vector<vector<double>> Matrix;

#if 0

int main()
{
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

	// IP Model
	clock_t startIP = clock(), finishIP; //Start clock for IP

	IloEnv MPenv;
	IloModel MPmodel(MPenv);

	IloNumVar theta (MPenv, 0, IloInfinity);

	NumVarArray2 I = CreateNumVarArray2(MPenv, P, T, "I", 0, IloInfinity);
	NumVarArray2 q = CreateNumVarArray2(MPenv, P, T, "q", 0, IloInfinity);

	BoolVarArray2 x = CreateBoolVarArray2(MPenv, P, T, "y");

	IloExpr MPobj1(MPenv);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			MPobj1 += (f[p] * x[p][t] + c[p] * q[p][t]);

	IloExpr MPobj2(MPenv);
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			MPobj2 += h[p] * I[p][t];

	IloExpr total_MPobj(MPenv);
	total_MPobj = MPobj1 + MPobj2 + theta;

	IloObjective MPobj = IloMinimize(MPenv, total_MPobj);
	MPmodel.add(MPobj);

	// Constraint (2)
	for (int p = 0; p < P; ++p)
		for (int t = 1; t < T; ++t) {
			MPmodel.add(I[p][t] == I[p][t - 1] + q[p][t] - d[p][t]);
		}

	// Constraint (2.1)
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < 1; ++t) {
			MPmodel.add(I[p][t] == q[p][t] - d[p][t]);
		}

	// Constraint (3)
	for (int t = 0; t < T; ++t) {
		IloExpr Cap(MPenv);
		for (int p = 0; p < P; ++p)
			Cap += a[p] * q[p][t];
		MPmodel.add(Cap <= K[t]);
	}

	// Constraint (4)
	for (int t = 0; t < T; ++t)
		for (int p = 0; p < P; ++p)
			MPmodel.add(q[p][t] * a[p] <= K[t] * x[p][t]);

	IloCplex MPcplex(MPmodel);

	// Subproblem
	IloEnv SPenv;
	IloModel SPmodel(SPenv);

	int J = IloSum(x);

	vector<int> a;
	a.resize(J + 1);
	for (int j = 1; j < J; ++j)
		a[j] = 1;
	a[0] = 0;
	a[4] = 0;

	vector<int> b;
	b.resize(J + 1);
	//for (int j = 1; j < J; ++j)
		//b[j] = 1;
	b[0] = 0;
	b[1] = 1;
	b[2] = 2;
	b[3] = 4;
	b[4] = 99;

	vector<int> h;
	h.resize(J + 1);
	for (int j = 1; j < J; ++j)
		h[j] = 1;
	h[0] = 0;
	h[4] = 0;

	Matrix setup;
	setup.resize(J + 1);
	for (int j = 0; j < J + 1; ++j)
		setup[j].resize(J + 1);

	for (int j = 0; j < J + 1; ++j)
		for (int l = 0; l < J + 1; ++l) {
			setup[j][l] = 0;
		}

	IloNumVarArray C = CreateNumVarArray(SPenv, J + 1, "C", 0, IloInfinity);

	BoolVarArray2 e = CreateBoolVarArray2(SPenv, J + 1, J + 1, "e");

	IloExpr SPobj2(SPenv);
	for (int j = 0; j < J; ++j)
		for (int l = 1; l < J; ++l)
			if (j != l)
				SPobj2 += setup[j][l] * e[j][l];

	IloObjective SPobj = IloMinimize(SPenv, SPobj2);
	SPmodel.add(SPobj);

	IloCplex SPcplex(SPmodel);

	// Constraint (2)
	for (int j = 0; j < J; ++j) {
		IloExpr Each1(SPenv);
		for (int l = 1; l < J + 1; ++l)
			if (j != l)
				Each1 += e[j][l];
		SPmodel.add(Each1 == 1);
	}

	// Constraint (3)
	for (int j = 1; j < J + 1; ++j) {
		IloExpr Each2(SPenv);
		for (int l = 1; l < J + 1; ++l)
			if (j != l)
				Each2 += e[l][j];
		SPmodel.add(Each2 == 1);
	}

	// Constraint (4)
	for (int j = 0; j < J; ++j)
		for (int l = 1; l < J + 1; ++l)
			SPmodel.add(C[j] + q_hat[l] * a[l] <= C[l] + M * (1 - e[j][l]));

	// Constraint (5)
	for (int j = 1; j < J; ++j)
		SPmodel.add(C[j] <= b[j]);

	SPmodel.add(C[0] == 0);

	//Deactivate Presolve operations and use Primal 
	//SPcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	//SPcplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	//Benders Algorithm
	clock_t startB = clock(), finishB; //Start clock for Benders Decomposition
	int iter = 0;
	//int feas_cut = 0;
	int opt_cut = 0;
	double UB = 0;
	//double UB_best = w * Rectangles.size() + b_sum;// 9999999;
	while (true)
	{
		iter += 1;
		// Solve Master Problem
		MPcplex.exportModel("MP_GLSP.lp");
		bool result = MPcplex.solve();

		double theta_hat = MPcplex.getValue(theta);
		double x_hat = 0;
		for (int p = 0; p < P; ++p)
			for (int t = 0; t < T; ++t)
				x_hat += MPcplex.getValue(x[p][t]);

		if (result && MPcplex.isPrimalFeasible())
		{

			//Solve subproblem
			SPcplex.exportModel("SP_GLSP.lp");
			bool resultSP = SPcplex.solve();

			double nu_hat = SPcplex.getObjValue();

			cout << SPcplex.getCplexStatus() << endl;
			// Optimal solution
			if (nu_hat <= theta_hat)
			{
				finishB = clock(); //Stop clock for Benders Decomposition
				cout << "Optimal solution found " << endl;
				cout << MPcplex.getObjValue() << endl;
				cout << "Number of iterations " << iter << endl;
				cout << "Number of optimality cuts " << opt_cut << endl;
				cout << "Time Spent for Benders Decomposition = " << (float)(finishB - startB) / CLOCKS_PER_SEC << endl;
				cout << "Benders objective " << MPcplex.getObjValue() << endl;

				//cout << "IP bound " << IPcplex.getBestObjValue() << endl;
				//cout << "IP objective " << IPcplex.getObjValue() << endl;
				//cout << "IP gap " << IPcplex.getMIPRelativeGap() << endl;
				//cout << "Time Spent for IP = " << (float)(finishIP - startIP) / CLOCKS_PER_SEC << endl;

					break;
			}
			//Add optimality cut to Master Problem
			else 
			{
				opt_cut += 1;
				IloExpr Sum1(MPenv);
				for (int p = 0; p < P; ++p)
					for (int t = 0; t < T; ++t)
					{
						Sum1 += x[p][t];
					}
				IloExpr Sum2(MPenv);
				for (int p = 0; p < P; ++p)
					for (int t = 0; t < T; ++t)
					{
						Sum2 += x[p][t];
					}

				IloExpr OptCut(MPenv);
				OptCut = Sum1 - Sum2;

				MPmodel.add(theta >= theta_hat * OptCut - theta_hat * (F - 1));
			}
		}
		finishB = clock(); //Stop clock for Benders Decomposition
		// Time limit exceeded... STOP
		if ((float)(finishB - startB) / CLOCKS_PER_SEC > 900)
		{
			bool result1 = MPcplex.solve();
			cout << "Time limit exceeded..." << endl;
			cout << "Number of iterations " << iter << endl;
			cout << "Number of optimality cuts " << opt_cut << endl;
			cout << "Time Spent for Benders Decomposition = " << (float)(finishB - startB) / CLOCKS_PER_SEC << endl;
			cout << "Lower bound " << MPcplex.getObjValue() << endl;

			//cout << "IP bound " << IPcplex.getBestObjValue() << endl;
			//cout << "IP objective " << IPcplex.getObjValue() << endl;
			//cout << "IP gap " << IPcplex.getMIPRelativeGap() << endl;
			//cout << "Time Spent for IP = " << (float)(finishIP - startIP) / CLOCKS_PER_SEC << endl;

			break;
		}

	}
	return 0;
}
#endif