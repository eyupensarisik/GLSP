
#include "TwoPhase.h"
#include "DataPrep.h"
#include "setupDP.h"
#include "Subproblems.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;

TwoPhase::TwoPhase(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), model(env), cplex(env)
{
	BestBoundThreshold = DBL_MAX;
}

TwoPhase::~TwoPhase()
{
	env.end();
}

void TwoPhase::SetupModel()
{
	// Setup Master Problem
	Matrix setup_pr;
	setup_pr.resize(PP.P);
	for (int p = 0; p < PP.P; ++p)
		setup_pr[p].resize(PP.P);

	for (int p = 0; p < PP.P; ++p)
		for (int r = 0; r < PP.P; ++r)
			setup_pr[p][r] = PP.Products[p].s_pr[r];

	theta = IloNumVar(env, 0, IloInfinity, "theta");
	gamma = IloNumVar(env, 0, IloInfinity, "gamma");
	theta_t = CreateNumVarArray(env, PP.T, "theta_t", 0, IloInfinity);

	I = CreateNumVarArray2(env, PP.P, PP.T, "I", 0, IloInfinity);
	q = CreateNumVarArray2(env, PP.P, PP.T, "q", 0, IloInfinity);

	x = CreateBoolVarArray2(env, PP.P, PP.T, "x");

	theta.setStringProperty("Type", "theta");
	theta.setIntProperty("Product", -1);

	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < PP.T; ++t)
		{
			x[p][t].setStringProperty("Type", "x");
			x[p][t].setIntProperty("Product", p);
		}
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < PP.T; ++t) {
			model.add(I[p][t]);
			model.add(q[p][t]);
			model.add(x[p][t]);
		}
	model.add(theta);
	model.add(gamma);

	for (int t = 0; t < PP.T; ++t)
		model.add(theta_t[t]);

	IloExpr MPobj1(env);
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < PP.T; ++t)
			MPobj1 += (PP.Products[p].f * x[p][t] + PP.Products[p].c * q[p][t]);

	IloExpr MPobj2(env);
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < PP.T; ++t)
			MPobj2 += PP.Products[p].h * I[p][t];

	IloExpr total_MPobj(env);
	total_MPobj = MPobj1 + MPobj2 + theta;

	obj = IloMinimize(env, total_MPobj);
	model.add(obj);

	// Constraint (2)
	for (int p = 0; p < PP.P; ++p)
		for (int t = 1; t < PP.T; ++t) {
			model.add(I[p][t] == I[p][t - 1] + q[p][t] - PP.Products[p].d[t]);
		}

	// Constraint (2.1)
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < 1; ++t) {
			model.add(I[p][t] == q[p][t] - PP.Products[p].d[t]);
		}

	// Constraint (3)
	for (int t = 0; t < PP.T; ++t) {
		IloExpr Cap(env);
		for (int p = 0; p < PP.P; ++p)
			Cap += PP.Products[p].a * q[p][t];
		model.add(Cap <= PP.K[t]);
	}

	// Constraint (4)
	for (int t = 0; t < PP.T; ++t)
		for (int p = 0; p < PP.P; ++p)
			model.add(q[p][t] * PP.Products[p].a <= PP.K[t] * x[p][t]);

	Subproblems SP(PP, Parameters);

	int W_In = PP.P;

	vector<int> wp_In;
	wp_In.resize(W_In + 2);
	vector<int> wt_In;
	wt_In.resize(W_In + 2);
	vector<int> wop_In;
	wop_In.resize(PP.P - W_In + 2);
	vector<int> wot_In;
	wot_In.resize(PP.P - W_In + 2);
	int j = 0;
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < 1; ++t) {
			j += 1;
			wp_In[j] = p;
			wt_In[j] = t;
		}

	SP.SetupBSPModel(W_In, wp_In, wt_In, wop_In, wot_In);
	SP.BSP_Solve();

	LB_theta = SP.GetBSP_UB();
	model.add(theta >= LB_theta);

	double ValIn_setup = 0;
	IloExpr ValIn(env);

	setupDP sDP;
	map <pair<int, set<int>>, int> Cache;
	set<int> current;
	for (int t = 0; t < PP.T; ++t){
		for (int p = 0; p < PP.P; ++p)
			if (PP.Products[p].d[t] > 0)
				current.insert(p);

		if (!current.empty()) {
			cout << "Finding minimum setup" << endl;
			ValIn_setup = sDP.MinSetup(Cache, setup_pr, -1, current);
			cout << "Done. Minimum setup: " << ValIn_setup << endl;
		}

		ValIn += theta_t[t];

		model.add(ValIn >= ValIn_setup);
	}

	IloExpr thetaCut(env);
	for (int t = 0; t < PP.T; ++t)
		thetaCut += theta_t[t];

	model.add(theta >= thetaCut + gamma);

	cplex = IloCplex(model);

}

bool TwoPhase::Solve(double timeLimit, double* TwoPhase_Iter, double* TwoPhase_Cut, double* TwoPhase_CPU, double* TwoPhase_UB, double* TwoPhase_LB, double* SP_Cons_CPU, double* SP_Solve_CPU, double* MP_CPU)
{	
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		//Algorithm

		auto startTP = chrono::high_resolution_clock::now();

		double iter = 0;
		double opt_cut = 0;
		double UB = 0;
		double UB_best = INT_MAX;
		double SPtime = 0;
		double TotalSPtime = 0;
		double TotalSPcons_time = 0;

		while (true)
		{
			iter += 1;

			// Solve Master Problem
			cplex.exportModel("MP_GLSP.lp");
			bool result = cplex.solve();

			if (result && cplex.isPrimalFeasible())
			{
				double theta_hat = cplex.getValue(theta);
				int W = 0;
				for (int p = 0; p < PP.P; ++p) {
					for (int t = 0; t < PP.T; ++t) {
						if (cplex.getValue(x[p][t]) >= 1 - 0.0001) {
							W += 1;
						}
					}
				}

				Matrix x_val;
				x_val.resize(PP.P);
				for (int p = 0; p < PP.P; ++p)
					x_val[p].resize(PP.T);

				for (int p = 0; p < PP.P; ++p) {
					for (int t = 0; t < PP.T; ++t) {
						if (cplex.getValue(x[p][t]) >= 1 - 0.0001) {
							x_val[p][t] = 1;
						}
						else
							x_val[p][t] = 0;
					}
				}

				vector<int> q_hat;
				q_hat.resize(W + 2);
				int j = 0;
				for (int p = 0; p < PP.P; ++p)
					for (int t = 0; t < PP.T; ++t)
						if (x_val[p][t] >= 1 - 0.0001) {
							j += 1;
							q_hat[j] = cplex.getValue(q[p][t]);
						}

				//Construct Subproblem
				auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP

				vector<int> wp;
				wp.resize(W + 2);
				vector<int> wt;
				wt.resize(W + 2);
				vector<int> wop;
				wop.resize(PP.P * PP.T - W + 2);
				vector<int> wot;
				wot.resize(PP.P * PP.T - W + 2);
				j = 0;
				int k = 0;
				for (int p = 0; p < PP.P; ++p)
					for (int t = 0; t < PP.T; ++t)
						if (x_val[p][t] >= 1 - 0.0001) {
							j += 1;
							wp[j] = p;
							wt[j] = t;
						}
						else {
							k += 1;
							wop[k] = p;
							wot[k] = t;
						}

				Subproblems SP(PP, Parameters);

				SP.SetupBSPModel(W, wp, wt, wop, wot);
				

				auto finishSP_cons = chrono::high_resolution_clock::now();
				auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
				TotalSPcons_time += SPcons_time.count();

				auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP
				SP.BSP_Solve();

				auto finishSP = chrono::high_resolution_clock::now();
				auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
				TotalSPtime += SPtime.count();

				double nu_hat = SP.GetBSP_UB();
				SP.GetBSP_Solutions(W);

				UB = cplex.getObjValue() - theta_hat + nu_hat;
				if (UB <= UB_best)
					UB_best = UB;

				for (int t = 0; t < PP.T; ++t)
					cout << cplex.getValue(theta_t[t]) << endl;
				cout << cplex.getValue(gamma) << endl;

				// Optimal solution
				if (nu_hat <= theta_hat + 0.0001)
				{
					Solved = 1;
					cplex.exportModel("MP_GLSP.lp");
					for (int t = 0; t < PP.T; ++t)
						cout << cplex.getValue(theta_t[t]) << endl;
					cout << cplex.getValue(gamma) << endl;
					auto finishTP = chrono::high_resolution_clock::now();
					auto TPtime = chrono::duration_cast<chrono::milliseconds>(finishTP - startTP);
					cout << "Optimal solution found " << endl;
					cout << cplex.getObjValue() << endl;
					cout << "Number of iterations " << iter << endl;
					cout << "Number of optimality cuts " << opt_cut << endl;
					cout << "Time Spent for TwoPhase Decomposition = " << (float)TPtime.count() / 1000 << endl;
					cout << "Time Spent for SP Construction = " << (float)TotalSPcons_time / 1000 << endl;
					cout << "Time Spent for SP Solve = " << (float)TotalSPtime / 1000 << endl;
					cout << "TwoPhase objective " << cplex.getObjValue() << endl;
					cout << "Value of theta is " << cplex.getValue(theta) << endl;

					for (int j = 1; j < W + 1; ++j)
						for (int l = 1; l < W + 1; ++l) {
							if (SP.e_val[j][l] > 0.5)
								cout << "Value of e[" << j << "][" << l << "] is " << SP.e_val[j][l] << endl;
						}

					for (int p = 0; p < PP.P; ++p) {
						for (int t = 0; t < PP.T; ++t) {
							if (x_val[p][t] >= 0.5)
								cout << "Value of x[" << p << "][" << t << "] is " << x_val[p][t] << endl;
						}
					}

					*TwoPhase_Iter = iter;
					*TwoPhase_Cut = opt_cut;
					*TwoPhase_CPU = (float)TPtime.count() / 1000;
					*TwoPhase_UB = cplex.getObjValue();
					*TwoPhase_LB = cplex.getBestObjValue();
					*SP_Cons_CPU = (float)TotalSPcons_time / 1000;
					*SP_Solve_CPU = (float)TotalSPtime / 1000;
					*MP_CPU = *TwoPhase_CPU - *SP_Cons_CPU - *SP_Solve_CPU;

					break;
				}
				//Add optimality cut to Master Problem
				else
				{
					opt_cut += 1;

					vector<int> SP_Sol;
					SP_Sol.resize(W);

					vector<int> SP_Sol2;

					int j = 0;
					int i = 0;
					for (int l = 1; l < W + 1; ++l)
						if (SP.e_val[j][l] > 0.5) {
							cout << "Value of e[" << j << "][" << l << "] is " << SP.e_val[j][l] << endl;
							j = l;
							l = 0;
							SP_Sol[i] = j;
							i += 1;
						}

					int LB = 0;
					for (int i = 0; i < W - 1; ++i)
						LB += SP.SP_setup[SP_Sol[i]][SP_Sol[i + 1]];

					for (int k = 0; k < W; ++k) {
						SP_Sol2 = SP_Sol;
						SP_Sol2.erase(SP_Sol2.begin() + k);
						int LBs = 0;
						for (int i = 0; i < W - 2; ++i)
							LBs += SP.SP_setup[SP_Sol2[i]][SP_Sol2[i + 1]];

						if (LBs < LB)
							LB = LBs;
					}

					int set_mins = 0;
					int set_min = 0;
					vector<int> set_min_p;
					set_min_p.resize(PP.P);

					for (int p = 0; p < PP.P; ++p) {
						int set_mins = 999999;
						int set_min = 999999;
						for (int r = 0; r < PP.P; ++r) {
							if (p != r)
								set_mins = PP.Products[p].s_pr[r];

							if (set_mins < set_min)
								set_min = set_mins;
						}
						set_min_p[p] = set_min;
					}

					vector<int> set_min_total;
					set_min_total.resize(W);

					for (int k = 0; k < W; ++k) {
						int kk = 0;
						for (int t = 0; t < PP.T; ++t) {
							int set_min_max = 0;
							for (int p = 0; p < PP.P; ++p) {
								if (cplex.getValue(x[p][t]) >= 0.5) {
									if (kk != k) {
										set_min_total[k] += set_min_p[p];

										if (set_min_p[p] > set_min_max)
											set_min_max = set_min_p[p];
									}
									kk += 1;
								}
							}
							set_min_total[k] -= set_min_max;
						}
					}

					LB = *min_element(set_min_total.begin(), set_min_total.end());

					//DP implementation
					setupDP sDP;
					map <pair<int, set<int>>, int> Cache;
					vector<double> LBsetup;
					LBsetup.resize(PP.T);

					Matrix LBsetupPer;
					LBsetupPer.resize(PP.P);
					for (int p = 0; p < PP.P; ++p)
						LBsetupPer[p].resize(PP.T);

					int LBsetupDP = 0;

					int per = 0;
					int pre_per = 0;

					for (int t = 0; t < PP.T; ++t) {
						for (int p = 0; p < PP.P; ++p) {
							if (x_val[p][t] >= 0.5)
								per += 1;
						}

						set<int> current;
						for (int i = pre_per; i < per; ++i)
							current.insert(SP_Sol[i]);

						if (!current.empty()) {
							cout << "Finding minimum setup" << endl;
							LBsetupDP = sDP.MinSetup(Cache, SP.SP_setup, -1, current);
							cout << "Done. Minimum setup: " << LBsetupDP << endl;
						}
						else
							LBsetupDP = 0;

						LBsetup[t] = LBsetupDP;

						auto currentCopy = current;
					if (pre_per != per){
						for (int i = pre_per; i < per; ++i) {
							current.erase(SP_Sol[i]);

							if (!current.empty()) {
								cout << "Finding minimum setup" << endl;
								LBsetupDP = sDP.MinSetup(Cache, SP.SP_setup, -1, current);
								cout << "Done. Minimum setup: " << LBsetupDP << endl;
							}
							else
								LBsetupDP = 0;

							current = currentCopy;

							LBsetupPer[wp[SP_Sol[i]]][t] = LBsetupDP;
						}
					}
					else
						for (int p = 0; p < PP.P; ++p)
							LBsetupPer[p][t] = 0;
						pre_per = per;
					}

					vector<int> LBsetupPermin;
					LBsetupPermin.resize(PP.T);
					for (int t = 0; t < PP.T; ++t)
						LBsetupPermin[t] = INT_MAX;

					for (int t = 0; t < PP.T; ++t)
						for (int p = 0; p < PP.P; ++p)
							if (LBsetupPer[p][t] < LBsetupPermin[t])
								LBsetupPermin[t] = LBsetupPer[p][t];

					if (LB < LB_theta)
						LB = LB_theta;

					IloExpr Sum1(env);
					for (int j = 1; j < W + 1; ++j)
						Sum1 += x[wp[j]][wt[j]];
					cout << Sum1 << endl;

					IloExpr Sum2(env);
					for (int k = 1; k < (PP.P * PP.T - W + 1); ++k)
						Sum2 += x[wop[k]][wot[k]];
					cout << Sum2 << endl;

					IloExpr OptCut(env);
					OptCut = Sum1 - Sum2;
					cout << OptCut << endl;

					for (int t = 0; t < PP.T; ++t) {
						IloExpr DPCut(env);
						for (int p = 0; p < PP.P; ++p) {
							if (x_val[p][t] >= 0.5) 
								DPCut += (LBsetup[t]) * (1 - x[p][t]);
						}
						model.add(theta_t[t] >= LBsetup[t] - DPCut);
					}

					double LBsetuptotal = 0;
					for (int t = 0; t < PP.T; ++t)
						LBsetuptotal += LBsetup[t];

					double LBresidual = nu_hat - LBsetuptotal;

					IloExpr GMCut(env);
					for (int t = 0; t < PP.T; ++t) {
						for (int p = 0; p < PP.P; ++p) {
							if (x_val[p][t] >= 0.5) 
								GMCut += LBresidual * (1 - x[p][t]);
							else
								GMCut += LBresidual * (x[p][t]);
						}
					}
					
					model.add(gamma >= LBresidual - GMCut);

					model.add(theta >= (nu_hat - LB) * Sum1 - (nu_hat - LB) * (W - 1) + LB);

				}
				SP.SPmodel.end();
			}
			else {
				Solved = 1;
				cout << "The problem is infeasible..." << endl;
			}
				
			auto finishB = chrono::high_resolution_clock::now();
			auto Btime = chrono::duration_cast<chrono::milliseconds>(finishB - startTP);
			// Time limit exceeded... STOP
			if (Btime.count() / 1000 > timeLimit)
			{
				bool result1 = cplex.solve();
				cout << "Time limit exceeded..." << endl;
				cout << "Number of iterations " << iter << endl;
				cout << "Number of optimality cuts " << opt_cut << endl;
				cout << "Time Spent for TwoPhase Decomposition = " << (float)Btime.count() / 1000 << endl;
				cout << "Lower bound " << cplex.getBestObjValue() << endl;
				cout << "Upper bound " << UB_best << endl;

				*TwoPhase_Iter = iter;
				*TwoPhase_Cut = opt_cut;
				*TwoPhase_CPU = (float)Btime.count() / 1000;
				*TwoPhase_UB = UB_best;
				*TwoPhase_LB = cplex.getBestObjValue();
				*SP_Cons_CPU = (float)TotalSPcons_time / 1000;
				*SP_Solve_CPU = (float)TotalSPtime / 1000;
				*MP_CPU = *TwoPhase_CPU - *SP_Cons_CPU - *SP_Solve_CPU;

				break;
			}
		}
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return Solved;
}

double TwoPhase::GetLB()
{
	return Solved ? cplex.getObjValue() : 0;
}

double TwoPhase::GetUB()
{
	return Solved ? cplex.getBestObjValue() : DBL_MAX;
}