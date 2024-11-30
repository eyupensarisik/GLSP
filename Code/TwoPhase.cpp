
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

void TwoPhase::SetupModel(double timeLimit)
{
	// Setup Master Problem
	Matrix setup_pr;
	setup_pr.resize(PP.P);
	for (int p = 0; p < PP.P; ++p)
		setup_pr[p].resize(PP.P);

	for (int p = 0; p < PP.P; ++p)
		for (int r = 0; r < PP.P; ++r)
			setup_pr[p][r] = PP.Products[p].s_pr[r];

	int TotalCap = 0;
	for (int t = 0; t < PP.T; ++t)
		TotalCap += PP.K[t];

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

	vector<ProductSet> ProductPPs_In;

	for (int t = 0; t < PP.T; ++t) {
		ProductSet MPcurrent;
		for (int p = 0; p < PP.P; ++p) {
			MPcurrent.insert(p);
		}
		ProductPPs_In.push_back(MPcurrent);
	}

	vector<ProductSet> ProductPPsCopy;
	for (int i = 0; i < ProductPPs_In.size(); ++i)
	{
		if (!ProductPPs_In[i].empty())
			ProductPPsCopy.push_back(ProductPPs_In[i]);
	}
	ProductPPs_In = ProductPPsCopy;

	cout << "Finding minimum setup multiperiod" << endl;
	cout << "Done. Minimum setup: " << sDP.MinSetupMultiPeriod(MPcache, setup_pr, -1, ProductPPs_In, -1) << endl;

	if (!ProductPPs_In[0].empty())
		LB_theta = sDP.MinSetupMiddle(MPcache, setup_pr, -1, ProductPPs_In[0], -1);

	model.add(theta >= LB_theta);

	double ValIn_setup = 0;
	IloExpr ValIn(env);

	set<int> current;
	for (int t = 0; t < PP.T; ++t) {
		for (int p = 0; p < PP.P; ++p)
			if (PP.Products[p].d[t] > 0)
				current.insert(p);

		if (!current.empty()) {
			ValIn_setup = sDP.MinSetupMiddle(MPcache, setup_pr, -1, current, -1);
		}

		ValIn += theta_t[t];

		model.add(ValIn >= ValIn_setup);
	}

	IloExpr thetaCut(env);
	for (int t = 0; t < PP.T; ++t)
		thetaCut += theta_t[t];

	model.add(theta >= thetaCut + gamma);

	for (int t = 1; t < PP.T; ++t)
		for (int p = 0; p < PP.P; ++p)
			model.add(I[p][t - 1] >= PP.Products[p].d[t] * (1 - x[p][t]));

	for (int p = 0; p < PP.P; ++p) {
		for (int t = 1; t < PP.T; ++t) {
			double eta;
			double alpha;
			for (int k = t; k < PP.T; ++k) {
				IloExpr CapVal(env);
				int d_pkt = 0;
				for (int j = t; j < k + 1; ++j) {
					CapVal += x[p][j];
					d_pkt += PP.Products[p].a * PP.Products[p].d[j];
				}
				//model.add(I[p][t - 1] + TotalCap * CapVal >= d_pkt);
				eta = ceil(d_pkt / TotalCap);
				alpha = d_pkt - TotalCap * (eta - 1);
				model.add(I[p][t - 1] >= alpha * (eta - CapVal));
			}
		}
	}

	SPtimelimit = timeLimit;

	cplex = IloCplex(model);

}

bool TwoPhase::Solve(double timeLimit, double* TwoPhase_Iter, double* TwoPhase_Cut, double* TwoPhase_CPU, double* TwoPhase_UB, double* TwoPhase_LB, double* SP_Cons_CPU, double* SP_Solve_CPU, double* MP_CPU, int SPtype)
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
		double TotalSPtime = 0;
		double TotalSPcons_time = 0;

		Matrix setup_pr;
		setup_pr.resize(PP.P);
		for (int p = 0; p < PP.P; ++p)
			setup_pr[p].resize(PP.P);

		for (int p = 0; p < PP.P; ++p)
			for (int r = 0; r < PP.P; ++r)
				setup_pr[p][r] = PP.Products[p].s_pr[r];

		setupDP sDP;
		map <pair<int, set<int>>, int> Cache;
		CacheMP MPcache;

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
				vector<int> W_t;
				W_t.resize(PP.T);
				for (int t = 0; t < PP.T; ++t) {
					int w = 0;
					for (int p = 0; p < PP.P; ++p) {
						if (cplex.getValue(x[p][t]) >= 1 - 0.0001) {
							W += 1;
							w += 1;
						}
					}
					W_t[t] = w;
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
				
				double nu_hat = 0;

				if (SPtype == 0){
					//Construct Subproblem
					auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP

					Subproblems SP(PP, Parameters);

					SP.SetupBSPModel(W, wp, wt, wop, wot);

					auto finishSP_cons = chrono::high_resolution_clock::now();
					auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
					TotalSPcons_time += SPcons_time.count();

					auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP
					SP.BSP_Solve(SPtimelimit);

					auto finishSP = chrono::high_resolution_clock::now();
					auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
					TotalSPtime += SPtime.count();

					nu_hat = SP.GetBSP_UB();
					SP.GetBSP_Solutions(W);

					vector<ProductSet> ProductPPs;

					for (int t = 0; t < PP.T; ++t) {
						ProductSet MPcurrent;
						for (int p = 0; p < PP.P; ++p) {
							if (x_val[p][t] >= 1 - 0.0001) {
								MPcurrent.insert(p);
							}
						}
						ProductPPs.push_back(MPcurrent);
					}
					vector<ProductSet> ProductPPsCopy;
					for (int i = 0; i < ProductPPs.size(); ++i)
					{
						if (!ProductPPs[i].empty())
							ProductPPsCopy.push_back(ProductPPs[i]);
					}
					ProductPPs = ProductPPsCopy;

					cout << "Finding minimum setup multiperiod" << endl;
					int nu_hatt = sDP.MinSetupMultiPeriod(MPcache, setup_pr, -1, ProductPPs, -1);
					if (nu_hatt - nu_hat < -1)
						return 0;

				}
				else if (SPtype == 1){
					auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP
					
					vector<ProductSet> ProductPPs;

					for (int t = 0; t < PP.T; ++t){
						ProductSet MPcurrent;
						for (int p = 0; p < PP.P; ++p){
							if (x_val[p][t] >= 1 - 0.0001){
								MPcurrent.insert(p);
							}
						}
						ProductPPs.push_back(MPcurrent);
					}

					auto finishSP_cons = chrono::high_resolution_clock::now();
					auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
					TotalSPcons_time += SPcons_time.count();

					auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP

					vector<ProductSet> ProductPPsCopy;
					for (int i = 0; i < ProductPPs.size(); ++i)
					{
						if (!ProductPPs[i].empty())
							ProductPPsCopy.push_back(ProductPPs[i]);
					}
					ProductPPs = ProductPPsCopy;

					nu_hat = sDP.MinSetupMultiPeriod(MPcache, setup_pr, -1, ProductPPs, -1);

					auto finishSP = chrono::high_resolution_clock::now();
					auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
					TotalSPtime += SPtime.count();
				}
				else if (SPtype == 2) {
					//Construct Subproblem
					auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP

					Subproblems CP(PP, Parameters);

					CP.SetupCPModel(W, wp, wt, wop, wot);

					auto finishSP_cons = chrono::high_resolution_clock::now();
					auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
					TotalSPcons_time += SPcons_time.count();

					auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP
					CP.CP_Solve(SPtimelimit);

					auto finishSP = chrono::high_resolution_clock::now();
					auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
					TotalSPtime += SPtime.count();

					nu_hat = CP.GetCP_Obj();
				}

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
					int LB = 0;
					
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
					
					vector<double> LBsetup;
					LBsetup.resize(PP.T);

					Matrix LBsetupPer;
					LBsetupPer.resize(PP.P);
					for (int p = 0; p < PP.P; ++p)
						LBsetupPer[p].resize(PP.T);

					int LBsetupDP = 0;

					for (int t = 0; t < PP.T; ++t) {
						set<int> current;
						for (int p = 0; p < PP.P; ++p) {
							if (x_val[p][t] >= 0.5)
								current.insert(p);
						}

						if (!current.empty()) {
							LBsetupDP = sDP.MinSetupMiddle(MPcache, setup_pr, -1, current, -1);
						}
						else
							LBsetupDP = 0;

						LBsetup[t] = LBsetupDP;

						auto currentCopy = current;

						if (!current.empty()) {
							for (int p = 0; p < PP.P; ++p) {
								if (x_val[p][t] >= 0.5)
									current.erase(p);

								if (!current.empty()) {
									LBsetupDP = sDP.MinSetupMiddle(MPcache, setup_pr, -1, current, -1);
								}
								else
									LBsetupDP = 0;

								current = currentCopy;

								LBsetupPer[p][t] = LBsetupDP;
							}
						}
						else
							for (int p = 0; p < PP.P; ++p)
								LBsetupPer[p][t] = 0;
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
					IloExpr Sum2(env);
					for (int p = 0; p < PP.P; ++p)
						for (int t = 0; t < PP.T; ++t)
							if (x_val[p][t] >= 0.5)
								Sum1 += x[p][t];
							else
								Sum2 += x[p][t];

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
				//SP.SPmodel.end();
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