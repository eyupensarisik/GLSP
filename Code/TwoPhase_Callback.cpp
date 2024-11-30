
#include "TwoPhase_Callback.h"
#include "DataPrep.h"
#include "setupDP.h"
#include "Subproblems.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;

TwoPhaseC::TwoPhaseC(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), model(env), cplex(env)
{
	BestBoundThreshold = DBL_MAX;
}

TwoPhaseC::~TwoPhaseC()
{
	env.end();
}

void TwoPhaseC::SetupModel(double timeLimit)
{
	// Setup Master Problem

	//Sequence dependent setup costs
	setup_pr.resize(PP.P);
	for (int p = 0; p < PP.P; ++p)
		setup_pr[p].resize(PP.P);

	for (int p = 0; p < PP.P; ++p)
		for (int r = 0; r < PP.P; ++r)
			setup_pr[p][r] = PP.Products[p].s_pr[r];

	int TotalCap = 0;
	for (int t = 0; t < PP.T; ++t)
		TotalCap += PP.K[t];

	theta = IloNumVar(env, 0, IloInfinity, "theta"); //Total estimated setup cost
	gamma = IloNumVar(env, 0, IloInfinity, "gamma"); //Total estimated setup cost of crossovers
	theta_t = CreateNumVarArray(env, PP.T, "theta_t", 0, IloInfinity); //Estimated setup cost for period t

	I = CreateNumVarArray2(env, PP.P, PP.T, "I", 0, IloInfinity); //Inventory level of product p at the end of period t
	q = CreateNumVarArray2(env, PP.P, PP.T, "q", 0, IloInfinity); //The amount of product p to be produced in period t

	x = CreateBoolVarArray2(env, PP.P, PP.T, "x"); //1 if product p is completed in macro - period t; and, 0, otherwise

	theta.setStringProperty("Type", "theta");
	theta.setIntProperty("Product", -1);

	gamma.setStringProperty("Type", "gamma");
	gamma.setIntProperty("Product", -1);

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
	model.add(theta_t);

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

	//Initialization for DP (Fill the Cache)
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

	//Set the Global Lower Bound Value
	if (!ProductPPs_In[0].empty())
		LB_theta = sDP.MinSetupMiddle(MPcache, setup_pr, -1, ProductPPs_In[0], -1);

	model.add(theta >= LB_theta);
	
	//Set the Global Lower Bound Value for all periods
	double ValIn_setup = 0;
	IloExpr ValIn(env);

	set<int> current;
	for (int t = 0; t < PP.T; ++t) {
		for (int p = 0; p < PP.P; ++p)
			if (PP.Products[p].d[t] > 0)
				current.insert(p);

		if (!current.empty()) 
			ValIn_setup = sDP.MinSetupMiddle(MPcache, setup_pr, -1, current, -1);

		ValIn += theta_t[t];

		model.add(ValIn >= ValIn_setup);
	}

	IloExpr thetaCut(env);
	for (int t = 0; t < PP.T; ++t)
		thetaCut += theta_t[t];

	model.add(theta >= thetaCut + gamma);

	//Valid Inequality (1)
	for (int t = 1; t < PP.T; ++t)
		for (int p = 0; p < PP.P; ++p)
			model.add(I[p][t - 1] >= PP.Products[p].d[t] * (1 - x[p][t]));

	//Valid Inequality (2) : Capacity releated
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

bool TwoPhaseC::Solve(double timeLimit)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	cplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));

	try
	{
		pCallback = make_unique<TwoPhaseCallback>(this);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::GlobalProgress;
		cplex.use(pCallback.get(), contextmask);

		cplex.setParam(IloCplex::IntSolLim, IntegerSolutionLimit);
		cout << "Started solving Master problem" << endl;
		Solved = cplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
	
	return Solved;
}

double TwoPhaseC::GetLB()
{
	return Solved ? cplex.getBestObjValue() : 0;
}

double TwoPhaseC::GetUB()
{
	return Solved ? cplex.getObjValue() : DBL_MAX;
}

double TwoPhaseC::GetGap()
{
	return Solved ? cplex.getMIPRelativeGap() : DBL_MAX;
}

double TwoPhaseC::GetCallbackCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CPU;
	return Total;
}

double TwoPhaseC::GetSPsolveCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += (float)Worker->TotalSPtime / 1000;
	return Total;
}

double TwoPhaseC::GetSPconsCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += (float)Worker->TotalSPcons_time / 1000;
	return Total;
}

size_t TwoPhaseC::GetCallCount()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CallCount;
	return Total;
}

size_t TwoPhaseC::GetCutCount()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CutCount;
	return Total;
}

int TwoPhaseC::GetConsts()
{
	return Solved ? cplex.getNrows() : 0;
}

int TwoPhaseC::GetVars()
{
	return Solved ? cplex.getNcols() : 0;
}

TwoPhaseCallback::TwoPhaseCallback(TwoPhaseC* pB) : pTwoPhaseC(pB)
{
	int nThreads = GetParameterValue(pTwoPhaseC->Parameters, "THREAD_COUNT");
	workers.resize(nThreads);

	for (int i = 0; i < nThreads; ++i)
		workers[i] = new WorkerWW(pTwoPhaseC);
}

TwoPhaseCallback::~TwoPhaseCallback()
{}

void TwoPhaseCallback::invoke(const IloCplex::Callback::Context& context)
{
	auto startTime = chrono::high_resolution_clock::now();

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

	// Get the right worker
	auto& worker = workers[threadNo];

	IloEnv env = context.getEnv();
	NumArray2 xVal = CreateNumArray2(env, pTwoPhaseC->PP.P, pTwoPhaseC->PP.T);
	IloNum thetaVal;

	// Get the current x solution
	switch (context.getId()) {
	case IloCplex::Callback::Context::Id::Candidate:
		if (!context.isCandidatePoint()) // The model is always bounded
			throw IloCplex::Exception(-1, "Unbounded solution");
		for (int p = 0; p < pTwoPhaseC->PP.P; ++p)
			for (int t = 0; t < pTwoPhaseC->PP.T; ++t) {
				xVal[p][t] = context.getCandidatePoint(pTwoPhaseC->x[p][t]);
			}
		thetaVal = context.getCandidatePoint(pTwoPhaseC->theta);
		break;
	case IloCplex::Callback::Context::Id::GlobalProgress:
	{
		double BestBound = context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound);
		if (abs(BestBound) >= pTwoPhaseC->BestBoundThreshold)
		{
			cout << "Aborting callback. BestBound: " << BestBound << " Threshold: " << pTwoPhaseC->BestBoundThreshold << endl;
			context.abort();
		}
		return;
	}
	default:
		// Free memory
		xVal.end();
		throw IloCplex::Exception(-1, "Unexpected contextID");
	}

	IloExprArray cuts(env, pTwoPhaseC->PP.T + 2);
	for (int t = 0; t < pTwoPhaseC->PP.T + 2; ++t)
		cuts[t] = IloExpr(env);

	double OptimalCost;

	IloBool sepStat = worker->separate(pTwoPhaseC->PP, pTwoPhaseC->Parameters, pTwoPhaseC->LB_theta, pTwoPhaseC->MPcache, thetaVal, xVal, OptimalCost, cuts);

	if (context.getId() == IloCplex::Callback::Context::Id::Candidate)
	{
		IloNum CandidateObjective = context.getCandidateObjective();
		IloNum IncumbentObjective = context.getIncumbentObjective();
		IloNum ActualCandidateObjective = CandidateObjective - thetaVal + OptimalCost;
		cout << "IncumbentObjective: " << IncumbentObjective << endl;
		cout << "ActualCandidateObjective: " << ActualCandidateObjective << endl;
	}

	xVal.end();

	if (sepStat) {
		// Add the cut
		IloRangeArray r(env, 0, cuts, IloInfinity);

		switch (context.getId()) {
		case IloCplex::Callback::Context::Id::Candidate:
			context.rejectCandidate(r);

			break;
		default:
			r.end();
			throw IloCplex::Exception(-1, "Unexpected contextID");
		}

		{
			std::unique_lock<std::mutex> lock(pTwoPhaseC->Mutex);
			for (int t = 0; t < pTwoPhaseC->PP.T + 2; ++t)
				pTwoPhaseC->GeneratedCuts.push_back(r[t]);
		}
	}
	++worker->CallCount;
	if (sepStat)
		++worker->CutCount;
	worker->CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
	worker->TotalSPtime += pTwoPhaseC->W_SPtime;
	worker->TotalSPcons_time += pTwoPhaseC->W_SPcons_time;
}

Worker::Worker(TwoPhaseC* pB) : pTwoPhaseC(pB)
{
	CallCount = 0;
	CutCount = 0;
	CPU = 0;
	TotalSPtime = 0;
	TotalSPcons_time = 0;
}

WorkerWW::WorkerWW(TwoPhaseC* pB) : Worker(pB)
{
}

bool WorkerWW::separate(ProductPeriods& PP, ParameterMap& Parameters, int& LB_theta, CacheMP& MPcache, const IloNum thetaVal, const NumArray2& xVal, double& OptimalCost, IloExprArray& cuts)
{
	int W = 0;
	for (int p = 0; p < PP.P; ++p)
		for (int t = 0; t < PP.T; ++t) {
			if (xVal[p][t] > 0.999)
				W += 1;
		}

	Matrix x_val;
	x_val.resize(PP.P);
	for (int p = 0; p < PP.P; ++p)
		x_val[p].resize(PP.T);

	for (int p = 0; p < PP.P; ++p) {
		for (int t = 0; t < PP.T; ++t) {
			if (xVal[p][t] >= 1 - 0.0001) {
				x_val[p][t] = 1;
			}
			else
				x_val[p][t] = 0;
		}
	}

		vector<int> wp;
		wp.resize(W + 2);
		vector<int> wt;
		wt.resize(W + 2);
		vector<int> wop;
		wop.resize(PP.P * PP.T - W + 2);
		vector<int> wot;
		wot.resize(PP.P * PP.T - W + 2);
		int j = 0;
		int k = 0;
		for (int p = 0; p < PP.P; ++p)
			for (int t = 0; t < PP.T; ++t)
				if (xVal[p][t] > 0.999) {
					j += 1;
					wp[j] = p;
					wt[j] = t;
				}
				else {
					k += 1;
					wop[k] = p;
					wot[k] = t;
				}
	if (GetParameterValue(Parameters, "SUBPROBLEM_TYPE") == 0) {
		auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP
		Subproblems SP(PP, Parameters);

		SP.SetupBSPModel(W, wp, wt, wop, wot);
		auto finishSP_cons = chrono::high_resolution_clock::now();
		auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
		pTwoPhaseC->W_SPcons_time = SPcons_time.count();

		cout << "Started solving Subproblem" << endl;
		auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP

		bool resultSP = SP.BSP_Solve(pTwoPhaseC->SPtimelimit);

		auto finishSP = chrono::high_resolution_clock::now();
		auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
		pTwoPhaseC->W_SPtime = SPtime.count();

		OptimalCost = SP.GetBSP_UB();
		SP.GetBSP_Solutions(W);

	}
	else if (GetParameterValue(Parameters, "SUBPROBLEM_TYPE") == 1) {
		auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP
		
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

		OptimalCost = pTwoPhaseC->sDP.MinSetupMultiPeriod(MPcache, pTwoPhaseC->setup_pr, -1, ProductPPs, -1);
		cout << "Done. Minimum setup: " << OptimalCost << endl;
		cout << "Theta: " << thetaVal << endl;

		auto finishSP = chrono::high_resolution_clock::now();
		auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
		TotalSPtime += SPtime.count();
	}
	else if (GetParameterValue(Parameters, "SUBPROBLEM_TYPE") == 2) {
		//Construct Subproblem
		auto startSP_cons = chrono::high_resolution_clock::now(); //Start clock for SP

		Subproblems CP(PP, Parameters);

		CP.SetupCPModel(W, wp, wt, wop, wot);

		auto finishSP_cons = chrono::high_resolution_clock::now();
		auto SPcons_time = chrono::duration_cast<chrono::milliseconds>(finishSP_cons - startSP_cons);
		TotalSPcons_time += SPcons_time.count();

		auto startSP = chrono::high_resolution_clock::now(); //Start clock for SP
		CP.CP_Solve(pTwoPhaseC->SPtimelimit);

		auto finishSP = chrono::high_resolution_clock::now();
		auto SPtime = chrono::duration_cast<chrono::milliseconds>(finishSP - startSP);
		TotalSPtime += SPtime.count();

		OptimalCost = CP.GetCP_Obj();
	}

	double check = 0;
	if (OptimalCost > thetaVal + 0.0001)
	{
		int LB = 0;
		
		int set_mins = 0;
		int set_min = 0;
		vector<int> set_min_p;
		set_min_p.resize(PP.P);

		for (int p = 0; p < PP.P; ++p) {
			int set_mins = INT_MAX;
			int set_min = INT_MAX;
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
					if (x_val[p][t] >= 0.5) {
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

			if (!current.empty()) 
				LBsetupDP = pTwoPhaseC->sDP.MinSetupMiddle(MPcache, pTwoPhaseC->setup_pr, -1, current, -1);
			else
				LBsetupDP = 0;

			LBsetup[t] = LBsetupDP;

			auto currentCopy = current;

			if (!current.empty()) {
				for (int p = 0; p < PP.P; ++p) {
					if (x_val[p][t] >= 0.5)
						current.erase(p);

					if (!current.empty())
						LBsetupDP = pTwoPhaseC->sDP.MinSetupMiddle(MPcache, pTwoPhaseC->setup_pr, -1, current, -1);
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
		
		IloExpr Sum1(pTwoPhaseC->env);
		IloExpr Sum2(pTwoPhaseC->env);

		for (int p = 0; p < PP.P; ++p)
			for (int t = 0; t < PP.T; ++t)
				if (x_val[p][t] >= 0.5)
					Sum1 += pTwoPhaseC->x[p][t];
				else
					Sum2 += pTwoPhaseC->x[p][t];

		IloExpr cutLhs(pTwoPhaseC->env);
		cutLhs -= (OptimalCost - LB) * (Sum1)-(OptimalCost - LB) * (W - 1) + LB - pTwoPhaseC->theta;
		
		
		IloExprArray cutDPs(pTwoPhaseC->env, pTwoPhaseC->PP.T);
		for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
			cutDPs[t] = IloExpr(pTwoPhaseC->env);
		
		IloExprArray DPCut(pTwoPhaseC->env, pTwoPhaseC->PP.T);
		for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
			DPCut[t] = IloExpr(pTwoPhaseC->env);

		for (int t = 0; t < pTwoPhaseC->PP.T; ++t) {
			for (int p = 0; p < pTwoPhaseC->PP.P; ++p) {
				if (x_val[p][t] >= 0.5)
					DPCut[t] += (LBsetup[t] * (1 - pTwoPhaseC->x[p][t]));
			}
			cutDPs[t] -= LBsetup[t] - DPCut[t] - pTwoPhaseC->theta_t[t];
		}

		double LBsetuptotal = 0;
		for (int t = 0; t < PP.T; ++t)
			LBsetuptotal += LBsetup[t];

		double LBresidual = OptimalCost - LBsetuptotal;

		IloExpr cutGMs(pTwoPhaseC->env);
		IloExpr GMCut(pTwoPhaseC->env);

		for (int t = 0; t < PP.T; ++t) {
			for (int p = 0; p < PP.P; ++p) {
				if (x_val[p][t] >= 0.5)
					GMCut += LBresidual * (1 - pTwoPhaseC->x[p][t]);
				else
					GMCut += LBresidual * (pTwoPhaseC->x[p][t]);
			}
		}

		cutGMs -= LBresidual - GMCut - pTwoPhaseC->gamma;

		cuts[0] = cutLhs;
		cuts[1] = cutGMs;
		for (int t = 2; t < PP.T + 2; ++t)
			cuts[t] = cutDPs[t-2];

		return true;
	}
	
	//SP.SPmodel.end();
	return false;
}