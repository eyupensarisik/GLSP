
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

void TwoPhaseC::SetupModel()
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
	for (int t = 0; t < PP.T; ++t) {
		for (int p = 0; p < PP.P; ++p)
			if (PP.Products[p].d[t] > 0)
				current.insert(p);

		if (!current.empty()) {
			//cout << "Finding minimum setup" << endl;
			ValIn_setup = sDP.MinSetup(Cache, setup_pr, -1, current);
			//cout << "Done. Minimum setup: " << ValIn_setup << endl;
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

double TwoPhaseC::GetCallbackCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CPU;
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

	IloExpr cutLhs(env);
	IloExpr cutGMs(env);
	IloExprArray cutDPs(env, pTwoPhaseC->PP.T);
	for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
		cutDPs[t] = IloExpr(env);

	IloExpr Sum1(env);
	IloExpr Sum2(env);
	IloExpr GMCut(env);
	IloExprArray DPCut(env, pTwoPhaseC->PP.T);
	for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
		DPCut[t] = IloExpr(env);

	double OptimalCost;


	IloBool sepStat = worker->separate(pTwoPhaseC->PP, pTwoPhaseC->Parameters, pTwoPhaseC->LB_theta, thetaVal, xVal, OptimalCost, cutLhs, cutGMs, Sum1, Sum2, GMCut, DPCut, cutDPs);

	if (context.getId() == IloCplex::Callback::Context::Id::Candidate)
	{
		IloNum CandidateObjective = context.getCandidateObjective();
		IloNum IncumbentObjective = context.getIncumbentObjective();
		IloNum ActualCandidateObjective = CandidateObjective + thetaVal - OptimalCost;
	}

	xVal.end();

	if (sepStat) {
		// Add the cut
		IloRange Lshaped_r(env, 0, cutLhs, IloInfinity);
		IloRange GM_r(env, 0, cutGMs, IloInfinity);
		IloRangeArray DP_r(env, 0, cutDPs, IloInfinity);

		switch (context.getId()) {
		case IloCplex::Callback::Context::Id::Candidate:
			context.rejectCandidate(Lshaped_r);
			context.rejectCandidate(GM_r);
			context.rejectCandidate(DP_r);

			break;
		default:
			Lshaped_r.end();
			GM_r.end();
			DP_r.end();
			throw IloCplex::Exception(-1, "Unexpected contextID");
		}

		{
			std::unique_lock<std::mutex> lock(pTwoPhaseC->Mutex);
			pTwoPhaseC->GeneratedCuts.push_back(Lshaped_r);
			pTwoPhaseC->GeneratedCuts.push_back(GM_r);
			for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
				pTwoPhaseC->GeneratedCuts.push_back(DP_r[t]);
		}
	}
	++worker->CallCount;
	if (sepStat)
		++worker->CutCount;
	worker->CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
}

Worker::Worker(TwoPhaseC* pB) : pTwoPhaseC(pB)
{
	CallCount = 0;
	CutCount = 0;
	CPU = 0;
}

WorkerWW::WorkerWW(TwoPhaseC* pB) : Worker(pB)
{
}

bool WorkerWW::separate(ProductPeriods& PP, ParameterMap& Parameters, int& LB_theta, const IloNum thetaVal, const NumArray2& xVal, double& OptimalCost, IloExpr& cutLhs, IloExpr& cutGMs, IloExpr& Sum1, IloExpr& Sum2, IloExpr& GMCut, IloExprArray& DPCut, IloExprArray& cutDPs)
{
	int W = 0;
	for (int p = 0; p < pTwoPhaseC->PP.P; ++p)
		for (int t = 0; t < pTwoPhaseC->PP.T; ++t) {
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
	wop.resize(pTwoPhaseC->PP.P * pTwoPhaseC->PP.T - W + 2);
	vector<int> wot;
	wot.resize(pTwoPhaseC->PP.P * pTwoPhaseC->PP.T - W + 2);
	int j = 0;
	int k = 0;
	for (int p = 0; p < pTwoPhaseC->PP.P; ++p)
		for (int t = 0; t < pTwoPhaseC->PP.T; ++t)
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

	Subproblems SP(PP, Parameters);

	SP.SetupBSPModel(W, wp, wt, wop, wot);

	bool resultSP = SP.BSP_Solve();

	OptimalCost = SP.GetBSP_UB();
	SP.GetBSP_Solutions(W);

	double check = 0;
	if (OptimalCost > thetaVal + 0.0001)
	{
		cutLhs.clear();
		cutGMs.clear();
		//cutDPs.clear();

		vector<int> SP_Sol;
		SP_Sol.resize(W);

		vector<int> SP_Sol2;

		int j = 0;
		int i = 0;
		for (int l = 1; l < W + 1; ++l)
			if (SP.e_val[j][l] > 0.5) {
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
				//cout << "Finding minimum setup" << endl;
				LBsetupDP = sDP.MinSetup(Cache, SP.SP_setup, -1, current);
				//cout << "Done. Minimum setup: " << LBsetupDP << endl;
			}
			else
				LBsetupDP = 0;

			LBsetup[t] = LBsetupDP;

			auto currentCopy = current;
			if (pre_per != per) {
				for (int i = pre_per; i < per; ++i) {
					current.erase(SP_Sol[i]);

					if (!current.empty()) {
						//cout << "Finding minimum setup" << endl;
						LBsetupDP = sDP.MinSetup(Cache, SP.SP_setup, -1, current);
						//cout << "Done. Minimum setup: " << LBsetupDP << endl;
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
		
		for (int j = 1; j < W + 1; ++j)
			Sum1 += pTwoPhaseC->x[wp[j]][wt[j]];

		for (int k = 1; k < (PP.P * PP.T - W + 1); ++k)
			Sum2 += pTwoPhaseC->x[wop[k]][wot[k]];

		cutLhs -= (OptimalCost - LB) * (Sum1)-(OptimalCost - LB) * (W - 1) + LB - pTwoPhaseC->theta;
		
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

		for (int t = 0; t < PP.T; ++t) {
			for (int p = 0; p < PP.P; ++p) {
				if (x_val[p][t] >= 0.5)
					GMCut += LBresidual * (1 - pTwoPhaseC->x[p][t]);
				else
					GMCut += LBresidual * (pTwoPhaseC->x[p][t]);
			}
		}

		cutGMs -= LBresidual - GMCut - pTwoPhaseC->gamma;
		return true;
	}
	SP.SPmodel.end();
	return false;
}