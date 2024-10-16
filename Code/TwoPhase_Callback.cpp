
#include "TwoPhase_Callback.h"
#include "DataPrep.h"
#include "setupDP.h"
#include "Subproblems.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;
typedef vector<IloExpr> Expr_vec;

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
	//MPmodel.add(theta >= LB_theta);

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
	return Solved ? cplex.getObjValue() : 0;
}

double TwoPhaseC::GetUB()
{
	return Solved ? cplex.getBestObjValue() : DBL_MAX;
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
	IloExpr Sum1(env);
	IloExpr Sum2(env);
	IloExpr GMCut(env);
	double OptimalCost;

	Expr_vec cutDPs;
	cutDPs.resize(pTwoPhaseC->PP.T);

	Expr_vec DPCut;
	DPCut.resize(pTwoPhaseC->PP.T);


	IloBool sepStat = worker->separate(pTwoPhaseC->PP, pTwoPhaseC->Parameters, pTwoPhaseC->LB_theta, thetaVal, xVal, OptimalCost, cutLhs, cutDPs, cutGMs, Sum1, Sum2, DPCut, GMCut);

	if (context.getId() == IloCplex::Callback::Context::Id::Candidate)
	{
		IloNum CandidateObjective = context.getCandidateObjective();
		IloNum IncumbentObjective = context.getIncumbentObjective();
		IloNum ActualCandidateObjective = CandidateObjective + thetaVal - OptimalCost;
	}

	xVal.end();

	if (sepStat) {
		// Add the cut
		IloRange r(env, 0, cutLhs, IloInfinity);
		//IloRange rr(env, 0, cutGMs, IloInfinity);

		switch (context.getId()) {
		case IloCplex::Callback::Context::Id::Candidate:
			context.rejectCandidate(r);
			//context.rejectCandidate(rr);

			break;
		default:
			r.end();
			//rr.end();
			throw IloCplex::Exception(-1, "Unexpected contextID");
		}

		{
			std::unique_lock<std::mutex> lock(pTwoPhaseC->Mutex);
			pTwoPhaseC->GeneratedCuts.push_back(r);
			//pTwoPhaseC->GeneratedCuts.push_back(rr);
		}
		/*
		for (int t = 0; t < pTwoPhaseC->PP.T; ++t){
			IloRange rrr(env, 0, cutDPs[t], IloInfinity);

			switch (context.getId()) {
			case IloCplex::Callback::Context::Id::Candidate:
				context.rejectCandidate(rrr);

				break;
			default:
				rrr.end();
				throw IloCplex::Exception(-1, "Unexpected contextID");
			}

			{
				std::unique_lock<std::mutex> lock(pTwoPhaseC->Mutex);
				pTwoPhaseC->GeneratedCuts.push_back(rrr);
			}
		}*/
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

bool WorkerWW::separate(ProductPeriods& PP, ParameterMap& Parameters, int& LB_theta, const IloNum thetaVal, const NumArray2& xVal, double& OptimalCost, IloExpr& cutLhs, Expr_vec& cutDPs, IloExpr& cutGMs, IloExpr& Sum1, IloExpr& Sum2, Expr_vec& DPCut, IloExpr& GMCut)
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
		//for (int t = 0; t < PP.T; ++t)
			//cutDPs[t].clear();

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
		
		for (int j = 1; j < W + 1; ++j)
			Sum1 += pTwoPhaseC->x[wp[j]][wt[j]];

		for (int k = 1; k < (PP.P * PP.T - W + 1); ++k)
			Sum2 += pTwoPhaseC->x[wop[k]][wot[k]];

		cutLhs -= (OptimalCost - LB_theta) * (Sum1)-(OptimalCost - LB_theta) * (W - 1) + LB_theta - pTwoPhaseC->theta;
		
		for (int t = 0; t < PP.T; ++t) {
			for (int p = 0; p < PP.P; ++p) {
				if (x_val[p][t] >= 0.5)
					LBsetup[t];//DPCut[t] += (LBsetup[t] - LBsetupPermin[t]) * (1 - pTwoPhaseC->x[p][t]);
			}
			//cutDPs[t] -= (pTwoPhaseC->theta_t[t] >= LBsetup[t] - DPCut[t]);
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

		//cutGMs -= (pTwoPhaseC->gamma >= LBresidual - GMCut);

		return true;
	}
	SP.SPmodel.end();
	return false;
}