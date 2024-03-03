#include "TwoPhase.h"
#include "DataPrep.h"
#include <chrono>
#include <numeric>

typedef vector<vector<double>> Matrix;

TwoPhase::TwoPhase(ProductPeriods& PPIn) : PP(PPIn), model(env), cplex(env)
{
	BestBoundThreshold = DBL_MAX;
}

TwoPhase::~TwoPhase()
{
	env.end();
}

void TwoPhase::SetupModel()
{
	int P = PP.P;
	int T = PP.T;
	int S = PP.S;
	int M = PP.M;

	I = CreateNumVarArray2(env, P, T, "I", 0, IloInfinity);
	q = CreateNumVarArray2(env, P, S, "q", 0, IloInfinity);

	y = CreateBoolVarArray2(env, P, S, "y");
	z = CreateBoolVarArray3(env, P, P, S, "z");


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
				model.add(z[p][r][s] >= y[p][s - 1] + y[r][s] - 1);

	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			model.add(z[p][r][0] == 0);

	// Constraint (7.1)
	/*for (int p = 0; p < P; ++p)
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
		}*/

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

bool TwoPhase::Solve(double timeLimit)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	//cplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));

	try
	{
		pCallback = make_unique<TwoPhaseCallback>(this);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::GlobalProgress;
		cplex.use(pCallback.get(), contextmask);

		/* Initial solution
		if (GetParameterValue(Parameters, "USE_INITIAL_SOLUTION") && !InitialSolution.empty())
		{
			IloNumVarArray vars(env);
			IloNumArray values(env);
			for (int m = 0; m < MS.M; ++m)
			{
				vars.add(z[m]);
				values.add(InitialSolution[m]);
			}

			if (cplex.getNMIPStarts())
				cplex.deleteMIPStarts(0);
			cplex.addMIPStart(vars, values, IloCplex::MIPStartEffort::MIPStartSolveFixed);

			vars.end();
			values.end();
		}

		cplex.setParam(IloCplex::IntSolLim, IntegerSolutionLimit);*/
		Solved = cplex.solve();
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

double TwoPhase::GetCallbackCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CPU;
	return Total;
}
/*
TwoPhaseCallback::TwoPhaseCallback(TwoPhase* pB) : pTwoPhase(pB)
{
	int nThreads = GetParameterValue(pTwoPhase->Parameters, "THREAD_COUNT");
	workers.resize(nThreads);

	for (int i = 0; i < nThreads; ++i)
		workers[i] = new WorkerWW(pTwoPhase);
}

TwoPhaseCallback::~DTwoPhaseCallback()
{}

void TwoPhaseCallback::invoke(const IloCplex::Callback::Context& context)
{
	auto startTime = chrono::high_resolution_clock::now();

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

	// Get the right worker
	auto& worker = workers[threadNo];

	IloEnv env = context.getEnv();
	IloNumArray zVal(env, pTwoPhase->PP.P);
	IloNum thetaVal;

	// Get the current z solution
	switch (context.getId()) {
	case IloCplex::Callback::Context::Id::Candidate:
		if (!context.isCandidatePoint()) // The model is always bounded
			throw IloCplex::Exception(-1, "Unbounded solution");
		context.getCandidatePoint(pTwoPhase->z, zVal);
		thetaVal = context.getCandidatePoint(pTwoPhase->Cost);
		break;
	case IloCplex::Callback::Context::Id::Relaxation:
		context.getRelaxationPoint(pTwoPhase->z, zVal);
		thetaVal = context.getRelaxationPoint(pDecomposition->Cost);
		break;
	case IloCplex::Callback::Context::Id::GlobalProgress:
	{
		double BestBound = context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound);
		if (abs(BestBound) >= pDecomposition->BestBoundThreshold)
		{
			cout << "Aborting callback. BestBound: " << BestBound << " Threshold: " << pDecomposition->BestBoundThreshold << endl;
			context.abort();
		}
		return;
	}
	default:
		// Free memory
		zVal.end();
		throw IloCplex::Exception(-1, "Unexpected contextID");
	}

	IloExpr cutLhs(env);
	double OptimalCost;
	bool applyLifting = context.isCandidatePoint();
	IloBool sepStat = worker->separate(pTwoPhase->PP, thetaVal, zVal, OptimalCost, cutLhs, applyLifting);

	if (context.getId() == IloCplex::Callback::Context::Id::Candidate)
	{
		IloNum CandidateObjective = context.getCandidateObjective();
		IloNum IncumbentObjective = context.getIncumbentObjective();
		IloNum ActualCandidateObjective = CandidateObjective + thetaVal - OptimalCost;
	}

	zVal.end();

	if (sepStat) {
		// Add the cut
		IloRange r(env, 0, cutLhs, IloInfinity);

		switch (context.getId()) {
		case IloCplex::Callback::Context::Id::Candidate:
			context.rejectCandidate(r);
			break;
		case IloCplex::Callback::Context::Id::Relaxation:
			if (OptimalCost - thetaVal > 1)
			{
				cout << r << endl;
				context.addUserCut(r,
					IloCplex::UseCutPurge,
					IloFalse);
			}
			break;
		default:
			r.end();
			throw IloCplex::Exception(-1, "Unexpected contextID");
		}

		{
			unique_lock lock(pTwoPhase->Mutex);
			pTwoPhase->GeneratedCuts.push_back(r);
		}
	}
	++worker->CallCount;
	if (sepStat)
		++worker->CutCount;
	worker->CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
}

Worker::Worker(TwoPhase* pB) : pTwoPhase(pB)
{
	CallCount = 0;
	CutCount = 0;
	CPU = 0;
}

WorkerWW::WorkerWW(TwoPhase* pB) : Worker(pB)
{
}

bool WorkerWW::separate(ProductPeriods& PP, const IloNum thetaVal, const IloNumArray& zVal, double& OptimalCost, IloExpr& cutLhs, bool applyLifting)
{
	//int nSelected = 0;
	DoubleVector Demand(PP.T, 0);
	double TotalIndividualCost = 0;
	double TotalConstantH = 0;
	vector<bool> IsSelected(PP.P);
	for (int m = 0; m < PP.P; ++m)
	{
		auto& Market = PP.Markets[m];

		for (int t = 0; t < PP.T; ++t)
			Demand[t] += Market.Demand[t] * zVal[m];
	}

	IntegerVector Ones(PP.T);
	DoubleVector Duals(PP.T, 0);
	double ConstantH;

	int counter;

	PP.WWgeneralPD(PP.T, Demand.data(), PP.SetupCost.data(), PP.ProductionCost.data(), PP.HoldingCost.data(), OptimalCost, Ones.data(), Duals.data(), ConstantH, counter);

	for (auto d : Duals)
	{
		if (isnan(d) || isnan(ConstantH))
		{
			cout << "Numerical trouble detected" << endl;
			PP.WWgeneralPD(MS.T, Demand.data(), PP.SetupCost.data(), PP.ProductionCost.data(), PP.HoldingCost.data(), OptimalCost, Ones.data(), Duals.data(), ConstantH, counter);
			return false;
		}
	}

	double check = 0;
	if (OptimalCost > thetaVal + pTwoPhase->Epsilon)
	{
		cutLhs.clear();
		cutLhs += pTwoPhase->Cost;

		for (int m = 0; m < PP.P; ++m)
		{
			auto& Market = PP.Markets[m];
			if (zVal[m] > pTwoPhase->Epsilon)
			{
				double Coefficient = 0;
				for (int t = 0; t < PP.T; ++t)
					Coefficient += Duals[t] * Market.Demand[t] * zVal[m];

				cutLhs -= pTwoPhase->z[m] * Coefficient;
				check += (zVal[m] * Coefficient);
			}
			else if (applyLifting)
			{
				// Apply lifting
				cutLhs -= pTwoPhase->z[m] * (Market.VariableCost);
			}
		}

		return true;
	}
	return false;
}


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
	*/
