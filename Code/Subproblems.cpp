
#include "Subproblems.h"
#include "DataPrep.h"
#include <chrono>
#include <numeric>


typedef vector<vector<double>> Matrix;

Subproblems::Subproblems(ProductPeriods& PPIn, ParameterMap& PM) : PP(PPIn), Parameters(PM), SPmodel(SPenv), SPcplex(SPenv), CPmodel(CPenv), CPcplex(CPenv), TSPmodel(TSPenv), TSPcplex(TSPenv), TSP_CPmodel(TSP_CPenv), TSP_CPcplex(TSP_CPenv)
{
	BestBoundThreshold = DBL_MAX;
}

Subproblems::~Subproblems()
{
	TSPenv.end();
	SPenv.end();
	CPenv.end();
	TSP_CPenv.end();
}

void Subproblems::SetupBSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot)
{
	// Setup Integer Subroblem (IP)
	SPcplex.setOut(SPenv.getNullStream());
	//Find due-dates and release dates
	vector<int> dd;
	dd.resize(PP.T);
	vector<int> rd;
	rd.resize(PP.T);
	int d_sum = 0;
	for (int t = 0; t < PP.T; ++t) {
		rd[t] = d_sum;
		d_sum += PP.K[t];
		dd[t] = d_sum;
	}

	//Calculate the Big M value
	int TotalCapacity = 0;
	for (int t = 0; t < PP.T; ++t)
		TotalCapacity += PP.K[t];

	int M = TotalCapacity;

	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_a[j] = PP.Products[wp[j]].a;
	SP_a[0] = 0;
	SP_a[W + 1] = 0;

	vector<int> SP_b;
	SP_b.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_b[j] = dd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = d_sum;

	vector<int> SP_r;
	SP_r.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_r[j] = rd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = 0;

	SP_setup.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		SP_setup[j].resize(W + 1);

	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SP_setup[j][l] = PP.Products[wp[j]].s_pr[wp[l]];
		}

	C = CreateNumVarArray(SPenv, W + 2, "C", 0, IloInfinity); //The completion time of job j
	e = CreateBoolVarArray2(SPenv, W + 2, W + 2, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise

	IloExpr SPobj1(SPenv);
	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (j != l)
				SPobj1 += SP_setup[j][l] * e[j][l];

	SPobj = IloMinimize(SPenv, SPobj1);
	SPmodel.add(SPobj);

	SPcplex = IloCplex(SPmodel);

	// Constraint (2)
	for (int j = 0; j < W + 1; ++j) {
		IloExpr Each1(SPenv);
		for (int l = 1; l < W + 2; ++l)
			if (j != l)
				Each1 += e[j][l];
		SPmodel.add(Each1 == 1);
	}

	// Constraint (3)
	for (int j = 1; j < W + 2; ++j) {
		IloExpr Each2(SPenv);
		for (int l = 0; l < W + 1; ++l)
			if (j != l)
				Each2 += e[l][j];
		SPmodel.add(Each2 == 1);
	}

	// Constraint (4)
	for (int j = 0; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SPmodel.add(C[j] + SP_a[l] <= C[l] + M * (1 - e[j][l]));
		}

	// Constraint (5)
	for (int j = 0; j < W + 1; ++j) {
		SPmodel.add(C[j] <= SP_b[j]);
	}
	
	// Constraint (6)
	for (int j = 0; j < W + 1; ++j) {
		SPmodel.add(C[j] >= SP_a[j] + SP_r[j]);
	}

}

bool Subproblems::BSP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		SPcplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		SPcplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		SPcplex.exportModel("IP.lp");
		Solved = SPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetBSP_UB()
{
	return Solved ? SPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetBSP_LB()
{
	return Solved ? SPcplex.getBestObjValue() : DBL_MAX;
}

Matrix Subproblems::GetBSP_Solutions(int W)
{
	e_val.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		e_val[j].resize(W + 1);

	for (int j = 0; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (SPcplex.getValue(e[j][l]) > 0.5)
				e_val[j][l] = SPcplex.getValue(e[j][l]);

	return e_val;
}

void Subproblems::SetupTSPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot)
{
	// Setup TSP Subroblem
	TSPcplex.setOut(TSPenv.getNullStream());
	//Find due-dates and release dates
	vector<int> dd;
	dd.resize(PP.T);
	vector<int> rd;
	rd.resize(PP.T);
	int d_sum = 0;
	for (int t = 0; t < PP.T; ++t) {
		rd[t] = d_sum;
		d_sum += PP.K[t];
		dd[t] = d_sum;
	}

	//Calculate the Big M value
	int TotalCapacity = 0;
	for (int t = 0; t < PP.T; ++t)
		TotalCapacity += PP.K[t];

	int M = TotalCapacity;

	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_a[j] = PP.Products[wp[j]].a;
	SP_a[0] = 0;
	SP_a[W + 1] = 0;

	vector<int> SP_b;
	SP_b.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_b[j] = dd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = d_sum;

	vector<int> SP_r;
	SP_r.resize(W + 2);
	for (int j = 1; j < W + 1; ++j)
		SP_r[j] = rd[wt[j]];
	SP_b[0] = 0;
	SP_b[W + 1] = 0;

	SP_setup.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		SP_setup[j].resize(W + 1);

	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l) {
			SP_setup[j][l] = PP.Products[wp[j]].s_pr[wp[l]];
		}

	vector<pair<int, int>> P;
	
	for (int j = 1; j < W + 1; ++j)
		if (wt[j] == 0)
			P.push_back(make_pair(0, j));

	for (int t = 0; t < PP.T - 1; ++t)
		for (int j = 1; j < W + 1; ++j)
			if (wt[j] == t)
				for (int l = 1; l < W + 1; ++l)
					if (wt[l] == t + 1)
						if (j != l)
							P.push_back(make_pair(j, l));

	for (int j = 1; j < W + 1; ++j)
		if (wt[j] == PP.T - 1)
			P.push_back(make_pair(j, W + 1));
	
	// u[i] for i in V, u[0] = 1 (node 1), u[i] in [2,n] for i != 0
	IloNumVarArray u(TSPenv, W + 2);
	{
		// u_0 = 0  (constraint 6)
		u[0] = IloNumVar(TSPenv, 0.0, 0.0, ILOFLOAT);
		char name[32];
		snprintf(name, sizeof(name), "u(0)");
		u[0].setName(name);
	}
	for (int j = 1; j < W + 2; j++) {
		char name[32];
		snprintf(name, sizeof(name), "u(%d)", j);
		// 2 <= u_i <= n  (constraint 7)
		u[j] = IloNumVar(TSPenv, 1.0, (double)(W + 1), ILOFLOAT, name);
	}
	/**/
	e = CreateBoolVarArray2(TSPenv, W + 2, W + 2, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise
	eta = CreateNumVarArray2(TSPenv, W + 2, W + 2, "eta", 0, IloInfinity);

	IloExpr TSPobj1(TSPenv);
	for (int j = 1; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (j != l)
				TSPobj1 += SP_setup[j][l] * e[j][l];

	TSPobj = IloMinimize(TSPenv, TSPobj1);
	TSPmodel.add(TSPobj);

	TSPcplex = IloCplex(TSPmodel);

	// Constraint (2)
	for (int j = 0; j < W + 1; ++j) {
		IloExpr Each1(TSPenv);
		for (int l = 1; l < W + 2; ++l)
			if (j != l)
				Each1 += e[j][l];
		TSPmodel.add(Each1 == 1);
	}

	// Constraint (3)
	for (int j = 1; j < W + 2; ++j) {
		IloExpr Each2(TSPenv);
		for (int l = 0; l < W + 1; ++l)
			if (j != l)
				Each2 += e[l][j];
		TSPmodel.add(Each2 == 1);
	}

	// MTZ Formulation
	
	for (int j = 1; j < W + 1; ++j) {
		for (int l = 1; l < W + 2; ++l)
			if (j != l)
				TSPmodel.add(u[j] - u[l] + (W + 1) * e[j][l] <= W);
	}

	for (auto& [pi, pj] : P) {
		int j = pi; // convert to 0-indexed
		int l = pj;
		IloExpr prec(TSPenv);
		prec += u[j] - u[l];
		TSPmodel.add(prec <= -1.0);
		prec.end();
	}
	
	//Single Commodity Flow Formulation
	/*
	for (int j = 0; j < W + 2; ++j) {
		for (int l = 0; l < W + 2; ++l)
			TSPmodel.add(eta[j][l] <= (W + 1) * e[j][l]);
	}

	IloExpr cExpr7(TSPenv);
	for (int l = 1; l < W + 1; ++l) {
		cExpr7 += eta[0][l];
	}
	TSPmodel.add(cExpr7 == (W + 1));

	for (int j = 1; j < W + 1; ++j) {
		IloExpr cExpr9(TSPenv);
		IloExpr cExpr10(TSPenv);

		for (int l = 0; l < W + 1; ++l)
			cExpr9 += eta[l][j];
		for (int l = 1; l < W + 2; ++l)
			cExpr10 += eta[j][l];

		TSPmodel.add(cExpr9 - cExpr10 == 1);
	}

	for (auto& [pi, pj] : P) {
		int j = pi; // convert to 0-indexed
		int l = pj;
		
		IloExpr cExpr11(TSPenv);
		IloExpr cExpr12(TSPenv);

		for (int k = 0; k < W + 1; ++k)
			if (k != j)
				cExpr11 += eta[k][j];
		for (int m = 0; m < W + 1; ++m)
			if (m != l)
				cExpr12 += eta[m][l];

		TSPmodel.add(cExpr11 - cExpr12 >= 1.0);
	}*/
}

bool Subproblems::TSP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		TSPcplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		TSPcplex.setParam(IloCplex::TiLim, timeLimit);

	try
	{
		//TSPcplex.exportModel("TSP.lp");
		Solved = TSPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetTSP_UB()
{
	return Solved ? TSPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetTSP_LB()
{
	return Solved ? TSPcplex.getBestObjValue() : DBL_MAX;
}

Matrix Subproblems::GetTSP_Solutions(int W)
{
	e_val.resize(W + 1);
	for (int j = 0; j < W + 1; ++j)
		e_val[j].resize(W + 1);

	for (int j = 0; j < W + 1; ++j)
		for (int l = 1; l < W + 1; ++l)
			if (j != l)
				if (TSPcplex.getValue(e[j][l]) > 0.5)
					e_val[j][l] = TSPcplex.getValue(e[j][l]);


	return e_val;
}

void Subproblems::SetupCPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot, Matrix x_val)
{
	// Setup Constraint Programming Subroblem (CP)

	//Find due-dates and release dates
	vector<int> dd;
	dd.resize(PP.T);
	vector<int> rd;
	rd.resize(PP.T);
	int d_sum = 0;
	for (int t = 0; t < PP.T; ++t){
		rd[t] = d_sum;
		for (int p = 0; p < PP.P; ++p) {
			d_sum += x_val[p][t];
		}
		dd[t] = d_sum;
	}

	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W);
	for (int j = 0; j < W; ++j)
		SP_a[j] = PP.Products[wp[j+1]].a;

	vector<int> SP_b;
	SP_b.resize(W);
	for (int j = 0; j < W; ++j)
		SP_b[j] = dd[wt[j+1]];

	vector<int> SP_r;
	SP_r.resize(W);
	for (int j = 0; j < W; ++j)
		SP_r[j] = rd[wt[j+1]];

	SP_setup.resize(W);
	for (int j = 0; j < W; ++j)
		SP_setup[j].resize(W);

	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			SP_setup[j][l] = PP.Products[wp[j+1]].s_pr[wp[l+1]];
		}

	IntArray2 iloSetup = CreateIntArray2(CPenv, W, W);
	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			iloSetup[j][l] = SP_setup[j][l];
		}

	IloIntArray type(CPenv, W);
	for (int j = 0; j < W; ++j)
		type[j] = j;

	IloIntervalVarArray Y(CPenv, W);
	for (int j = 0; j < W; j++)
	{
		Y[j] = IloIntervalVar(CPenv, SP_a[j], "Y");
		Y[j].setPresent();
	}

	IloIntervalSequenceVar V (CPenv, Y, type, "V");

	IloExpr CPobj1(CPenv);
	for (int j = 0; j < W; j++){
		IloIntExpr l = IloTypeOfNext(V, Y[j], j);
		CPobj1 += IloElement(iloSetup[j], l);
	}

	CPobj = IloMinimize(CPenv, CPobj1);
	CPmodel.add(CPobj);

	CPcplex = IloCP(CPmodel);

	CPmodel.add(IloNoOverlap(CPenv, V));

	for (int j = 0; j < W; j++)
		CPmodel.add(IloStartOf(Y[j], 0) >= SP_r[j]);

	for (int j = 0; j < W; j++)
		CPmodel.add(IloEndOf(Y[j], 0) <= SP_b[j]);
}

bool Subproblems::CP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		CPcplex.setParameter(IloCP::Workers, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		CPcplex.setParameter(IloCP::TimeLimit, timeLimit);

	try
	{
		CPcplex.exportModel("CP_GLSP.cpo");
		Solved = CPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetCP_Obj()
{
	return Solved ? CPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetCP_Bound()
{
	return Solved ? CPcplex.getObjBound() : 0;
}

double Subproblems::GetCP_Gap()
{
	return Solved ? CPcplex.getObjGap() : 100;
}

///////////////////////

void Subproblems::SetupTSP_CPModel(int W, vector<int> wp, vector<int> wt, vector<int> wop, vector<int> wot, Matrix x_val)
{
	// Setup TSP Constraint Programming Subroblem (CP)
	vector<pair<int, int>> P;
	/*
	for (int j = 1; j < W + 1; ++j)
		if (wt[j] == 0)
			P.push_back(make_pair(0, j));
*/
	for (int t = 0; t < PP.T - 1; ++t)
		for (int j = 1; j < W + 1; ++j)
			if (wt[j] == t)
				for (int l = 1; l < W + 1; ++l)
					if (wt[l] == t + 1)
						if (j != l)
							P.push_back(make_pair(j-1, l-1));
	/*
	for (int j = 1; j < W + 1; ++j)
		if (wt[j] == PP.T - 1)
			P.push_back(make_pair(j, W + 1));
	*/
	//Construct Subproblem

	vector<int> SP_a;
	SP_a.resize(W);
	for (int j = 0; j < W; ++j)
		SP_a[j] = PP.Products[wp[j + 1]].a;

	SP_setup.resize(W);
	for (int j = 0; j < W; ++j)
		SP_setup[j].resize(W);

	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			SP_setup[j][l] = PP.Products[wp[j + 1]].s_pr[wp[l + 1]];
		}

	IntArray2 iloSetup = CreateIntArray2(TSP_CPenv, W, W);
	for (int j = 0; j < W; ++j)
		for (int l = 0; l < W; ++l) {
			iloSetup[j][l] = SP_setup[j][l];
		}

	IloIntArray type(TSP_CPenv, W);
	for (int j = 0; j < W; ++j)
		type[j] = j;

	IloIntervalVarArray Y(TSP_CPenv, W);
	for (int j = 0; j < W; j++)
	{
		Y[j] = IloIntervalVar(TSP_CPenv, SP_a[j], "Y");
		Y[j].setPresent();
	}

	IloIntervalSequenceVar V(TSP_CPenv, Y, type, "V");
	IloExpr CPobj1(TSP_CPenv);
	for (int j = 0; j < W; j++) {
		IloIntExpr l = IloTypeOfNext(V, Y[j], j);
		CPobj1 += IloElement(iloSetup[j], l);
	}

	TSP_CPobj = IloMinimize(TSP_CPenv, CPobj1);
	TSP_CPmodel.add(TSP_CPobj);

	TSP_CPcplex = IloCP(TSP_CPmodel);

	TSP_CPmodel.add(IloNoOverlap(TSP_CPenv, V));

	for (auto& [pi, pj] : P) {
		int j = pi;
		int l = pj;
		TSP_CPmodel.add(IloEndBeforeStart(TSP_CPenv, Y[j], Y[l]));
	}

	for (int j = 0; j < W; j++) {
		TSP_CPmodel.add(IloEndOf(Y[j]) <= W);
	}
}

bool Subproblems::TSP_CP_Solve(double timeLimit)
{
	if (GetParameterValue(Parameters, "THREAD_COUNT"))
		TSP_CPcplex.setParameter(IloCP::Workers, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (timeLimit > 0)
		TSP_CPcplex.setParameter(IloCP::TimeLimit, timeLimit);

	try
	{
		TSP_CPcplex.exportModel("CP_GLSP.cpo");
		Solved = TSP_CPcplex.solve();
	}
	catch (IloException& ex)
	{
		cout << ex.getMessage() << endl;
	}

	return Solved;
}

double Subproblems::GetTSP_CP_Obj()
{
	return Solved ? TSP_CPcplex.getObjValue() : DBL_MAX;
}

double Subproblems::GetTSP_CP_Bound()
{
	return Solved ? TSP_CPcplex.getObjBound() : 0;
}

double Subproblems::GetTSP_CP_Gap()
{
	return Solved ? TSP_CPcplex.getObjGap() : 100;
}
/*
void SolveTSPWithPrecedence() {
	CpModelBuilder cp_model;
	const int num_nodes = 5;

	// Dummy 5x5 distance matrix
	std::vector<std::vector<int64_t>> costs(num_nodes, std::vector<int64_t>(num_nodes, 0));

	// FIX 1: Change BoolExpr to BoolVar
	std::vector<std::vector<BoolVar>> arcs(num_nodes, std::vector<BoolVar>(num_nodes));

	for (int i = 0; i < num_nodes; ++i) {
		for (int j = 0; j < num_nodes; ++j) {
			if (i == j) continue;
			// Create boolean variables tracking if the route goes directly i -> j
			arcs[i][j] = cp_model.NewBoolVar();
		}
	}

	// 2. Define Rank/Sequence Variables
	std::vector<IntVar> rank;
	for (int i = 0; i < num_nodes; ++i) {
		rank.push_back(cp_model.NewIntVar(Domain(0, num_nodes - 1)));
	}

	// 3. Add Circuit Constraint (Ensures a single Hamiltonian cycle)
	std::vector<std::vector<BoolVar>> circuit_arcs(num_nodes, std::vector<BoolVar>(num_nodes));
	// ... logic to populate circuit_arcs and call cp_model.AddCircuit ...

	// 4. PRECEDENCE CONSTRAINT: Node 1 must come before Node 3
	cp_model.AddLessThan(rank[1], rank[3]);

	// 5. Connect Arcs to Ranks (If node i goes to j, then rank[j] = rank[i] + 1)
	for (int i = 0; i < num_nodes; ++i) {
		for (int j = 0; j < num_nodes; ++j) {
			if (i == j) continue;
			if (j != 0) {
				cp_model.AddEquality(rank[j], rank[i] + 1).OnlyEnforceIf(arcs[i][j]);
			}
		}
	}

	// FIX 6: Build the total objective function using LinearExpr::WeightedSum
	std::vector<BoolVar> flat_arcs;
	std::vector<int64_t> flat_costs;
	for (int i = 0; i < num_nodes; ++i) {
		for (int j = 0; j < num_nodes; ++j) {
			if (i == j) continue;
			flat_arcs.push_back(arcs[i][j]);
			flat_costs.push_back(costs[i][j]);
		}
	}

	LinearExpr objective = LinearExpr::WeightedSum(flat_arcs, flat_costs);
	cp_model.Minimize(objective);

	// Solve
	const CpSolverResponse response = Solve(cp_model.Build());


	if (response.status() == CpSolverStatus::OPTIMAL) {
		std::cout << "Optimal TSP solution found!" << std::endl;
	}
}
*/