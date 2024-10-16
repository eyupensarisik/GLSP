#include "Common.h"
#include "DataPrep.h"

int main()
{
	IloEnv env;
	IloModel model;
	IloObjective obj;
	IloCplex cplex;
	bool Solved = false;

	NumVarArray2 I; //Inventory level of product p at the end of macro-period t
	NumVarArray2 q; //Amount of product p produced in micro-period s

	BoolVarArray2 x;
	BoolVarArray3 y;
	IloNumVar theta;
	IloNumVar gamma;
	IloNumVarArray theta_t;
	IloNumVarArray gamma_t;


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

	// MP Model

	theta = IloNumVar(env, 0, IloInfinity, "theta");
	gamma = IloNumVar(env, 0, IloInfinity, "gamma");
	theta_t = CreateNumVarArray(env, PP.T, "theta_t", 0, IloInfinity);
	//gamma_t = CreateNumVarArray(env, PP.T, "gamma_t", 0, IloInfinity);

	I = CreateNumVarArray2(env, PP.P, PP.T, "I", 0, IloInfinity);
	q = CreateNumVarArray2(env, PP.P, PP.T, "q", 0, IloInfinity);

	x = CreateBoolVarArray2(env, PP.P, PP.T, "x");
	//y = CreateBoolVarArray3(env, PP.P, PP.P, PP.T, "y");

	theta.setStringProperty("Type", "theta");
	theta.setIntProperty("Product", -1);

}