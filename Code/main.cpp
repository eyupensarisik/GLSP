#include "Common.h"
#include "DataPrep.h"
#include "GLSP_gen.h"
#include "TwoPhase.h"
#include "TwoPhase_Callback.h"
#include <chrono>

typedef vector<vector<double>> Matrix;
using namespace std::chrono;

#if 1

int main(int nParams, char* params[])
{
	string inputFileName = (nParams > 1 ? params[1] : "../../Data/SingleMachineCapVar2/Data1-15-15-0.6-0.5-100-100-100-0.dat");
	string inputFileNameOnly;

	auto fileNameBegin = inputFileName.find_last_of("/\\");

	if (fileNameBegin == string::npos)
		inputFileNameOnly = inputFileName;
	else
		inputFileNameOnly = inputFileName.substr(fileNameBegin + 1);

	auto fileNameEnd = inputFileNameOnly.find_last_of(".");

	if (fileNameEnd != string::npos)
		inputFileNameOnly = inputFileNameOnly.substr(0, fileNameEnd);

	string outputFileName = (nParams > 2 ? params[2] : inputFileNameOnly + "_result.csv");
	string summaryFileName = (nParams > 3 ? params[3] : "summary.csv");

	string parameterFileName = (nParams > 4 ? params[4] : "../Run/Parameters/TwoPhase.txt");

	int pOverride = 3;// stoi((nParams > 5 ? params[5] : "0"));
	int tOverride = 3;//stoi((nParams > 6 ? params[6] : "0"));

	ParameterMap Parameters;
	ReadParameterMapFromFile(Parameters, parameterFileName);

	ofstream summaryFile;

	if (nParams > 3)
	{
		summaryFile.open(summaryFileName.c_str(), ios::app);
		if (!summaryFile)
			cerr << "Unable to open summary file in append mode " << summaryFileName << endl;
	}
	else
	{
		summaryFile.open(summaryFileName.c_str(), ios::trunc);
		if (!summaryFile)
			cerr << "Unable to create summary file " << summaryFileName << endl;
		else if (GetParameterValue(Parameters, "GLSP_Standard") || GetParameterValue(Parameters, "GLSP_NF"))
		{
			summaryFile << "Name,P,T,GLSP_CPU,GLSP_LB,GLSP_UB,GLSP_Gap" << endl;
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE"))
		{
			summaryFile << "Name,P,T,TwoPhase_Iter,TwoPhase_Cut,TwoPhase_CPU,TwoPhase_Obj,MP_CPU,SP_Cons_CPU,SP_Solve_CPU" << endl;
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE_CALLBACK"))
		{
			summaryFile << "Name,P,T,TwoPhase_Cut,TwoPhase_CPU,TwoPhase_Obj,TwoPhase_Callback_CPU,TwoPhase_Callback_Cut" << endl;
		}
	}

	ProductPeriods PP;
	ifstream file(inputFileName);
	if (file)
	{
		cout << endl << endl << "Reading " << inputFileName << endl;
		PP.ReadData(file);
		if (pOverride > 0 && tOverride > 0)
			PP.Resize(pOverride, tOverride);

		cout << "Parameter " << parameterFileName << endl;

		int timeLimit = GetParameterValue(Parameters, "TIME_LIMIT");

		string paramFileName = parameterFileName;
		if (summaryFile)
		{
			auto pos = parameterFileName.find_last_of("/\\");
			if (pos != string::npos)
				paramFileName = parameterFileName.substr(pos + 1);
		}

		if (GetParameterValue(Parameters, "GLSP_Standard"))
		{
			double CPUTime_GLSP = 0;
			double SolveTime_GLSP = 0;
			double ObjVal = 0;
			double RelativeGap_GLSP = 0;

			GLSP glsp(PP, Parameters);


			cout << "Started solving GLSP_Standard" << endl;
			glsp.SetupModel_S(timeLimit);
			glsp.Solve(timeLimit);
			cout << "Finished solving GLSP_Standard. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;

			CPUTime_GLSP = glsp.GetCPUTime();
			ObjVal = glsp.GetUB();
			RelativeGap_GLSP = glsp.GetGap();
			glsp.GetSolutions(PP.P, PP.T, PP.S);

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< CPUTime_GLSP << ","
					<< glsp.GetLB() << ","
					<< glsp.GetUB() << ","
					<< ObjVal << ","
					<< RelativeGap_GLSP << endl;
			}
		}
		else if (GetParameterValue(Parameters, "GLSP_NF"))
		{
			double CPUTime_GLSP = 0;
			double SolveTime_GLSP = 0;
			double ObjVal = 0;
			double RelativeGap_GLSP = 0;

			GLSP glsp(PP, Parameters);

			cout << "Started solving GLSP_NF" << endl;
			glsp.SetupModel_NF(timeLimit);
			glsp.Solve(timeLimit);
			cout << "Finished solving GLSP_NF. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;

			CPUTime_GLSP = glsp.GetCPUTime();
			ObjVal = glsp.GetUB();
			RelativeGap_GLSP = glsp.GetGap();
			glsp.GetSolutions(PP.P, PP.T, PP.S);

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< CPUTime_GLSP << ","
					<< glsp.GetLB() << ","
					<< glsp.GetUB() << ","
					<< ObjVal << ","
					<< RelativeGap_GLSP << endl;
			}
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE")){
			TwoPhase tp(PP, Parameters);

			double TwoPhase_Iter;
			double MP_CPU = 0;
			double TwoPhase_CPU = 0;
			double TwoPhase_Obj = 0;
			double TwoPhase_Cut = 0;
			double SP_Cons_CPU = 0;
			double SP_Solve_CPU = 0;

			cout << "Started solving TwoPhase" << endl;
			tp.SetupModel();
			tp.Solve(timeLimit, &TwoPhase_Iter, &TwoPhase_Cut, &TwoPhase_CPU, &TwoPhase_Obj, &SP_Cons_CPU, &SP_Solve_CPU, &MP_CPU);
			cout << "Finished solving TwoPhase. UB: " << tp.GetUB() << " LB: " << tp.GetLB() << endl;

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< TwoPhase_Iter << ","
					<< TwoPhase_Cut << ","
					<< TwoPhase_CPU << ","
					<< TwoPhase_Obj << ","
					<< MP_CPU << ","
					<< SP_Cons_CPU << ","
					<< SP_Solve_CPU << endl;
			}
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE_CALLBACK"))
		{
			TwoPhaseC tpC(PP, Parameters);

			double TwoPhase_Callback_CPU = 0;
			double TwoPhase_Callback_Call = 0;
			double TwoPhase_CPU = 0;
			double TwoPhase_Obj = 0;
			double TwoPhase_Cut = 0;

			cout << "Started solving TwoPhase_Callback" << endl;
			tpC.SetupModel();
			tpC.Solve(timeLimit);
			cout << "Finished solving TwoPhase_Callback. UB: " << tpC.GetUB() << " LB: " << tpC.GetLB() << endl;

			TwoPhase_CPU = tpC.GetCPUTime();
			TwoPhase_Obj = tpC.GetLB();
			TwoPhase_Callback_CPU += tpC.GetCallbackCPU();
			TwoPhase_Callback_Call += tpC.GetCallCount();
			TwoPhase_Cut += tpC.GetCutCount();

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< TwoPhase_Cut << ","
					<< TwoPhase_CPU << ","
					<< TwoPhase_Obj << ","
					<< TwoPhase_Callback_CPU << ","
					<< TwoPhase_Callback_Call << endl;
			}

		}

	}
	else
	{
		cout << "Input file could not be opened!" << endl;
		return -1;
	}

	return 0;
}
#endif