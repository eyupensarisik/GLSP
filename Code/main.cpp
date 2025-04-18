#include "Common.h"
#include "DataPrep.h"
#include "GLSP_gen.h"
#include "TwoPhase_Callback.h"
#include <chrono>

typedef vector<vector<double>> Matrix;
using namespace std::chrono;

#if 1

int main(int nParams, char* params[])
{
	string inputFileName = (nParams > 1 ? params[1] : "../../Data/SingleMachineCapVar2/Data1-15-15-0.8-0.5-100-100-100-2.dat");
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

	string parameterFileName = (nParams > 4 ? params[4] : "../Run/Parameters/TwoPhase_Callback.txt");

	int pOverride = stoi((nParams > 5 ? params[5] : "0"));
	int tOverride = stoi((nParams > 6 ? params[6] : "0"));

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
		else if (GetParameterValue(Parameters, "GLSP_Standard") || GetParameterValue(Parameters, "GLSP_NF") || GetParameterValue(Parameters, "GLSP_CC") )
		{
			summaryFile << "Name,P,T,GLSP_CPU,GLSP_LB,GLSP_UB,GLSP_Gap,Nconsts,Nvars" << endl;
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE_CALLBACK") || GetParameterValue(Parameters, "TWO_PHASE_SO"))
		{
			summaryFile << "Name,P,T,TwoPhase_Cut,TwoPhase_CPU,TwoPhase_LB,TwoPhase_UB,TwoPhase_Obj,TwoPhase_Gap,TwoPhase_Callback_CPU,TwoPhase_Callback_Cut,SP_Cons_CPU,SP_Solve_CPU,Nconsts,Nvars" << endl;
		}
	}

	ProductPeriods PP;
	ifstream file(inputFileName);
	int SetupCostLevel = GetParameterValue(Parameters, "SETUP_COST_LEVEL");
	if (file)
	{
		cout << endl << endl << "Reading " << inputFileName << endl;
		PP.ReadData(file, SetupCostLevel);
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
			GLSP glsp(PP, Parameters);


			cout << "Started solving GLSP_Standard" << endl;
			glsp.SetupModel_S(timeLimit);
			glsp.Solve(timeLimit);
			cout << "Finished solving GLSP_Standard. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;

			double CPUTime_GLSP = glsp.GetCPUTime();
			double ObjVal = glsp.GetUB();
			double RelativeGap_GLSP = glsp.GetGap();
			int Nconsts = glsp.GetConsts();
			int Nvars = glsp.GetVars();
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
					<< RelativeGap_GLSP << ","
					<< Nconsts << ","
					<< Nvars << endl;
			}
		}
		else if (GetParameterValue(Parameters, "GLSP_NF"))
		{
			GLSP glsp(PP, Parameters);

			cout << "Started solving GLSP_NF" << endl;
			glsp.SetupModel_NF(timeLimit);
			glsp.Solve(timeLimit);
			cout << "Finished solving GLSP_NF. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;

			double CPUTime_GLSP = glsp.GetCPUTime();
			double ObjVal = glsp.GetUB();
			double RelativeGap_GLSP = glsp.GetGap();
			int Nconsts = glsp.GetConsts();
			int Nvars = glsp.GetVars();
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
					<< RelativeGap_GLSP << ","
					<< Nconsts << ","
					<< Nvars << endl;
			}
		}
		else if (GetParameterValue(Parameters, "GLSP_CC"))
		{
			GLSP glsp(PP, Parameters);

			cout << "Started solving GLSP_CC" << endl;
			glsp.SetupModel_CC(timeLimit);
			glsp.Solve(timeLimit);
			cout << "Finished solving GLSP_CC. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;

			double CPUTime_GLSP = glsp.GetCPUTime();
			double ObjVal = glsp.GetUB();
			double RelativeGap_GLSP = glsp.GetGap();
			int Nconsts = glsp.GetConsts();
			int Nvars = glsp.GetVars();
			glsp.GetSolutions_CC(PP.P, PP.T, PP.S);

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< CPUTime_GLSP << ","
					<< glsp.GetLB() << ","
					<< glsp.GetUB() << ","
					<< ObjVal << ","
					<< RelativeGap_GLSP << ","
					<< Nconsts << ","
					<< Nvars << endl;
			}
		}
		else if (GetParameterValue(Parameters, "TWO_PHASE_CALLBACK"))
		{
			TwoPhaseC tpC(PP, Parameters);

			cout << "Started solving TwoPhase_Callback" << endl;
			tpC.SetupModel(timeLimit);
			tpC.AddValidInequalities();
			tpC.Solve(timeLimit);
			
			cout << "Finished solving TwoPhase_Callback. UB: " << tpC.GetUB() << " LB: " << tpC.GetLB() << endl;

			double TwoPhase_CPU = tpC.GetCPUTime();
			double TwoPhase_LB = tpC.GetLB();
			double TwoPhase_UB = tpC.GetUB();
			double TwoPhase_Obj = tpC.GetUB();
			double TwoPhase_Gap = tpC.GetGap();
			double TwoPhase_Callback_CPU = tpC.GetCallbackCPU();
			double TwoPhase_Callback_Call = tpC.GetCallCount();
			double TwoPhase_Cut = tpC.GetCutCount();
			double SP_Cons_CPU = tpC.GetSPconsCPU();
			double SP_Solve_CPU = tpC.GetSPsolveCPU();
			int Nconsts = tpC.GetConsts();
			int Nvars = tpC.GetVars();

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< PP.P << ","
					<< PP.T << ","
					<< TwoPhase_Cut << ","
					<< TwoPhase_CPU << ","
					<< TwoPhase_LB << ","
					<< TwoPhase_UB << ","
					<< TwoPhase_Obj << ","
					<< TwoPhase_Gap << ","
					<< TwoPhase_Callback_CPU << ","
					<< TwoPhase_Callback_Call << ","
					<< SP_Cons_CPU << ","
					<< SP_Solve_CPU << ","
					<< Nconsts << ","
					<< Nvars << endl;
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