#include "Common.h"
#include "DataPrep.h"
#include "GLSP_gen.h"

#if 1

int main(int nParams, char* params[])
{
	string inputFileName = (nParams > 1 ? params[1] : "../../Data/SingleMachineCapVar2/Data1-15-15-0.6-0.5-50-100-100-0.dat");
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

	int pOverride = stoi((nParams > 4 ? params[4] : "0"));
	int tOverride = stoi((nParams > 5 ? params[5] : "0"));

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
		else
		{
			summaryFile << "Name,P,T,GLSP_CPU,GLSP_Solve,GLSP_LB,GLSP_UB,GLSP_Gap"	<< endl;
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

		int timeLimit = 600;

		double CPUTime_GLSP = 0;
		double SolveTime_GLSP = 0;
		double ObjVal = 0;	
		double RelativeGap_GLSP = 0;

		GLSP glsp(PP);

		cout << "Started solving GLSP" << endl;

		clock_t startGLSP = clock(), finishGLSP; //Start clock for GLSP

		glsp.SetupModel(timeLimit);
		glsp.Solve(timeLimit);
		
		cout << "Finished solving GLSP. UB: " << glsp.GetUB() << " LB: " << glsp.GetLB() << endl;
		
		finishGLSP = clock(); //Stop clock for GLSP

		CPUTime_GLSP = glsp.GetCPUTime();
		SolveTime_GLSP = (finishGLSP - startGLSP) / CLOCKS_PER_SEC;
		ObjVal = glsp.GetUB();
		RelativeGap_GLSP = glsp.GetGap();
		glsp.GetSolutions(PP.P, PP.T, PP.S);

		if (summaryFile)
		{
			summaryFile << outputFileName << ","
				<< PP.P << ","
				<< PP.T << ","
				<< CPUTime_GLSP << ","
				<< SolveTime_GLSP << ","
				<< glsp.GetLB() << ","
				<< glsp.GetUB() << ","
				<< ObjVal << ","
				<< RelativeGap_GLSP << endl;
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