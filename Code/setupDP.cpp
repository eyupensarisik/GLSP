#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <chrono> 
#include "setupDP.h"
#include "DataPrep.h"

using namespace std;

setupDP::setupDP()
{
}

setupDP::~setupDP()
{
	TSPenv.end();
}

int setupDP::SetupBetween(const Matrix& Setups, int i, int j)
{
    if (i >= 0 && i < Setups.size() && j >= 0 && j < Setups[i].size())
        return Setups[i][j];
    return 0;
}

int setupDP::MinSetupMiddle(CacheMP& cache, const Matrix& Setups, int First, ProductSet& Middle, int Last)
{
    auto [it, inserted] = cache.try_emplace(CacheKey(First, Middle, Last));
    if (inserted)
    {
        int Min = INT_MAX;
        if (First == -1)
        {
            if (Middle.empty())
                Min = 0;
            else
            {
                auto MiddleCopy = Middle;
                for (auto i : Middle)
                {
                    MiddleCopy.erase(i);
                    auto Setup = MinSetupMiddle(cache, Setups, i, MiddleCopy, Last);
                    if (Setup < Min)
                        Min = Setup;
                    MiddleCopy.insert(i);
                }
            }
        }
        else if (Last == -1)
        {
            if (Middle.empty())
                Min = 0;
            else
            {
                auto MiddleCopy = Middle;
                for (auto i : Middle)
                {
                    MiddleCopy.erase(i);
                    auto Setup = MinSetupMiddle(cache, Setups, First, MiddleCopy, i);
                    if (Setup < Min)
                        Min = Setup;
                    MiddleCopy.insert(i);
                }
            }
        }
        else if (Middle.empty())
        {
            Min = SetupBetween(Setups, First, Last);
        }
        else if (Middle.size() == 1)
        {
            int Mid = *Middle.begin();
            Min = SetupBetween(Setups, First, Mid) + SetupBetween(Setups, Mid, Last);
        }
        else
        {
            auto MiddleCopy = Middle;
            for (auto i : Middle)
            {
                MiddleCopy.erase(i);
                int Setup = SetupBetween(Setups, First, i);
                Setup += MinSetupMiddle(cache, Setups, i, MiddleCopy, Last);
                if (Setup < Min)
                    Min = Setup;
                MiddleCopy.insert(i);
            }
        }
        if (Min == INT_MAX)
            int i = 4;
        it->second = Min;
    }
    return it->second;
}

int setupDP::MinSetupMultiPeriod(ProductPeriods& PPs, CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last)
{
    if (Middle.empty())
        return SetupBetween(Setups, First, Last);
    else if (Middle.size() == 1){
        //cout << MinSetupMiddle(cache, Setups, First, Middle.front(), Last) << endl;
        //if (MinSetupMiddle(cache, Setups, First, Middle.front(), Last) == 0)
			//cout << "Here" << endl;
        //return MinSetupMiddle(cache, Setups, First, Middle.front(), Last);
		return SetupTSPModel(PPs, Setups, First, Middle.front(), Last);
    }
    else
    {
        auto Split = Middle.size() / 2;
        vector<ProductSet> Left(Middle.begin(), Middle.begin() + Split);

        map<int, int> LeftSetups;
        auto LeftLastCopy = Left.back();
        if (LeftLastCopy.size() > 1) {
        for (auto i : LeftLastCopy)
        {
            Left.back().erase(i);
            LeftSetups[i] = MinSetupMultiPeriod(PPs, cache, Setups, First, Left, i);
            Left.back().insert(i);
        }
        }
        else {
            for (auto i : LeftLastCopy)
                LeftSetups[i] = MinSetupMultiPeriod(PPs, cache, Setups, First, Left, i);
        }

        vector<ProductSet> Right(Middle.begin() + Split, Middle.end());
        map<int, int> RightSetups;
        auto RightFirstCopy = Right.front();
        if (RightFirstCopy.size() > 1){
            for (auto i : RightFirstCopy)
            {
                Right.front().erase(i);
                RightSetups[i] = MinSetupMultiPeriod(PPs, cache, Setups, i, Right, Last);
                Right.front().insert(i);
            }
        }
        else {
            for (auto i : RightFirstCopy)
                RightSetups[i] = MinSetupMultiPeriod(PPs, cache, Setups, i, Right, Last);
        }

        int MinSetup = INT_MAX;
        if (LeftSetups.empty() && RightSetups.empty())
            MinSetup = 0;
        else if (LeftSetups.empty())
        {
            for (auto r : RightSetups) 
                if (r.second < MinSetup)
                    MinSetup = r.second;
        }
        else if (RightSetups.empty())
        {
            for (auto l : LeftSetups)
                if (l.second < MinSetup)
                    MinSetup = l.second;
        }
        else
        {
            for (auto l : LeftSetups)
                for (auto r : RightSetups)
                {
                    int Setup = SetupBetween(Setups, l.first, r.first);
                    Setup += l.second;
                    Setup += r.second;
                    if (Setup < MinSetup)
                        MinSetup = Setup;
                }
        }
        return MinSetup;
    }
}

void setupDP::generateCombinations(const vector<int>& nums, set<int>& current, int index, vector<set<int>>& all) {
    // Base case: We've reached the end of the list
    if (index == nums.size()) {
        all.emplace_back(current);
        return;
    }

    // Exclude the current element
    generateCombinations(nums, current, index + 1, all);

    // Include the current element
    current.insert(nums[index]);
    generateCombinations(nums, current, index + 1, all);

    // Backtrack
    current.erase(nums[index]);
}

double setupDP::SetupTSPModel(ProductPeriods& PPs, const Matrix& Setups, int First, ProductSet& Middle, int Last)
{
	double Min = 0;

    if (Middle.size() == 1) {
        int Mid = *Middle.begin();
        Min = SetupBetween(Setups, First, Mid) + SetupBetween(Setups, Mid, Last);   
    }
    else {
        IloModel TSPmodel(TSPenv);
        IloCplex TSPcplex(TSPenv);

        //TSPcplex.setParam(IloCplex::Threads, 1);

        // Setup TSP Subroblem
        TSPcplex.setOut(TSPenv.getNullStream());

        if (First == -1) {
            //Construct Subproblem
            int W = Middle.size() + 1; // number of nodes in the TSP (including start and end)
            set<int> Middle_Last = Middle;
            Middle_Last.emplace(Last);

            // u[i] for i in V, u[0] = 1 (node 1), u[i] in [2,n] for i != 0
            IloNumVarArray u(TSPenv, PPs.P);
            for (auto& j : Middle_Last) {
                char name[32];
                snprintf(name, sizeof(name), "u(%d)", j);
                // 2 <= u_i <= n  (constraint 7)
                u[j] = IloNumVar(TSPenv, 1.0, (double)(W), ILOFLOAT, name);
            }

            e = CreateBoolVarArray2(TSPenv, PPs.P, PPs.P, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise

            IloExpr TSPobj1(TSPenv);
            for (auto& j : Middle_Last)
                for (auto& l : Middle_Last)
                    if (j != l)
                        TSPobj1 += Setups[j][l] * e[j][l];

            IloObjective TSPobj = IloMinimize(TSPenv, TSPobj1);
            TSPmodel.add(TSPobj);
			TSPobj1.end();

            TSPcplex = IloCplex(TSPmodel);

            // Constraint (2)
            for (auto& j : Middle) {
                IloExpr Each1(TSPenv);
                for (auto& l : Middle_Last)
                    if (j != l)
                        Each1 += e[j][l];
                TSPmodel.add(Each1 == 1);
				Each1.end();
            }

            // Constraint (3)
            for (auto& j : Middle_Last) {
                IloExpr Each2(TSPenv);
                for (auto& l : Middle_Last)
                    if (j != l)
                        Each2 += e[l][j];
                TSPmodel.add(Each2 <= 1);
                Each2.end();
            }

            // Constraint (4)
            for (auto& l : Middle) {
                IloExpr Each3(TSPenv);
                for (auto& j : Middle_Last)
                    if (j != l)
                        Each3 += e[j][l];
                IloExpr Each33(TSPenv);
                for (auto& k : Middle_Last)
                    if (l != k)
                        Each33 += e[l][k];

                TSPmodel.add(Each33 - Each3 >= 0);
                Each3.end();
                Each33.end();
            }

            // Constraint (5)
            IloExpr Each4(TSPenv);
            for (auto& l : Middle)
                Each4 += e[l][Last];
            TSPmodel.add(Each4 == 1);

            IloExpr Each44(TSPenv);
            for (auto& l : Middle)
                Each44 += e[Last][l];
            TSPmodel.add(Each44 == 0);
            Each4.end();
            Each44.end();

            // MTZ Formulation

            for (auto& j : Middle_Last) {
                for (auto& l : Middle_Last)
                    if (j != l)
                        TSPmodel.add(u[j] - u[l] + W * e[j][l] <= W - 1);
            }
        }
        else if (Last == -1) {
            //Construct Subproblem
            int W = Middle.size() + 1; // number of nodes in the TSP (including start and end)
            set<int> Middle_First = Middle;
            Middle_First.emplace(First);

            // u[i] for i in V, u[0] = 1 (node 1), u[i] in [2,n] for i != 0
            IloNumVarArray u(TSPenv, PPs.P);
            {
                // u_0 = 0  (constraint 6)
                u[First] = IloNumVar(TSPenv, 0.0, 0.0, ILOFLOAT);
                char name[32];
                snprintf(name, sizeof(name), "u(0)");
                u[First].setName(name);
            }
            for (auto& j : Middle) {
                char name[32];
                snprintf(name, sizeof(name), "u(%d)", j);
                // 2 <= u_i <= n  (constraint 7)
                u[j] = IloNumVar(TSPenv, 1.0, (double)(W), ILOFLOAT, name);
            }

            e = CreateBoolVarArray2(TSPenv, PPs.P, PPs.P, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise

            IloExpr TSPobj1(TSPenv);
            for (auto& j : Middle_First)
                for (auto& l : Middle_First)
                    if (j != l)
                        TSPobj1 += Setups[j][l] * e[j][l];

            IloObjective TSPobj = IloMinimize(TSPenv, TSPobj1);
            TSPmodel.add(TSPobj);
            TSPobj1.end();

            TSPcplex = IloCplex(TSPmodel);

            // Constraint (2)
            for (auto& j : Middle_First) {
                IloExpr Each1(TSPenv);
                for (auto& l : Middle_First)
                    if (j != l)
                        Each1 += e[j][l];
                TSPmodel.add(Each1 <= 1);
                Each1.end();
            }

            // Constraint (3)
            for (auto& j : Middle) {
                IloExpr Each2(TSPenv);
                for (auto& l : Middle_First)
                    if (j != l)
                        Each2 += e[l][j];
                TSPmodel.add(Each2 == 1);
                Each2.end();
            }

            // Constraint (4)
            for (auto& j : Middle) {
                IloExpr Each3(TSPenv);
                for (auto& l : Middle_First)
                    if (j != l)
                        Each3 += e[j][l];
                IloExpr Each33(TSPenv);
                for (auto& k : Middle_First)
                    if (j != k)
                        Each33 += e[k][j];

                TSPmodel.add(Each33 - Each3 >= 0);
                Each3.end();
                Each33.end();
            }

            // Constraint (5)
            IloExpr Each4(TSPenv);
            for (auto& l : Middle)
                Each4 += e[First][l];
            TSPmodel.add(Each4 == 1);
            Each4.end();

            IloExpr Each44(TSPenv);
            for (auto& l : Middle)
                Each44 += e[l][First];
            TSPmodel.add(Each44 == 0);
            Each44.end();

            // MTZ Formulation

            for (auto& j : Middle_First) {
                for (auto& l : Middle_First)
                    if (j != l)
                        TSPmodel.add(u[j] - u[l] + W * e[j][l] <= W - 1);
            }
        }
        else {
            //Construct Subproblem
            int W = Middle.size() + 2; // number of nodes in the TSP (including start and end)
            set<int> Middle_Last = Middle;
            Middle_Last.emplace(Last);

            set<int> Middle_First = Middle;
            Middle_First.emplace(First);

            set<int> Middle_Full = Middle;
            Middle_Full.emplace(First);
            Middle_Full.emplace(Last);

            // u[i] for i in V, u[0] = 1 (node 1), u[i] in [2,n] for i != 0
            IloNumVarArray u(TSPenv, PPs.P);
            {
                // u_0 = 0  (constraint 6)
                u[First] = IloNumVar(TSPenv, 0.0, 0.0, ILOFLOAT);
                char name[32];
                snprintf(name, sizeof(name), "u(0)");
                u[First].setName(name);
            }
            for (auto& j : Middle_Last) {
                char name[32];
                snprintf(name, sizeof(name), "u(%d)", j);
                // 2 <= u_i <= n  (constraint 7)
                u[j] = IloNumVar(TSPenv, 1.0, (double)(W), ILOFLOAT, name);
            }

            e = CreateBoolVarArray2(TSPenv, PPs.P, PPs.P, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise

            IloExpr TSPobj1(TSPenv);
            for (auto& j : Middle_Full)
                for (auto& l : Middle_Full)
                    if (j != l)
                        TSPobj1 += Setups[j][l] * e[j][l];

            IloObjective TSPobj = IloMinimize(TSPenv, TSPobj1);
            TSPmodel.add(TSPobj);
            TSPobj1.end();

            TSPcplex = IloCplex(TSPmodel);

            // Constraint (2)
            for (auto& j : Middle_First) {
                IloExpr Each1(TSPenv);
                for (auto& l : Middle_Full)
                    if (j != l)
                        Each1 += e[j][l];
                TSPmodel.add(Each1 == 1);
                Each1.end();
            }

            // Constraint (3)
            for (auto& j : Middle_Last) {
                IloExpr Each2(TSPenv);
                for (auto& l : Middle_Full)
                    if (j != l)
                        Each2 += e[l][j];
                TSPmodel.add(Each2 == 1);
                Each2.end();
            }
            

            // Constraint (4)
            IloExpr Each3(TSPenv);
            for (auto& l : Middle)
                Each3 += e[First][l];
            TSPmodel.add(Each3 == 1);
            Each3.end();

            // Constraint (5)
            IloExpr Each4(TSPenv);
            for (auto& l : Middle)
                Each4 += e[l][Last];
            TSPmodel.add(Each4 == 1);
            Each4.end();

            // MTZ Formulation

            for (auto& j : Middle) {
                for (auto& l : Middle_Last)
                    if (j != l)
                        TSPmodel.add(u[j] - u[l] + W * e[j][l] <= W - 1);
            }
        }

        try
        {
            //TSPcplex.exportModel("TSP.lp");
            Solved = TSPcplex.solve();
            Min = Solved ? TSPcplex.getObjValue() : DBL_MAX;
            TSPmodel.end();
            TSPcplex.end();
        }
        catch (IloException& ex)
        {
            cout << ex.getMessage() << endl;
        }
    }
    
    return Min;
}