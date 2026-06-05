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

void setupDP::BuildTSPModel(ProductPeriods& PPs, const Matrix& Setups)
{
    if (ModelBuilt)
        return;

    const int P = PPs.P;
    const int M = P + 1;            // constant big-M for MTZ (RHS = M - 1 = P)

    TSPmodel = IloModel(TSPenv);

    // u over all products, pinned inactive [0,0] initially.
    u = IloNumVarArray(TSPenv, P);
    for (int j = 0; j < P; ++j) {
        char name[32];
        snprintf(name, sizeof(name), "u(%d)", j);
        u[j] = IloNumVar(TSPenv, 0.0, 0.0, ILOFLOAT, name);
    }

    // e over all ordered product pairs, pinned inactive [0,0] initially.
    e = CreateBoolVarArray2(TSPenv, P, P, "e"); //1 if job j is scheduled just before job l; and, 0 otherwise
    for (int j = 0; j < P; ++j)
        for (int l = 0; l < P; ++l)
            if (j != l)
                e[j][l].setBounds(0.0, 0.0);

    // Fixed objective over the whole grid (inactive terms pinned to 0).
    IloExpr TSPobj1(TSPenv);
    for (int j = 0; j < P; ++j)
        for (int l = 0; l < P; ++l)
            if (j != l)
                TSPobj1 += Setups[j][l] * e[j][l];
    TSPobj = IloMinimize(TSPenv, TSPobj1);
    TSPmodel.add(TSPobj);
    TSPobj1.end();

    // Per-node degree / balance rows, built once, summed over all of P,
    // initially OFF (free rows). Driven by setBounds() each call.
    OutDeg = IloRangeArray(TSPenv, P);
    InDeg = IloRangeArray(TSPenv, P);
    Bal = IloRangeArray(TSPenv, P);
    for (int j = 0; j < P; ++j) {
        IloExpr outE(TSPenv), inE(TSPenv), balE(TSPenv);
        for (int l = 0; l < P; ++l)
            if (l != j) {
                outE += e[j][l];
                inE += e[l][j];
                balE += e[j][l];
                balE -= e[l][j];
            }
        OutDeg[j] = IloRange(TSPenv, -IloInfinity, outE, IloInfinity);
        InDeg[j] = IloRange(TSPenv, -IloInfinity, inE, IloInfinity);
        Bal[j] = IloRange(TSPenv, -IloInfinity, balE, IloInfinity);
        outE.end(); inE.end(); balE.end();
    }
    TSPmodel.add(OutDeg);
    TSPmodel.add(InDeg);
    TSPmodel.add(Bal);

    // Full P x P MTZ rows, built ONCE, constant big-M, never modified again.
    for (int j = 0; j < P; ++j)
        for (int l = 0; l < P; ++l)
            if (j != l) {
                IloExpr m(TSPenv);
                m += u[j] - u[l] + M * e[j][l];
                TSPmodel.add(IloRange(TSPenv, -IloInfinity, m, (double)(M - 1)));
                m.end();
            }

    // Solver built once; re-solved after each round of bound flips.
    TSPcplex = IloCplex(TSPmodel);
    TSPcplex.setOut(TSPenv.getNullStream());
    //TSPcplex.setParam(IloCplex::Threads, 1);

    ModelBuilt = true;
}

double setupDP::SetupTSPModel(ProductPeriods& PPs, const Matrix& Setups, int First, ProductSet& Middle, int Last)
{
    double Min = 0;

    if (Middle.size() == 1) {
        int Mid = *Middle.begin();
        Min = SetupBetween(Setups, First, Mid) + SetupBetween(Setups, Mid, Last);
        return Min;
    }

    BuildTSPModel(PPs, Setups);   // builds once; no-op on later calls

    // --- Determine case, active node set, and W ---
    int W;
    int caseId;                       // 1: open start, 2: open end, 3: both fixed
    set<int> active = Middle;
    if (First == -1) {                // open start, fixed end
        active.insert(Last);
        W = (int)Middle.size() + 1;
        caseId = 1;
    }
    else if (Last == -1) {            // fixed start, open end
        active.insert(First);
        W = (int)Middle.size() + 1;
        caseId = 2;
    }
    else {                            // fixed start and end
        active.insert(First);
        active.insert(Last);
        W = (int)Middle.size() + 2;
        caseId = 3;
    }

    // --- Reset everything the PREVIOUS call touched back to inactive/off ---
    for (int j : PrevActiveNodes) {
        u[j].setBounds(0.0, 0.0);
        OutDeg[j].setBounds(-IloInfinity, IloInfinity);
        InDeg[j].setBounds(-IloInfinity, IloInfinity);
        Bal[j].setBounds(-IloInfinity, IloInfinity);
        for (int l : PrevActiveNodes)
            if (j != l)
                e[j][l].setBounds(0.0, 0.0);
    }

    // --- Activate this call's e-vars: every ordered pair within `active` ---
    for (int j : active)
        for (int l : active)
            if (j != l)
                e[j][l].setBounds(0.0, 1.0);

    // --- Set u bounds and degree-row senses for this call's case ---
    if (caseId == 1) {
        // open start: no node anchored at 0; all active in [1, W]
        for (int j : active)
            u[j].setBounds(1.0, (double)W);

        for (int j : Middle) {
            OutDeg[j].setBounds(1.0, 1.0);            // exactly one successor
            InDeg[j].setBounds(-IloInfinity, 1.0);    // at most one predecessor
            Bal[j].setBounds(0.0, IloInfinity);       // out - in >= 0
        }
        OutDeg[Last].setBounds(0.0, 0.0);             // sink: no successor
        InDeg[Last].setBounds(1.0, 1.0);              // sink: exactly one predecessor
        // Bal[Last] stays off
    }
    else if (caseId == 2) {
        u[First].setBounds(0.0, 0.0);                 // source anchored at 0
        for (int j : Middle)
            u[j].setBounds(1.0, (double)W);

        for (int j : Middle) {
            OutDeg[j].setBounds(-IloInfinity, 1.0);   // at most one successor
            InDeg[j].setBounds(1.0, 1.0);             // exactly one predecessor
            Bal[j].setBounds(-IloInfinity, 0.0);      // in - out >= 0  (i.e. out - in <= 0)
        }
        OutDeg[First].setBounds(1.0, 1.0);            // source: exactly one successor
        InDeg[First].setBounds(0.0, 0.0);             // source: no predecessor
        // Bal[First] stays off
    }
    else {
        u[First].setBounds(0.0, 0.0);                 // source anchored at 0
        for (int j : Middle)
            u[j].setBounds(1.0, (double)W);
        u[Last].setBounds(1.0, (double)W);

        for (int j : Middle) {
            OutDeg[j].setBounds(1.0, 1.0);            // interior: one successor
            InDeg[j].setBounds(1.0, 1.0);             // interior: one predecessor
        }
        OutDeg[First].setBounds(1.0, 1.0);            // source: one successor
        InDeg[Last].setBounds(1.0, 1.0);             // sink: one predecessor
        // InDeg[First] / OutDeg[Last] stay off; MTZ + the source anchor force
        // in-degree(First) = 0, and the degree budget forces out-degree(Last)=0.
        e[First][Last].setBounds(0.0, 0.0);           // no direct First -> Last edge
    }

    // Remember the active set so the next call resets exactly these.
    PrevActiveNodes.assign(active.begin(), active.end());

    try
    {
        //TSPcplex.exportModel("TSP.lp");
        Solved = TSPcplex.solve();
        Min = Solved ? TSPcplex.getObjValue() : DBL_MAX;
    }
    catch (IloException& ex)
    {
        cout << ex.getMessage() << endl;
        Min = DBL_MAX;
    }

    return Min;
}
