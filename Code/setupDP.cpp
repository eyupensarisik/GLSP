#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <chrono> 
#include "setupDP.h"

using namespace std;

setupDP::setupDP()
{
}

int setupDP::MinSetup(map<pair<int, set<int>>, int>& Cache, const Matrix& Setups, int Prev, set<int>& Current)
{
    auto [it, inserted] = Cache.try_emplace(make_pair(Prev, Current));
    if (inserted)
    {
        int Min = INT_MAX;
        for (auto i : Current)
        {
            if (i == Prev)
                continue;
            else
            {
                auto CurrentCopy = Current;
                CurrentCopy.erase(i);
                int Setup = 0;
                if (Prev != -1) {
                    Setup = Setups[Prev][i];
                }
                if (!CurrentCopy.empty())
                    Setup += MinSetup(Cache, Setups, i, CurrentCopy);
                if (Setup < Min) {
                    Min = Setup;
                }
            }
        }
        it->second = Min;
    }
    else
        int i = 4;
    return it->second;
}

int setupDP::MinSetupMiddle(CacheMP& cache, const Matrix& Setups, int First, ProductSet& Middle, int Last)
{
    auto [it, inserted] = cache.try_emplace(CacheKey(First, Middle, Last));
    if (inserted)
    {
        int Min = INT_MAX;
        if (First == -1)
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
        else if (Last == -1)
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
        else if (Middle.empty())
        {
            Min = Setups[First][Last];
        }
        else if (Middle.size() == 1)
        {
            int Mid = *Middle.begin();
            Min = Setups[First][Mid] + Setups[Mid][Last];
        }
        else
        {
            auto MiddleCopy = Middle;
            for (auto i : Middle)
            {
                MiddleCopy.erase(i);
                int Setup = Setups[First][i];
                Setup += MinSetupMiddle(cache, Setups, i, MiddleCopy, Last);
                if (Setup < Min)
                    Min = Setup;
                MiddleCopy.insert(i);
            }
        }
        it->second = Min;
    }
    return it->second;
}

int setupDP::MinSetupMultiPeriod(CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last)
{
    if (Middle.empty())
        return Setups[First][Last];
    else if (Middle.size() == 1)
        return MinSetupMiddle(cache, Setups, First, Middle.front(), Last);
    else
    {
        auto Split = Middle.size() / 2;
        vector<ProductSet> Left(Middle.begin(), Middle.begin() + Split);

        map<int, int> LeftSetups;
        auto LeftLastCopy = Left.back();
        for (auto i : LeftLastCopy)
        {
            Left.back().erase(i);
            LeftSetups[i] = MinSetupMultiPeriod(cache, Setups, First, Left, i);
            Left.back().insert(i);
        }

        vector<ProductSet> Right(Middle.begin() + Split, Middle.end());
        map<int, int> RightSetups;
        auto RightFirstCopy = Right.front();
        for (auto i : RightFirstCopy)
        {
            Right.front().erase(i);
            RightSetups[i] = MinSetupMultiPeriod(cache, Setups, i, Right, Last);
            Right.front().insert(i);
        }

        int MinSetup = INT_MAX;
        for (auto l : LeftSetups)
            for (auto r : RightSetups)
            {
                int Setup = Setups[l.first][r.first];
                Setup += l.second;
                Setup += r.second;
                if (Setup < MinSetup)
                    MinSetup = Setup;
            }
        return MinSetup;
    }
}