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

int setupDP::MinSetupMultiPeriod(CacheMP& cache, const Matrix& Setups, int First, vector<ProductSet>& Middle, int Last)
{
    if (Middle.empty())
        return SetupBetween(Setups, First, Last);
    else if (Middle.size() == 1)
        return MinSetupMiddle(cache, Setups, First, Middle.front(), Last);
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
            LeftSetups[i] = MinSetupMultiPeriod(cache, Setups, First, Left, i);
            Left.back().insert(i);
        }
        }
        else {
            for (auto i : LeftLastCopy)
                LeftSetups[i] = MinSetupMultiPeriod(cache, Setups, First, Left, i);
        }

        vector<ProductSet> Right(Middle.begin() + Split, Middle.end());
        map<int, int> RightSetups;
        auto RightFirstCopy = Right.front();
        if (RightFirstCopy.size() > 1){
            for (auto i : RightFirstCopy)
            {
                Right.front().erase(i);
                RightSetups[i] = MinSetupMultiPeriod(cache, Setups, i, Right, Last);
                Right.front().insert(i);
            }
        }
        else {
            for (auto i : RightFirstCopy)
                RightSetups[i] = MinSetupMultiPeriod(cache, Setups, i, Right, Last);
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