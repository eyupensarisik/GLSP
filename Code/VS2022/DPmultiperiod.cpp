#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <chrono>     

using namespace std;

typedef vector<vector<int>> SetupMatrix;
typedef set<int> ProductSet;
typedef tuple<int, ProductSet, int> CacheKey;
typedef map<CacheKey, int> Cache;

int MinSetup(Cache& cache, const SetupMatrix& Setups, int First, ProductSet& Middle, int Last)
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
                auto Setup = MinSetup(cache, Setups, i, MiddleCopy, Last);
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
                auto Setup = MinSetup(cache, Setups, First, MiddleCopy, i);
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
                Setup += MinSetup(cache, Setups, i, MiddleCopy, Last);
                if (Setup < Min)
                    Min = Setup;
                MiddleCopy.insert(i);
            }
        }
        it->second = Min;
    }
    return it->second;
}


int MinSetupMultiPeriod(Cache& cache, const SetupMatrix& Setups, int First, vector<ProductSet>& Middle, int Last)
{
    if (Middle.empty())
        return Setups[First][Last];
    else if (Middle.size() == 1)
        return MinSetup(cache, Setups, First, Middle.front(), Last);
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

#if 0
int main()
{
    int n = 4;
    int T = 4;

    vector<int> nums;
    for (int i = 0; i < n; ++i)
        nums.emplace_back(i);

    SetupMatrix Setups;
    Setups.resize(n);
    for (int i = 0; i < n; ++i)
        Setups[i].resize(n, 0);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                Setups[i][j] = i + j;

    Cache cache;

    ProductSet current;
    for (int i = 0; i < n; ++i)
        current.insert(i);

    {
        cout << "Finding minimum setup" << endl;

        auto start = std::chrono::system_clock::now();
        int setup = MinSetup(cache, Setups, -1, current, -1);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end - start;

        cout << "Done. Minimum setup: " << setup << endl;
        cout << "Cache size: " << cache.size() << endl;
        cout << "Time: " << duration.count() << endl;
    }

    vector<ProductSet> ProductPPs;
    int Case = 0;


    switch (Case) {
    case 0:
        for (int t = 0; t < T; ++t)
            ProductPPs.push_back(current);
        break;
    case 1:
        for (int t = 0; t < T; ++t) {
            ProductSet current;
            for (int i = t; i < n; ++i) {
                current.insert(i);
            }
            ProductPPs.push_back(current);
        }
        break;
    case 2:
        for (int t = 0; t < T - 1; ++t) {
            ProductSet current;
            for (int i = t; i < n; ++i) {
                current.insert(i);
            }
            ProductPPs.push_back(current);
        }
        ProductPPs.push_back({});
        break;
    }

    {
        cout << "Finding minimum setup multiperiod" << endl;

        auto start = std::chrono::system_clock::now();
        int setup = MinSetupMultiPeriod(cache, Setups, -1, ProductPPs, -1);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end - start;

        cout << "Done. Minimum setup: " << setup << endl;
        cout << "Cache size: " << cache.size() << endl;
        cout << "Time: " << duration.count() << endl;
    }

    return 0;
}
#endif