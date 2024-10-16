#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
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