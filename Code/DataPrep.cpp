#include "DataPrep.h"
#include <iostream>
#include <algorithm>
#include <sstream>

void ProductPeriods::SetDimensions(int p, int t)
{
	S_b.resize(T);
	S_f.resize(T);
	L.resize(T);
	K.resize(T);

	Products.resize(P);
	for (int i = 0; i < p; ++i)
	{
		Products[i].Index = i;
		Products[i].d.resize(T);
		Products[i].s_pr.resize(P);
	}
}

void ProductPeriods::Resize(int NewP, int NewT)
{
	auto ResizeFunction = [](auto& Vect, int NewSize)
	{
		int CurrentSize = Vect.size();
		if (NewSize < CurrentSize)
			Vect.resize(NewSize);
		else if (NewSize > CurrentSize)
		{
			Vect.reserve(NewSize);
			for (int i = 0; i < NewSize - CurrentSize; ++i)
				Vect.push_back(Vect[i % CurrentSize]);
		}
	};

	ResizeFunction(Products, NewP);
	ResizeFunction(S_b, NewT);
	ResizeFunction(S_f, NewT);
	ResizeFunction(L, NewT);
	ResizeFunction(K, NewT);
	for (auto& p : Products){
		ResizeFunction(p.d, NewT);
		ResizeFunction(p.s_pr, NewP);
	}

	if (NewT > T)
	{
		for (auto& p : Products){
			p.a *= double(NewT) / T;
			p.m *= double(NewT) / T;
			p.h *= double(NewT) / T;
			p.f *= double(NewT) / T;
			p.c *= double(NewT) / T;
		}
	}

	P = NewP;
	T = NewT;
	S = NewP * NewT;

	for (int t = 0; t < T; ++t) {
		S_b[t] = P * t;
		S_f[t] = P * (t + 1);
		L[t] = (P * (t + 1)) - 1;
	}
}

void ProductPeriods::ReadData(ifstream& file, int& SetupCostLevel)
{
	string DummyToken;
	char digit;
	vector<int> digits;	

	file >> DummyToken >> DummyToken >> P >> DummyToken;
	file >> DummyToken >> DummyToken >> T >> DummyToken;
	S = P * T;

	file >> DummyToken >> DummyToken >> DummyToken;

	SetDimensions(P, T);

	file >> DummyToken >> DummyToken >> digit;	
	for (int p = 0; p < P; ++p)
		file >> Products[p].h; 
	file >> DummyToken >> DummyToken;

	file >> DummyToken >> digit >> digit;
	for (int p = 0; p < P; ++p)
		file >> Products[p].a;
	file >> DummyToken;

	getline(file, DummyToken, ';');

	file >> DummyToken >> DummyToken >> digit >> digit;

	for (int t = 0; t < T; ++t)
		file >> K[t];
	file >> DummyToken >> DummyToken;

	file >> DummyToken;
	for (int p = 0; p < P; ++p){
		file >> digit >> digit >> digit;
		for (int t = 0; t < T; ++t)
			file >> Products[p].d[t];
	}
	file >> DummyToken;

	getline(file, DummyToken, ';');

	file >> DummyToken >> DummyToken;
	for (int p = 0; p < P; ++p) {
		file >> digit >> digit >> digit;
		for (int r = 0; r < P; ++r)
			file >> Products[p].s_pr[r];
	}
	file >> DummyToken;
	
	file >> DummyToken >> DummyToken >> DummyToken >> digit >> digit;
	for (int p = 0; p < P; ++p)
		file >> Products[p].m;
	file >> DummyToken;

	for (int p = 0; p < P; ++p){
		Products[p].f = 0;
		Products[p].c = 0;
	}

	for (int t = 0; t < T; ++t) {
		S_b[t] = P * t;
		S_f[t] = P * (t + 1);
		L[t] = (P * (t + 1)) - 1;
	}

	if (SetupCostLevel > 1){
		for (int p = 0; p < P; ++p) {
			for (int r = 0; r < P; ++r)
				Products[p].s_pr[r] = ceil(Products[p].s_pr[r] / SetupCostLevel);
		}
	}

}

string to_string(const IntegerVector& vect)
{
	stringstream str;
	str << "[";
	for (auto i : vect)
		str << i << " ";
	str << "]";
	return str.str();
}