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
}

void ProductPeriods::ReadData(ifstream& file)
{
	const int MAX_LINE_LENGTH = 1000;
	char DummyLine[MAX_LINE_LENGTH];
	string DummyToken;
	char digit;
	vector<int> digits;

	file >> DummyToken >> DummyToken >> P >> DummyToken;
	file >> DummyToken >> DummyToken >> T >> DummyToken;
	S = P * T;
	M = 9999;
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

	file.getline(DummyLine, MAX_LINE_LENGTH); file.getline(DummyLine, MAX_LINE_LENGTH);

	file >> DummyToken >> DummyToken >> DummyToken >> digit >> digit;
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

	file.getline(DummyLine, MAX_LINE_LENGTH); file.getline(DummyLine, MAX_LINE_LENGTH);

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

	cout << "Value of DummyToken is " << DummyToken << endl;
	cout << "Value of digit is " << digit << endl;

	cout << "Value of P is " << P << endl;
	cout << "Value of T is " << T << endl;
	cout << "Value of S is " << S << endl;

	for (int p = 0; p < P; ++p)
		cout << "Value of h[" << p << "] is " << Products[p].h << endl;
	for (int p = 0; p < P; ++p)
		cout << "Value of a[" << p << "] is " << Products[p].a << endl;
	for (int p = 0; p < P; ++p)
		for (int t = 0; t < T; ++t)
			cout << "Value of d[" << p << "][" << t << "] is " << Products[p].d[t] << endl;
	for (int p = 0; p < P; ++p)
		for (int r = 0; r < P; ++r)
			cout << "Value of s_pr[" << p << "][" << r << "] is " << Products[p].s_pr[r] << endl;
	for (int p = 0; p < P; ++p)
		cout << "Value of m[" << p << "] is " << Products[p].m << endl;
	for (int t = 0; t < T; ++t) {
		cout << "Value of S_b[" << t << "] is " << S_b[t] << endl;
		cout << "Value of S_f[" << t << "] is " << S_f[t] << endl;
	}
	for (int t = 0; t < T; ++t)
		cout << "Value of L[" << t << "] is " << L[t] << endl;
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