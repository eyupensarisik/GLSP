#ifndef DATA_PREPARATION
#define DATA_PREPARATION

#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>

using namespace std;

typedef vector<int> IntegerVector;
typedef vector<double> DoubleVector;
typedef set<double> DoubleSet;
typedef vector<DoubleVector> DoubleMatrix;
typedef vector<IntegerVector> IntegerMatrix;

struct Product
{
	int Index;

	int a; //Total capacity (time) used to produce one unit of products
	int m; //Smallest lot-size for products
	int h; //Holding cost of one unit of products in inventory
	int f; //Fixed production cost of products
	int c; //Variable production cost of products

	DoubleVector s_pr; //Setup cost (time) required for transition from product p to r
	DoubleVector d; //Total demand for product p in macro-period t
};

struct ProductPeriods
{
	int P; //The set of products, denoted by p
	int T; //The set of macro-periods, denoted by t
	int S; //The set of micro-periods
	int M; //Sufficiently large number

	vector<Product> Products;
	IntegerVector S_b;	//Starting micro-period for macro-period t
	IntegerVector S_f;	//Ending micro-period for macro-period t
	IntegerVector L;	//The last position (micro-period) in macro-period t
	IntegerVector K;	//Total capacity in macro-period t
	
	void SetDimensions(int p, int t);
	void ReadData(ifstream& file);
	void Resize(int NewP, int NewT);
	
};

string to_string(const IntegerVector& vect);
#endif // !DATA_PREPARATION
