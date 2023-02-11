#include <iostream>
#include <fstream>
#include <vector>

#include "utils.h"

using namespace std;

int main()
{
   ifstream input;
   input.open("in.txt");
   int n;
   int m;
   double gamma;
   input >> n >> gamma;
   m = 2 * n + 1;

   vector<double> di(m);
   vector<int> ig(m + 1);
   vector<double> gg;
   vector<double> G(m);
   vector<double> F(m);

   vector<vector<double>> B(3, vector<double>(3));
   vector<vector<double>> C(3, vector<double>(3));

   vector<double> x;
   x.reserve(3 * n - 1);
   for (int i = 0; i <= n; i++)
   {
      double value;
      input >> value;
      x.push_back(value);
   }
   assembly(n, m, gamma, ig, di, gg, B, C, F, G, x);
   boundary(n, m, di, G, x);
   slau_llt(n, m, ig, di, gg, G);
   output(n, x, G);
   return 0;
}

//int main()
//{
//   ifstream input;
//   input.open("in.txt");
//   int n;
//   int m;
//   double gamma;
//   input >> n >> gamma;
//   m = 3 * n + 1;
//
//   vector<double> di(m);
//   vector<int> ig(m + 1);
//   vector<double> gg;
//   vector<double> G(m);
//   vector<double> F(m);
//
//   vector<vector<double>> B(4, vector<double>(4));
//   vector<vector<double>> C(4, vector<double>(4));
//
//   vector<double> x;
//   x.reserve(4 * n - 1);
//   for (int i = 0; i <= n; i++)
//   {
//      double value;
//      input >> value;
//      x.push_back(value);
//   }
//   assembly2(n, m, gamma, ig, di, gg, B, C, F, G, x);
//   boundary(n, m, di, G, x);
//   slau_llt2(n, m, ig, di, gg, G);
//   output2(n, x, G);
//   return 0;
//}