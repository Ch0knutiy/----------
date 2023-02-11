#pragma once
#include <vector>
#include <fstream>

using namespace std;

double func(double x);

double lambda(double x);

void assembly(
   int n, int m,
   double gamma,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& G,
   vector<double>& x
);

void local_build(
   int n,
   int k,
   double gamma,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& x
);

void assembly2(
   int n, int m,
   double gamma,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& G,
   vector<double>& x
);

void local_build2(
   int n,
   int k,
   double gamma,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& x
);

void llt_decompose(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
);

void forward_prop(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
);

void backward_prop(
   int n,
   int m,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
);

void slau_llt(
   int n, int m,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<double>& G
);

void llt_decompose2(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
);

void gauss(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
);

void slau_llt2(
   int n, int m,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<double>& G
);

void boundary(
   int n,
   int m,
   vector<double>& di,
   vector<double>& G,
   vector<double>& x
   );

void output(
   int n,
   vector<double>& x,
   vector<double>& G
   );

void output2(
   int n,
   vector<double>& x,
   vector<double>& G
);