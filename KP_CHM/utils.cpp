#include "utils.h"

double func(double x)
{
   //return -4;
   //return x * x;
   //return -6 * x + x * x;
   return -2 + x;
   //return -12 * x * x + x * x * x;
   //return -8 * x * x + x * x;
   //return -10 * x * x * x + x * x;
}

double lambda(double x)
{
   //return 1;
   //return 0;
   return x;
   //return x * x;
   //return x * x * x;
}

double u(double x)
{
   //return x * x;
   return x;
   //return x * x * x;
}

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
)
{
   //профиль
   ig[0] = 1;
   ig[1] = 1;
   for (int i = 2; i < m + 1; i++)
      ig[i] = ig[i - 1] + i - 1;
   gg.resize(ig[m] - 1);
   for (int i = 0; i < m; i++)
   {
      G[i] = 0.0;
      di[i] = 0.0;
   }

   for (int i = 1; i <= n; i++)
   {
      local_build(n, i - 1, gamma, B, C, F, x);
      for (int j = 0; j < 3; j++)
      {
         di[2 * (i - 1) + j] += B[j][j] + gamma * C[j][j];
         G[2 * (i - 1) + j] += F[j];
      }
      gg[ig[2 * i - 1] - 1] = B[1][0] + gamma * C[1][0];
      gg[ig[2 * i] - 1] = B[2][0] + gamma * C[2][0];
      gg[ig[2 * i]] = B[2][1] + gamma * C[2][1];
   }
}

void local_build(
   int n,
   int k,
   double gamma,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& x
)
{
   vector<vector<double>> G1 = {{37. / 30, -22. / 15, 7. / 30}, {-22. / 15, 8. / 5, -2. / 15}, {7. / 30, -2. / 15, -1. / 10} };
   vector<vector<double>> G2 = { {6. / 5, -16. / 15, -2. / 15}, {-16. / 15, 32. / 15, -16. / 15}, {-2. / 15, -16. / 15, 6. / 5} };
   vector<vector<double>> G3 = { {-1. / 10, -2. / 15, 7. / 30}, {-2. / 15, 8. / 5, -22. / 15}, {7. / 30, -22. / 15, 37. / 30} };
   vector<vector<double>> M = { {2. / 15, 1. / 15, -1. / 30}, {1. / 15, 8. / 15, 1. / 15}, {-1. / 30, 1. / 15, 2. / 15} };
   double h = x[k + 1] - x[k];
   double xk = x[k];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         B[i][j] = (lambda(x[k]) * G1[i][j] + lambda(x[k] + x[k + 1] / 2.) * G2[i][j] + lambda(x[k + 1]) * G3[i][j]) / h;
         C[i][j] = M[i][j] * h;
      }
   vector<double> tmp(3);
   tmp[0] = func(x[k]);
   tmp[1] = func(((x[k] + x[k + 1]) / 2.));
   tmp[2] = func(x[k + 1]);
   for (int i = 0; i < 3; i++)
   {
      F[i] = 0;
      for (int j = 0; j < 3; j++)
         F[i] += C[i][j] * tmp[j];
   }
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
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
)
{
   //профиль
   ig[0] = 1;
   ig[1] = 1;
   for (int i = 2; i < m + 1; i++)
      ig[i] = ig[i - 1] + i - 1;
   gg.resize(ig[m] - 1);
   for (int i = 0; i < m; i++)
   {
      G[i] = 0.0;
      di[i] = 0.0;
   }

   for (int i = 1; i <= n; i++)
   {
      local_build2(n, i - 1, gamma, B, C, F, x);
      for (int j = 0; j < 4; j++)
      {
         di[2 * (i - 1) + j] += B[j][j] + gamma * C[j][j];
         G[2 * (i - 1) + j] += F[j];
      }
      gg[ig[3 * (i - 1) + 2] - 2] = B[1][0] + gamma * C[1][0];
      gg[ig[3 * (i - 1) + 2] - 1] = B[2][0] + gamma * C[2][0];
      gg[ig[3 * (i - 1) + 2]] = B[2][1] + gamma * C[2][1];
      gg[ig[3 * (i - 1) + 3] - 1] = B[3][0] + gamma * C[3][0];
      gg[ig[3 * (i - 1) + 3]] = B[3][1] + gamma * C[3][1];
      gg[ig[3 * (i - 1) + 3] + 1] = B[3][2] + gamma * C[3][2];
   }
}

void local_build2(
   int n,
   int k,
   double gamma,
   vector<vector<double>>& B,
   vector<vector<double>>& C,
   vector<double>& F,
   vector<double>& x
)
{
   vector<vector<double>> G1 = { {2237. / 840, -1947. / 560, 291. / 280, -379. / 1680}, {-1947. / 560, 1377. / 280, -1053. / 560, 123. / 280}, {291. / 280, -1053. / 560, 243. / 280, -3. / 112}, {-379. / 1680, 123. / 280, -3. / 112, -157. / 840} };
   vector<vector<double>> G2 = { {257. / 210, -171. / 140, -9. / 70, 53. / 420}, {-171. / 140, 351. / 70, -513. / 140, -9. / 70}, {-9. / 70, -513. / 140, 351. / 70, -171. / 140}, {53. / 420, -9. / 70, -171. / 140, 257. / 210} };
   vector<vector<double>> G3 = { {-157. / 840, -3. / 112, 123. / 280, -379. / 1680}, {-3. / 112, 243. / 280, -1053. / 560, 291. / 280}, {123. / 280, -1053. / 560, 1377. / 280, -1947. / 560}, {-379. / 1680, 291. / 280, -1947. / 560, 2237. / 840} };
   vector<vector<double>> M = { {8. / 105, 33. / 560, -3. / 140, 19. / 1680}, {33. / 560, 27. / 70, -27. / 560, -3. / 140}, {-3. / 140, -27. / 560, 27. / 70, 33. / 560}, {19. / 1680, -3. / 140, 33. / 560, 8. / 105} };
   double h = x[k + 1] - x[k];
   double xk = x[k];
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
         B[i][j] = (lambda(x[k]) * G1[i][j] + lambda(x[k] + x[k + 1] / 2.) * G2[i][j] + lambda(x[k + 1]) * G3[i][j]) / h;
         C[i][j] = M[i][j] * h;
      }
   vector<double> tmp(4);
   tmp[0] = func(x[k]);
   tmp[1] = func(((x[k] + x[k + 1])/ 3));
   tmp[2] = func(2 * ((x[k] + x[k + 1]) / 3));
   tmp[3] = func(x[k + 1]);
   for (int i = 0; i < 4; i++)
   {
      F[i] = 0;
      for (int j = 0; j < 4; j++)
         F[i] += C[i][j] * tmp[j];
   }
}

void llt_decompose2(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
)
{
   int i, j, k;
   double per;
   int a;
   int b;
   di[0] = sqrt(di[0]);
   for (i = 1; i < n; i++)
   {
      a = i - ig[i + 1] + ig[i]; 
      for (j = 0; j < ig[i + 1] - ig[i]; j++)
      {
         b = a + j - ig[a + j + 1] + ig[a + j]; 
         per = gg[ig[i] + j - 1];

         if (a < b)
         {
         for (k = ig[a + j + 1] - ig[a + j] - 1; k >= 0; k--)
         {
         per -= gg[ig[a + j] + k - 1] * gg[ig[i] + b - a + k - 1];
         }
         }
         else
         {
         for (k = a - b; k < ig[a + j + 1] - ig[a + j]; k++)
         {
         per -= gg[ig[i] + k - 1 - (a - b)] * gg[ig[a + j] + k - 1];
         }
         }
         gg[ig[i] + j - 1] = per / di[a + j];
      }
      per = di[i];
      for (k = 0; k < ig[i + 1] - ig[i]; k++)
         per -= gg[ig[i] + k - 1] * gg[ig[i] + k - 1];
      di[i] = sqrt(per);
   }
}

void gauss(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
)
{
   vector<double>q(n*4);
   int i, j;
   //прямой ход метода Гаусса
   q[0] = G[0] / di[0];
   for (i = 1; i < n; i++)
   {
      q[i] = G[i];
      for (j = 0; j < ig[i + 1] - ig[i]; j++)
      {
         q[i] -= gg[ig[i] + j - 1] * q[i - ig[i + 1] + ig[i] + j];
      }
      q[i] = q[i] / di[i];
   }
   //Заносим в вектор правой части решение обратного хода метода Гаусса
   for (i = 0; i < n; i++)
      G[i] = q[i];
   //обратный ход метода Гаусса
   for (i = n - 1; i >= 0; i--)
   {
      q[i] = G[i] / di[i];
      for (j = 0; j < ig[i + 1] - ig[i]; j++)
      {
         G[i - ig[i + 1] + ig[i] + j] -= q[i] * gg[ig[i] + j - 1];
      }
   }
}

void slau_llt2(
   int n, int m,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<double>& G
)
{
   llt_decompose2(n, G, gg, ig, di);
   gauss(n, G, gg, ig, di);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void llt_decompose(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
)
{
   di[0] = sqrt(di[0]);
   for (size_t i = 1; i <= n; i++)
   {
      gg[ig[2 * i - 1] - 1] /= di[2 * i - 2];
      di[2 * i - 1] = sqrt(di[2 * i - 1] - (gg[ig[2 * i - 1] - 1] * gg[ig[2 * i - 1] - 1]));
      gg[ig[2 * i] - 1] /= di[2 * i - 2];
      gg[ig[2 * i]] = (gg[ig[2 * i]] - gg[ig[2 * i - 1] - 1] * gg[ig[2 * i] - 1]) / di[2 * i - 1];
      di[2 * i] = sqrt(di[2 * i] - gg[ig[2 * i] - 1] * gg[ig[2 * i] - 1] - gg[ig[2 * i]] * gg[ig[2 * i]]);
   }
}

void forward_prop(
   int n,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di
)
{
   G[0] /= di[0];
   for (int i = 1; i <= n; i++)
   {
      G[2 * i - 1] = (G[2 * i - 1] - gg[ig[2 * i - 1] - 1] * G[2 * i - 2]) / di[2 * i - 1];
      G[2 * i] = (G[2 * i] - gg[ig[2 * i] - 1] * G[2 * i - 2] - gg[ig[2 * i]] * G[2 * i - 1]) / di[2 * i];
   }
}

void backward_prop(
   int n,
   int m,
   vector<double>& G,
   vector<double>& gg,
   vector<int>& ig,
   vector<double>& di)
{
   G[m - 1] /= di[m - 1];
   for (int i = n; i > 0; i--)
   {
      G[2 * i - 1] = (G[2 * i - 1] - gg[ig[2 * i]] * G[2 * i]) / di[2 * i - 1];
      G[2 * i - 2] = (G[2 * i - 2] - gg[ig[2 * i] - 1] * G[2 * i] - gg[ig[2 * i - 1] - 1] * G[2 * i - 1]) / di[2 * i - 2];
   }
}

void slau_llt(
   int n, int m,
   vector<int>& ig,
   vector<double>& di,
   vector<double>& gg,
   vector<double>& G
)
{
   //разложение
   llt_decompose(n, G, gg, ig, di);
   //прямой ход
   forward_prop(n, G, gg, ig, di);
   //обратный ход
   backward_prop(n, m, G, gg, ig, di);
}

void boundary(
   int n,
   int m,
   vector<double>& di,
   vector<double>& G,
   vector<double>& x
   )
{
   bool e1[2] = { 1,1 };
   bool e2[2] = { 0,0 };
   bool e3[2] = { 0,0 };
   double ubetta1 = 1.5;
   double ubetta2 = 5.0;
   double betta1 = 2.0;
   double betta2 = 2.0;
   double tetta1 = 2.0;
   double tetta2 = -1;
   if (e3[0] == 1)
   {
      di[0] += ubetta1 + betta1;
      G[0] += betta1 * u(x[0]);
   }
   if (e3[1] == 1)
   {
      di[m - 1] += ubetta2 + betta2;
      G[m - 1] += betta2 * u(x[n]);
   }
   if (e2[0] == 1)
   {
      G[0] += tetta1;
   }
   if (e2[1] == 1)
   {
      G[n - 1] += tetta2;
   }
   if (e1[0] == 1)
   {
      di[0] = 1.e+30;
      G[0] = 1.e+30 * u(x[0]);
   }
   if (e1[1] == 1)
   {
      di[m - 1] = 1.e+30;
      G[m - 1] = 1.e+30 * u(x[n]);
   }
}

void output(
   int n,
   vector<double>& x,
   vector<double>& G
   )
{
   ofstream output;
   output.open("out.txt");
   output << x[0] << " " << u(x[0]) << " " << G[0] << endl;
   for (int i = 1; i <= n; i++)
   {
      double x2 = (x[i - 1] + x[i]) / 2.0;
      output << x2 << " " << u(x2) << " " << G[2 * i - 1] << endl;
      output << x[i] << " " << u(x[i]) << " " << G[2 * i] << endl;
   }
}

void output2(
   int n,
   vector<double>& x,
   vector<double>& G
)
{
   ofstream output;
   output.open("out.txt");
   output << x[0] << " " << u(x[0]) << " " << G[0] << endl;
   for (int i = 1; i <= n; i++)
   {
      double x2 = (x[i - 1] + x[i]) / 3;
      double x3 = 2 * (x[i - 1] + x[i]) / 3;
      output << x2 << " " << u(x2) << " " << G[3 * i - 2] << endl;
      output << x3 << " " << u(x3) << " " << G[3 * i - 1] << endl;
      output << x[i] << " " << u(x[i]) << " " << G[3 * i] << endl;
   }
}