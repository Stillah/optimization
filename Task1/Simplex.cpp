#include <bits/stdc++.h>
#include "LinearAlgebra.h"
using namespace std;
/*
Function_name(C, A, b, eps = eps_default)
Input:
- C: A vector of coefficients of the objective function
- A: A matrix of coefficients of the constraint functions
- b: A vector of right-hand side values
- eps: Approximation accuracy (optional, default = eps_default)

Steps:
1. Print the optimization problem:
   - max (or min) z = C[0] * x1 + C[1] * x2 + ... + C[n] * xn
   - subject to the constraints:
     - A[0] * x <= b[0]
     - A[1] * x <= b[1]
     - ...
     - A[m] * x <= b[m]

2. Initialize:
   - Form the initial tableau by introducing slack variables to convert inequalities into equalities.

3. Iteratively apply the Simplex method:
   - Step 1: Identify the entering variable (most negative coefficient in the objective row).
   - Step 2: Identify the leaving variable (smallest positive ratio of RHS to pivot column).
   - Step 3: Perform pivot operations to update the tableau.

4. Check for optimality or unboundedness:
   - If all coefficients in the objective function row are non-negative, the solution is optimal.
   - If no leaving variable exists, the problem is unbounded.

5. Return:
   - solver_state: {solved, unbounded}
   - x*: Optimal vector of decision variables (if solved)
   - z: Maximum (or minimum) value of the objective function (if solved)
 */


pair<long double, ColumnVector<long double>>
simplex(const ColumnVector<long double> &C, const Matrix<long double> &A,
        const ColumnVector<long double> &b, bool maximization, const long double eps = 1e-9) {
  // Printing z function
  cout << (maximization ? "max " : "min ") << "z = ";
  for (int i = 0; i < C.getM(); i++) {
    cout << "(" << C[i] << " * x" << i + 1 << ") ";
    if (i != C.getM() - 1) cout << "+ ";
  }
  cout << endl << "subject to the constraints:" << endl;

  for (int i = 0; i < A.getM(); i++) {
    for (int j = 0; j < A.getN(); j++) {
      cout << "(" << A[i][j] << " * x" << j + 1 << ") ";
      if (j != A.getN() - 1) cout << "+ ";
    }
    cout << "<= " << b[i] << endl;
  }
  cout << endl;

  Matrix<long double> tableau(A.getM() + 1, A.getM() + A.getN() + 1);
  for (int i = 0; i < C.getM(); i++) tableau[0][i] = maximization? -C[i]: C[i];
  for (int i = 0; i < b.getM(); i++) tableau[i + 1].back() = b[i];
  for (int i = 0; i < A.getM(); i++) {
    for (int j = 0; j < A.getN(); j++)
      tableau[i + 1][j] = A[i][j];
    // slack variable
    tableau[i + 1][A.getN() + i] = 1;
  }


  vector<int> nonBasicVariables(A.getM() + A.getN() + 1);
  vector<int> basicVariables(A.getM() + 1);
  iota(nonBasicVariables.begin(), nonBasicVariables.end(), 0);
  iota(basicVariables.begin() + 1, basicVariables.end(), A.getN());
  // starting to iterate
  long double solution = 0;
  while (true) {
    long double minZ = 2e9;
    int indCol = -1;

    for (int i = 0; i < tableau.getN() - 1; i++) {
      if (minZ > tableau[0][i]) {
        minZ = tableau[0][i];
        indCol = i;
      }
    }
    // No variables to improve, optimal solution reached
    if (minZ >= -eps) break;

    long double minRation = 2e9;
    int indRow = -1;

    for (int i = 0; i < tableau.getM(); i++) {
      if (tableau[i][indCol] <= eps) continue;
      if (minRation > 1.0 * tableau[i].back() / tableau[i][indCol]) {
        minRation = 1.0 * tableau[i].back() / tableau[i][indCol];
        indRow = i;
      }
    }
    if (indRow == -1) throw invalid_argument("Unbounded solution\n");

    // some variable enters, some variable leaves
    swap(nonBasicVariables[indCol], basicVariables[indRow]);

    for (int i = 0; i < tableau.getM(); i++) {
      if (i == indRow) continue;
      long double multiplier = tableau[i][indCol] / tableau[indRow][indCol];
      for (int j = 0; j < tableau.getN(); j++)
        tableau[i][j] -= tableau[indRow][j] * multiplier;
    }
    long double divideBy = tableau[indRow][indCol];
    for (int j = 0; j < tableau.getN(); j++)
      tableau[indRow][j] /= divideBy;

    if (abs(solution - tableau[0].back()) < eps) break;
    solution = tableau[0].back();
  }
  ColumnVector<long double> X(A.getN());
  for (int i = 1; i < basicVariables.size(); i++)
    if (basicVariables[i] < X.getM()) X[basicVariables[i]] = tableau[i].back();

  return {solution, X};
}

int main() {
  int n, m;
  cin >> n >> m;
  ColumnVector<long double> C(n);
  Matrix<long double> A(m, n);
  ColumnVector<long double> b(m);

  for (int i = 0; i < n; i++)
    cin >> C[i];

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      cin >> A[i][j];
    cin >> b[i];
  }
  try {
    auto ans = simplex(C, A, b, true);
    cout << ans.first << "\n";
    for (int i = 0; i < ans.second.getM(); i++)
      cout << ans.second[i] << " ";
    cout << "\n";
  }
  catch (exception &e) {
    cout << e.what() << '\n';
  }
  return 0;
}
