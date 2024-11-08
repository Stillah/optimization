#include <bits/stdc++.h>
#include "LinearAlgebra.h"
using namespace std;

bool balanced(const ColumnVector<double> &S, const ColumnVector<double> &D) {
  double supply = 0, demand = 0;
  for (int i = 0; i < S.getM(); i++) supply += S[i];
  for (int i = 0; i < D.getM(); i++) demand += D[i];
  return supply == demand;
}

void printTable(const Matrix<double> &C, const ColumnVector<double> &S, const ColumnVector<double> &D) {
  const int m = C.getM(), n = C.getN();

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      cout << C[i][j] << " ";
    }
    cout << "| " << S[i] << "\n";
  }
  cout << "---------------\n";
  for (int i = 0; i < n; i++)
    cout << D[i] << " ";
  cout << "\n\n";
}


void NorthWest(const Matrix<double> &C, ColumnVector<double> S, ColumnVector<double> D) {
  cout << "North-West corner:\n";
  if (!balanced(S, D)) { cout << "The problem is not balanced\n"; return; }
  const int m = C.getM(), n = C.getN();
  int x = 0, y = 0;
  Matrix<double> X(m, n);
  while (y != n) {
    double mn = min(S[x], D[y]);
    X[x][y] = mn;
    S[x] -= mn;
    D[y] -= mn;
    while (x < m && S[x] == 0) x++;
    while (y < n && D[y] == 0) y++;
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) cout << X[i][j] << " ";
    cout << "\n";
  }
  cout << "\n";
}

void Vogel(const Matrix<double> &C, ColumnVector<double> S, ColumnVector<double> D) {
  cout << "Vogel approximation:\n";
  if (!balanced(S, D)) { cout << "The problem is not balanced\n"; return; }
  const int m = C.getM(), n = C.getN();
  Matrix<double> X(m, n);

  set<pair<int, bool>> removed;
  while (true) {
    vector<pair<double, int>> diffRows, diffCols;
    vector<priority_queue<double>> minCol(n), minRow(m);
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (removed.contains({i, true}) || removed.contains({j, false})) continue;
        minRow[i].push(C[i][j]);
        minCol[j].push(C[i][j]);

        if (minRow[i].size() > 2) minRow[i].pop();
        if (minCol[j].size() > 2) minCol[j].pop();
      }
    }

    for (int i = 0; i < minRow.size(); i++) {
      if (minRow[i].empty()) continue;
      double a = minRow[i].top();
      if (minRow[i].size() > 1) minRow[i].pop();
      double b = minRow[i].top();
      diffRows.emplace_back(a - b, i);
    }
    for (int i = 0; i < minCol.size(); i++) {
      if (minCol[i].empty()) continue;
      double a = minCol[i].top();
      if (minCol[i].size() > 1) minCol[i].pop();
      double b = minCol[i].top();
      diffCols.emplace_back(a - b, i);
    }

    if (diffCols.empty() && diffRows.empty()) break;
    sort(diffRows.rbegin(), diffRows.rend());
    sort(diffCols.rbegin(), diffCols.rend());

    if (diffCols.empty() || diffRows[0] >= diffCols[0]) {
      double mn = 1e9;
      int j = 0;
      for (int i = 0; i < n; i++)
        if (mn > C[diffRows[0].second][i] && !removed.contains({i, false}))
          mn = C[diffRows[0].second][i], j = i;

      mn = min(D[j], S[diffRows[0].second]);
      D[j] -= mn;
      S[diffRows[0].second] -= mn;
      X[diffRows[0].second][j] += mn;
      if (D[j] == 0) removed.insert({j, false});
      if (S[diffRows[0].second] == 0) removed.insert({diffRows[0].second, true});
    }

    else {
      double mn = 1e9;
      int j = 0;
      for (int i = 0; i < m; i++)
        if (mn > C[i][diffCols[0].second] && !removed.contains({i, true}))
          mn = C[i][diffCols[0].second], j = i;

      mn = min(D[diffCols[0].second], S[j]);
      D[diffCols[0].second] -= mn;
      S[j] -= mn;
      X[j][diffCols[0].second] += mn;
      if (D[diffCols[0].second] == 0) removed.insert({diffCols[0].second, false});
      if (S[j] == 0) removed.insert({j, true});
    }
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) cout << X[i][j] << " ";
    cout << "\n";
  }
  cout << "\n";
}

void Russel(const Matrix<double> &C, ColumnVector<double> S, ColumnVector<double> D) {
  cout << "Russel approximation:\n";
  if (!balanced(S, D)) { cout << "The problem is not balanced\n"; return; }
  const int m = C.getM(), n = C.getN();
  Matrix<double> X(m, n);

  set<pair<int, bool>> removed;
  vector<double> u(n), v(m);
  while (true) {
    map<double, pair<int, int>> delta;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (removed.contains({i, true}) || removed.contains({j, false})) continue;
        u[i] = max(u[i], C[i][j]);
        v[j] = max(v[j], C[i][j]);
      }
    }

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (removed.contains({i, true}) || removed.contains({j, false})) continue;
        delta[C[i][j] - u[i] - v[j]] = {i, j};
      }
    }

    if (delta.empty()) break;

    int i = delta.begin()->second.first, j = delta.begin()->second.second;
    double mn = min(S[i], D[j]);
    X[i][j] = mn;
    S[i] -= mn;
    D[j] -= mn;
    if (S[i] == 0) removed.insert({i, true});
    if (D[j] == 0) removed.insert({j, false});
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) cout << X[i][j] << " ";
    cout << "\n";
  }
  cout << "\n";
}

int main() {
  const int m = 3, n = 4;
  Matrix<double> C(m, n);
  ColumnVector<double> S(m);
  ColumnVector<double> D(n);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      cin >> C[i][j];
  for (int i = 0; i < m; i++)
    cin >> S[i];
  for (int i = 0; i < n; i++)
    cin >> D[i];

  try {
    printTable(C, S, D);
    NorthWest(C, S, D);
    Vogel(C, S, D);
    Russel(C, S, D);
  } catch (exception &e) { cout << e.what() << "\n"; }
  return 0;
}
