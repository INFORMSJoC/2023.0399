#ifndef UTIL_H
#define UTIL_H

#include <bits/stdc++.h>
#include "gurobi_c++.h"

using namespace std;
using vi = vector<int>;
using vvi = vector<vi>;
using vc = vector<char>;
using vvc = vector<vc>;
using ii = pair<int, int>;
using ll = long long;
using ld = double;
using vd = vector<ld>;
using vll = vector<ll>;
using vvd = vector<vd>;
using si = set<int>;
#define all(a) a.begin(), a.end()
#define rall(a) a.rbegin(), a.rend()
#define SZ(v) ((int)v.size())
#define max(a, b) ((a) > (b) ? (a) : (b))
#define FOR(i, b)   for(int i = 0; i < (b); i++)
#define FORI(i, b)  for(int i = 1; i <= (b); i++)
#define REV(i, b)   for(int i = (b) - 1; i >= 0; i--)


const ll BASIS = 1LL << 49;
const ld BASISD = 1LL << 49;
const ll INFLL = LONG_LONG_MAX - BASIS - 10;
const ll eps10 = BASIS / 1e10;
const ll epsM = BASIS / (1LL << 38);
const ll eps13 = BASIS / 1e13;

#define DD(val) ((val) / BASISD)
#define MF 400.0
#define EPS8 (5e-9)
#define EPS9 (1e-9 + 1e-14)
#define EPS10 (1e-10)
#define FAC(a) (fabs((a) - round(a)))
#define FLOORLL(a) (((a) / epsM) * epsM)
#define CEILLL(a) ((((a) + epsM - 1) / epsM) * epsM)
#define ROUND(a) (((a) + BASIS - 1) / BASIS)
//#define MAX_CUT_ALLOWED 100

extern int binCapacity;
extern unordered_set<string> features;

namespace Util{
  template<class T>
  void Unique(vector<T>& v){
    sort(all(v));
    v.resize(unique(all(v))  - v.begin());
  }

  template<class T>
  void Append(vector<T>& a, const vector<T>& b){
    a.insert(a.end(), all(b));
  }

  vi mask(const vi& v, const vi& M);
}

#endif