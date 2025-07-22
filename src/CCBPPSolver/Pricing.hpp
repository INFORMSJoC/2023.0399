#ifndef PRICING_H
#define PRICING_H

#include "CutPlaneManager.hpp"

using namespace std;

struct Pos{
  inline static vector<vll>* dp;
  int i;
  Pattern* p;
  vi seen;
  Pos(int i, Pattern* p, vi seen): i(i), p(p), seen(seen){}

  bool operator < (const Pos& p2) const{
    return val() == p2.val() ? (i > p2.i) : (val() < p2.val());
  }

  ll val() const{
    return (*dp)[i][p->residue()] + p->val;
  }
};

struct pll{
  ll v[2];
  
  pll(ll v1 = 0, ll v2 = 0){
    v[0] = v1;
    v[1] = v2;
  }

  ll& operator [] (int i){
    return v[i];
  }
};

class Pricing {
  private:
    const ll INFMAX = LONG_LONG_MAX / 2;
    set<Pattern*, Pattern::lessPatternEQ> pool;
    int W, trysDP, mxFreqRet = 3;
    vw &items;
    const vi &freq, &conflictIts;
    const vll &v, &vSR;
    const vector<vector<IM>>& itsToSR;
    const vector<Cut>& setsSR;
    const vector<vc>& compatible;
    vi& seenSR;
    vector<vector<vector<pll>>> dp;
    vi its, freqRet;
    vector<Pattern*> retDP;
    Pattern* incumb = NULL;
    bool goodAll, abortDP;
    vi seenColors;

    vector<Pattern*> getFromPool();
    void dryBests(vector<Pattern*>& bests, bool del);

    bool isValid(Pattern* pat, ll val2 = INFLL);
    void recoverBest(int i, int b, int c, bool u, Pattern* pat);
    void recoverDPLexic(int i, int b, int c, bool u, Pattern* pat);

    bool updateDPsize(int N, ll INFMAX);
    vector<Pattern*> DP(bool& opt, bool& feasible);

    bool canPack(Pattern* pat, int it);
    void addItVal(Pattern* pat, int it);
    void addIt(Pattern* pat, int it);
    void remIt(Pattern* pat, int it);

  public:
    int residueAccept, cnt = 0, cntDP = 0;
    ll generatedCols = 0;
    bool unique = false, &multiColor;

    Pricing(CPM* cpm, vi& freq, vi& conflictIts, vll& v, vector<vc>& compatible, bool& multiColor);
    ~Pricing();

    vector<Pattern*> Solve(bool& opt, bool& feasible);

    void AddPool(Pattern* bin, bool force = false);
    void RemPool(Pattern* bin, bool force = false);
    void ClearPool(int residue);

    bool IsCompatible(Pattern* pat);
    bool IsFeasible(Pattern* pat);
    bool IsColorFeasible(Pattern* pat);

    bool CanBeAdd(Pattern* pat);
    // ll MinReducedCost();
    ll UpdateVal(Pattern* pat);

    ll RecursiveDP(int i, int b, int c, bool u);
};

#endif
