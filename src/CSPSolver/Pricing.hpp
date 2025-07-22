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

class Pricing {
  private:
    set<Pattern*, Pattern::lessPatternEQ> pool;
    int W, mxFreqRet = 3;
    ll trysDP;
    const vi &weights, &freq, &wToIt, &conflictIts;
    const vll &v, &vSR;
    const vector<vector<IM>>& itsToSR;
    const vector<Cut>& setsSR;
    const vector<vc>& compatible;
    vi& seenSR;
    vector<vll> dp;
    vi its, freqRet;
    vector<Pattern*> retDP;
    Pattern* incumb = NULL;
    bool goodAll, abortDP;

    vector<Pattern*> getFromPool();
    void dryBests(vector<Pattern*>& bests, bool del);

    bool isValid(Pattern* pat, ll val2 = INFLL);
    void recoverBest(int i, Pattern* pat);
    void recoverDPLexic(int i, Pattern* pat);

    bool updateDPsize(int N, ll INFMAX);
    vector<Pattern*> DP(bool& opt, bool& feasible);

    bool canPack(Pattern* pat, int it);
    void addItVal(Pattern* pat, int it);
    void addIt(Pattern* pat, int it);
    void remIt(Pattern* pat, int it);

  public:
    int residueAccept, cnt = 0, cntDP = 0;
    ll generatedCols = 0;
    bool unique = false;

    Pricing(CPM* cpm, vi& weights, vi& freq, vi& wToIt, vi& conflictIts, vll& v, vector<vc>& compatible);
    ~Pricing();

    vector<Pattern*> Solve(bool& opt, bool& feasible);

    void AddPool(Pattern* bin, bool force = false);
    void RemPool(Pattern* bin, bool force = false);
    void ClearPool(int residue);

    bool IsCompatible(Pattern* pat);
    bool IsFeasible(Pattern* pat);

    bool CanBeAdd(Pattern* pat);
    ll MinReducedCost();
    ll UpdateVal(Pattern* pat);
};

#endif
