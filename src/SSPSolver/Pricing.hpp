#ifndef PRICING_H
#define PRICING_H

#include "CutPlaneManager.hpp"

using namespace std;


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
    bool& insideGT;

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
    int capAccepted, cnt = 0, cntDP = 0;
    ll generatedCols = 0;
    bool unique = false;
    si gtSet;

    Pricing(CPM* cpm, vi& weights, vi& freq, vi& wToIt, vi& conflictIts, vll& v, vector<vc>& compatible, bool& insideGT);
    ~Pricing();

    vector<Pattern*> Solve(bool& opt, bool& feasible);

    void AddPool(Pattern* bin, bool force = false);
    void RemPool(Pattern* bin, bool force = false);
    void ClearPool(int residue);

    bool IsCompatible(Pattern* pat);
    bool IsFeasible(Pattern* pat);

    bool CanBeAdd(Pattern* pat);
    ll UpdateVal(Pattern* pat);
};

#endif
