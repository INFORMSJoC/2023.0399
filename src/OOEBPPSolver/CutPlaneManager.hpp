#ifndef CPM_H
#define CPM_H

#include "Pattern.hpp"

struct IM{
  int i, m;
  bool operator <(const IM& im2) const{
    return i != im2.i ? i < im2.i : m < im2.m;
  }
};

struct Cut{
  vector<IM> its;
  int RHS, d, id = -1;

  Cut(vector<IM> its, int RHS, int d): its(its), RHS(RHS), d(d){}

  bool operator < (const Cut& c) const{
    assert(is_sorted(all(its)));
    assert(is_sorted(all(c.its)));
    return d != c.d ? d < c.d : its < c.its;
  }

  typedef vector<IM>::iterator iterator;
  typedef vector<IM>::const_iterator const_iterator;

  iterator begin(){
    return its.begin();
  }

  iterator end(){
    return its.end();
  }

  const_iterator begin() const{
    return its.cbegin();
  }

  const_iterator end() const{
    return its.cend();
  }
};

class CPM{
  GRBModel*& model;
  GRBVar*& vars;
  const int& N;
  bool& bounded;
  const vw &items;
  const vi &freq, &ids;
  const vd &primal;
  const vll &dual;
  const vector<vd>& afim;
  int actCuts = 0, mxActCuts = 0;
  const vector<set<int>> &itToIds;
  map<int, int>& invIds;
  vector<vector<IM>> itsToSRAll;
  vector<Cut> setsSRAll;
  vi countSR, mapSR;
  set<Cut> cutsPool;

  void updateItemList();
  vector<IM> getIdsCut(Cut cut);
  pair<ld, vector<IM>> calcCut(Cut cut);
  multimap<ld, pair<Cut, vector<IM>>, greater<ld>> calcCuts();
  multimap<ld, pair<Cut, vector<IM>>, greater<ld>> calcCutsUp();
  void rank4(multimap<ld, pair<Cut, vector<IM>>, greater<ld>>& cuts, vi& itsCut, set<vi>& cutsVI, int& cnt);
  void rank5(multimap<ld, pair<Cut, vector<IM>>, greater<ld>>& cuts, vi& itsCut, set<vi>& cutsVI, int& cnt);
  int addCutPool(Cut cut);

  public:
  vector<vector<IM>> itsToSR;
  vector<Cut> setsSR;
  vi seenSR;
  vll dualSR;

  CPM(GRBModel*& model, GRBVar*& vars, const int& N, bool& bounded, const vi& freq,
   const vi& ids, const vd& primal, const vll& dual, const vector<vd>& afim,
   const vector<set<int>> &itToIds, map<int, int>& invIds);

  void AddSRCoef(GRBColumn& col, Pattern* pat);
  void AddIt();

  bool GenCuts();
  int SetDualSR(ll val, int it);
  void ResetSR();
  void CheckCuts();

  int GetActCuts();
  int GetCuts();
  int GetMaxActCuts();
  void UpdMaxActCuts();
};

#endif
