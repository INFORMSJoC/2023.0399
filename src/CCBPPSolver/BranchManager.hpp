#ifndef BRANCH_MANAGER_H
#define BRANCH_MANAGER_H

#include "Model.hpp"

struct Node {
  int it, it1, it2;
  vi join, split, backtrack;
  bool state, fixed;

  Node(int it1, int it2): it(-1), it1(it1), it2(it2), state(false), fixed(false){}
};

class BranchManager{
  private:
    Model* model;
    const int& N, &bestInteger;
    int lvl = 0, numBranches = 0;
    bool& bounded, &optimal;
    const vw &items;
    const vi &freq, &ids;
    vector<vd>& afim;
    vvc& compatible;
    const vector<set<int>> &itToIds;
    const vector<Pattern*>& idToPat;
    ConflictDAG<Node> conflictDAG;
    vector<Node> st, oldST;
    vector<ii> globalPrior;

    bool mergeIts(int it1 = -1, int it2 = -1);
    Pattern* mergeItPat(Pattern* pat, int it1, int it2);

    void idsMergedPats(int it1, int it2, vi& idsPats);
    void idsSplitPats(int it1, int it2, vi& idsPats);
    void findSplitAndBT(const vi& from, vi& split, vi& bk, int it1, int it2);

    vi reCalcAndPropag();
    vi conflictsPropagation(const vector<ii>& diff);
    void getConflictList(vi& co, int it1, int it2);

    void swapState();
    int backtracking();
    void remPosStack(int p);

  public:
    BranchManager(Model* model): model(model), N(model->N), bestInteger(model->bestInteger),
    bounded(model->bounded), optimal(model->optimal), items(Pattern::items), freq(model->freq),
    ids(model->ids), afim(model->afim), compatible(model->compatible),
    itToIds(model->itToIds), idToPat(model->idToPat), conflictDAG(compatible, st){}
    ~BranchManager(){}
    void BranchToLeft(ll& curObj);
    bool FixST(ll& curObj, bool print);
    
    ll CleanStack(bool hardReset, bool run);
    void RecoverST();
    void Print(bool force = false);

    bool AllRight();
    int GetLvl();
    int GetNumBranches();
};

#endif
