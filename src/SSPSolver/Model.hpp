#ifndef MODEL_H
#define MODEL_H

#include "Pattern.hpp"
#include "Pricing.hpp"
#include "DynamicDAG.hpp"
#include "RoundingHeuristic.hpp"
#include "ConflictDAG.hpp"

extern int TIMELIMIT;
extern const int K;

#define PRUNE(val) (ROUND(val) <= bestInteger)
#define BETTER_EQUAL_THAN(cur, old) (ROUND(cur) >= ROUND(old))

enum WasteRule {ROOT, OPTIMIZE, NONE };

struct ConstrGT{
  GRBConstr grbConstr;
  si ids;
  int rhs, cnt = 0;
  bool enable;
  ConstrGT(si ids, int rhs): ids(ids), rhs(rhs), enable(false){}

  bool operator < (const ConstrGT& c2) const{
    assert(!enable && !c2.enable);
    return ii(-rhs, SZ(ids)) < ii(-c2.rhs, SZ(c2.ids));
  }
};

class Model{
  private:
    Instance* instance;
    GRBModel* model;
    Pricing* pricing;
    CPM* cpm;
    RoundingHeuristic* roundingHeuristic;
    
    chrono::steady_clock::time_point begin;
    GRBVar* vars = NULL;
    GRBConstr* constrs = NULL;
    int lastPrint = 0, numDualIneq = 0, lastPrintBP = 0, numIssues = 0;
    ld pricingTime = 0, LPTime = 0, timeCurrent = 0;
    ll lastVal, obj, safeObj, lastObj;

    vi conflictIts;
    vector<GRBVar> dualInq;

    void initModel();

    void getPrimal();
    void getDual();
    bool colGen();
    
    void addPattern(Pattern* pat);
    void remPattern(Pattern* pat);

    ll lagrangianDualBound();

    void useInitialHeuristic();
    void storeSolution(Solution* sol);

    void rootRuleOptim();
    bool costOptim(ll curObj);

    void computeInvIds();
    void updateConstrs();
    void updateConflictIts();
    void updateIntegers();

    void checkModelValid();
    
  public:
    Solution* bestSol, *partial;
    ld L0;
    int L1, L2, N, NC, M, bestInteger, integers, nonZero;
    bool bounded = false, optimal = false, masterCol = false, insideRF = false, insideGT = false;
    bool rfHeuristic = false, optPricing = false;

    vi weights, freq, ids, wToIt;
    vll dual;
    vd primal;
    vector<vd> afim;
    vvc compatible;

    vector<Pattern*>& idToPat;
    vector<set<int>> itToIds;
    map<int, int> invIds;
    multimap<int, ii> mergedIts;
    vector<ConstrGT> constrsGT;
    vector<Solution*> toAddGT;
    set<Pattern*, Pattern::lessPatternEQ> actPatterns;

    Model(Instance* instance);
    ~Model();

    ll GetObj();
    ld UpdateTime();
    void Print(bool force = false, vi st = vi());
    void Update(bool updateInv = false);

    GRBConstr AddConstr(vi patIds, char sense, int rhs, string name);
    void RemConstr(GRBConstr constr);
    
    int GetSzMergeIts();
    int GetNumVars();
    int GetIntegers();
    int GetNumIssues(){ return numIssues; }

    int GetNumCuts();
    int GetGenCuts();
    int GetMaxActCuts();

    int GetCA();
    int GetGenCols();
    ld GetLPTime();
    ld GetPricingTime();
    
    bool ReachTL();
    void ShrinkRA();
    
    vi GetFreqInst();
    vi GetFreqCur();
    
    bool AddCuts(ll& curObj);
    ll Optimize(WasteRule rule, ll oldObj = INFLL);
    ll RunPL();
    void FastPL();

    void CalcAfim();
    vi MergeIts(int it1, int it2, int it, vi idsToAdd);
    vi SplitIt(int it1, int it2, int it, vi idsToAdd);
    
    void AddConstrGT(Solution* sol);
    bool PH(ll lb);

    void AddIts(const vi& its, vi pats, bool force = false, bool improveIncumb = false, bool gt = false);
    vi RemoveIts(vi toRem);
    void AddPatternList(vi& l, bool force = false, bool improveIncumb = false, bool gt = false);
    bool RemPatternList(vi& l);
    
    int AppendIt(int wn);
    void AddPoolInactive(Pattern* pat);

    void EnableDualIneq();
    bool ResidueOptim(WasteRule rule);

    ii GetNumPricingCalls();

    bool CanBeAdd(Pattern* pat);
};

#endif
