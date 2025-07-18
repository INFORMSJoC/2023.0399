#ifndef GT_HEURISTIC
#define GT_HEURISTIC

#include "BranchManager.hpp"
#include "RFHeuristic.hpp"

class GTHeuristic{
  private:
  Model* model;
  BranchManager* bm;
  RFHeuristic* relaxAndFix;
  vector<ConstrGT>& constrsGT;
  vector<Solution*>& toAddGT;
  map<int, int>& invIds;
  bool& insideGT;
  ld initialTime;

  void addPoolSols();
  void disableConstrsGT();
  void enableConstrsGT();

  bool reachTL();

  public:
  bool Run(ll& curObj);

  GTHeuristic(Model* model, BranchManager* bm, RFHeuristic* relaxAndFix): model(model),
  bm(bm), relaxAndFix(relaxAndFix), constrsGT(model->constrsGT), toAddGT(model->toAddGT),
  invIds(model->invIds), insideGT(model->insideGT){}
  ~GTHeuristic(){}
};

#endif