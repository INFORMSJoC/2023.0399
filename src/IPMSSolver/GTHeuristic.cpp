#include "GTHeuristic.hpp"

const int KM = 12;

void GTHeuristic::addPoolSols(){
  for(Solution* sol : toAddGT){
    if(sol->Cost() > KM){
      set<int> patIds;
      for(auto orig : sol->GetBins()){
        Pattern* pat = new Pattern(*orig);
        Pattern::updateId(pat);
        patIds.insert(pat->id);
        model->AddPoolInactive(pat);
      }

      constrsGT.push_back(ConstrGT(patIds, sol->Cost()));
    }
    delete sol;
  }
  toAddGT.clear();
}

void GTHeuristic::disableConstrsGT(){
  for(auto& constr : constrsGT)
    if(constr.enable){
      model->RemConstr(constr.grbConstr);
      constr.enable = false;
    }

  model->Update();
}

void GTHeuristic::enableConstrsGT(){
  assert(bm->AllRight());

  REV(i, SZ(constrsGT)){
    bool good = true;
    si ids = constrsGT[i].ids;
    for(int id : ids){
      bool keep = Pattern::idToPat[id]->residue() <= model->GetRA() && model->CanBeAdd(Pattern::idToPat[id]);
      if(!keep){
        good = false;
        constrsGT[i].ids.erase(id);
      }
    }
    if(!good){
      if(constrsGT[i].ids.empty()){
        swap(constrsGT[i], constrsGT.back());
        constrsGT.pop_back();
      }
      else
        constrsGT[i].rhs = SZ(constrsGT[i].ids);
    }
  }

  assert(!constrsGT.empty());

  stable_sort(all(constrsGT));
  while(SZ(constrsGT) > 10)
    constrsGT.pop_back();

  auto& constr = constrsGT[0];
  vi allPats(all(constr.ids));
  model->AddPatternList(allPats, true);

  constr.enable = true;
  model->Update(true);

  constr.grbConstr = model->AddConstr(allPats, GRB_GREATER_EQUAL, constr.rhs - 6, "H" + to_string(0));

  model->Update();
}

bool GTHeuristic::reachTL(){
  int timeCurrent = model->UpdateTime();
  return (timeCurrent - initialTime > TIMELIMIT / 10) || timeCurrent >= TIMELIMIT;
}

bool GTHeuristic::Run(ll& curObj){
  static int cnt = 0;
  cnt++;
  addPoolSols();
  if(constrsGT.empty()){
    cout << "Skip GT Heuristic" << endl;
    return false;
  }

  cout << "Begin GT Heuristic" << endl;

  ll rootObj = bm->CleanStack(HARD, true);

  assert(model->bounded);

  enableConstrsGT();
  insideGT = true;
  initialTime = model->UpdateTime();
  curObj = model->Optimize(OPTIMIZE);

  while(!model->optimal && !reachTL() && !BETTER_EQUAL_THAN(curObj, rootObj)){
    assert(constrsGT[0].enable);
    constrsGT[0].rhs -= 2;
    constrsGT[0].grbConstr.set(GRB_DoubleAttr_RHS, constrsGT[0].rhs);
    curObj = model->Optimize(OPTIMIZE);
  }
  cout << "ConstGT RHS: " << constrsGT.front().rhs << "\n";

  if(model->optimal){
    cout << "Find Optimal Solution in GT before run Relax And Fix" << endl;
    return true;
  }

  relaxAndFix->Run(curObj);

  disableConstrsGT();

  constrsGT[0].cnt++;
  addPoolSols();
  if(constrsGT[0].cnt >= 5){
    swap(constrsGT[0], constrsGT.back());
    constrsGT.pop_back();
  }

  cout << "End GT Heuristic - Fails: " << constrsGT[0].cnt << "\n" << endl;
  insideGT = false;

  bm->RecoverST();
  curObj = model->Optimize(OPTIMIZE);

  return model->optimal;
}
