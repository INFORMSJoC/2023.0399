#include "Solver.hpp"

Solver::Solver(Instance* instance): instance(instance), model(new Model(instance)), 
  bm(new BranchManager(model)), relaxAndFix(new RFHeuristic(model)),
  gtHeuristic(new GTHeuristic(model, bm, relaxAndFix)), bestSol(model->bestSol),
  L0(model->L0), bestInteger(model->bestInteger), L1(model->L1), L2(model->L2), 
  optimal(model->optimal), bounded(model->bounded){}

Solver::~Solver(){
  delete model;
  delete bm;
  delete relaxAndFix;
  delete gtHeuristic;
}

bool Solver::solveRoot(){
  cout << "Initial Heuristic Solution: " << bestInteger << " - LB: " << L0 << endl;

  if(bestInteger <= max(L1, numMachines)){
    cout << "Greedy solution is optimal" << endl;
    L2 = L1;
    return optimal = true;
  }

  const bool useIneq = features.count("Ineq");
  if(useIneq)
    model->EnableDualIneq();

  ll rootObj = model->Optimize(ROOT);

  if(optimal){
    cout << "Optimal solution founded while the root was being resolved" << endl;
    L2 = L1;
    return true;
  }

  L2 = ROUND(rootObj);
  cout << "Solved Root" << endl;
  model->Print(true);

  if(PRUNE(rootObj))
    return optimal = true;

  if(model->PH(rootObj))
    return true;

  rootObj = model->Optimize(OPTIMIZE);
  model->Print(true);

  if(PRUNE(rootObj))
    return true;

  cout << "Cutting at Root" << endl;
  if(model->AddCuts(rootObj))
    return optimal = true;
  model->Print(true);

  return model->PH(rootObj);
}

bool Solver::relax(ll& curObj, bool Print){
  while(!bm->FixST(curObj, Print)){
    if(model->AddCuts(curObj)){
      if(optimal)
        break;
      continue;
    }
    return false;
  }
  return true;
}

void Solver::branchAndPrice(){
  try{
    int toRelax = 0, toCRF = 0;
    cout << "Begin Branch and Price" << endl;
    model->ShrinkRA();
    ll curObj = model->GetObj();

    while(!model->ReachTL()){
      assert(bm->GetLvl() >= 0);
      int D = bm->GetLvl() < 30 ? 30 : 20;
      if(features.count("CRF") && toCRF >= D && gtHeuristic->Run(curObj))
        return;
      if(features.count("RF") && toRelax % 10 == 0 && relaxAndFix->Run(curObj))
        return;

      iterationsBP++;
      toRelax++;
      toCRF %= D;
      toCRF++;
      bm->BranchToLeft(curObj);

      if(optimal || relax(curObj, true))
        return;

      if(model->PH(curObj))
        if(optimal || (model->Optimize(OPTIMIZE), relax(curObj, false)))
          return;
    }
    assert(optimal == false);
    return;
  } catch(GRBException e) {
    cout << "Error in Solver\nError code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(const char* msg) {
    if(strcmp(msg, "TLE") == 0)
      return;
    cout << "Exception during optimization of the LP" << endl;
  }
  exit(1);
}

Solution* Solver::Graham(){
  vi freq = instance->freq;
  int numBins = min(numMachines, accumulate(all(freq), 0));
  multiset<Pattern*, Pattern::lessPatternEQ> pats;

  REV(i, SZ(freq)){
    FOR(f, freq[i]){
      if(SZ(pats) < numBins){
        auto pat = new Pattern(instance->weights[i]);
        pats.insert(pat);
        continue;
      }
      auto iter = prev(pats.end());
      auto pat = *iter;
      pats.erase(iter);
      pat->add(instance->weights[i]);
      pat->update();
      pats.insert(pat);
    }
  }

  auto sol = new Solution(instance);
  for(auto pat : pats)
    sol->InsertBin(pat);
  assert(sol->Cost() <= numMachines);
  
  int span = (*sol->GetBins().begin())->cap;
  assert(instance->B1 <= span);

  binCapacity = span;
  sol->Check();

  return sol;
}

void Solver::Solve(){
  Solution* incSol = Graham();

  int LB = instance->B1, UB = (*incSol->GetBins().begin())->cap;
  if(LB == UB)
    optimal = true;

  bool lifting = true;
  int nxtPow = 0;

  while(LB != UB){
    int MID = lifting ? min(UB - 1, LB + nxtPow) : (LB + UB) / 2;
    bm->CleanStack(ALL, false);
    model->Reset(MID);

    if(!optimal){
      delete gtHeuristic;
      gtHeuristic = new GTHeuristic(model, bm, relaxAndFix);

      if(!solveRoot())
        branchAndPrice();
    }
    if(!optimal)
      break;
    if(bestInteger <= numMachines){
      delete incSol;
      incSol = new Solution(*bestSol);
      int mxCap = (*bestSol->GetBins().begin())->cap;
      assert(mxCap <= MID);
      UB = mxCap;
      lifting = false;
    }
    else{
      LB = MID + 1;
      if(lifting)
        nxtPow = nxtPow == 0 ? 1 : nxtPow << 1;
    }
  }

  binCapacity = UB;
  swap(bestSol, incSol);
  delete incSol;

  StatusSolution statusSol;
  statusSol.status =  optimal ? "Optimal" : "Time_Limit";
  statusSol.Z = UB;
  statusSol.finalTime = model->UpdateTime();

  bestSol->PrintSolution(instance->pathOutSol, statusSol);
  if(!optimal)
    bm->Print(true);
  model->Print(true);

  cout << "BinCap: " << statusSol.Z  << " - Machines: " << bestSol->Cost() << " - Status: " << statusSol.status;
  cout << " - Time: " << statusSol.finalTime <<  " - genCol: " << model->GetGenCols();
  cout << " - genCuts: " << model->GetGenCuts() << " - mxActCuts: " << model->GetMaxActCuts() << endl;
  cout << "Pricing Time: " << model->GetPricingTime();
  cout << " - LP Time: " << model->GetLPTime() << endl;
  cout << endl;

  ofstream logStream;
  logStream.open(instance->prefix + ".log");
  logStream << statusSol.Z << "," << bestSol->Cost() << "," << statusSol.finalTime << "," << (optimal ? "OPT" : "TLE") << ",";
  logStream << UB - instance->B1 << "," << (UB == instance->B1 ? "B1" : "NT") << "," << model->GetGenCols() << "," << model->GetGenCuts() << "," << iterationsBP;
  logStream << endl;
  logStream.close();

  delete bestSol;
  return;
}
