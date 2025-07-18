#include "Solver.hpp"

Solver::Solver(Instance* instance): instance(instance), model(new Model(instance)), 
  bm(new BranchManager(model)), relaxAndFix(new RFHeuristic(model)), bestSol(model->bestSol),
  L0(model->L0), bestInteger(model->bestInteger), L1(model->L1), L2(model->L2), 
  optimal(model->optimal), bounded(model->bounded){}

Solver::~Solver(){
  delete model;
  delete bm;
  delete relaxAndFix;
}

bool Solver::solveRoot(){
  cout << "Initial Heuristic Solution: " << bestInteger << " - LB: " << L0 << endl;

  if(bestInteger <= L1){
    cout << "Greedy solution is optimal" << endl;
    statusSol.rootObj = statusSol.L2 = L2 = L1;
    statusSol.integers = statusSol.nonZero = 1;
    statusSol.rootTime = model->UpdateTime();
    statusSol.where = 0;
    return optimal = true;
  }

  const bool useIneq = features.count("Ineq") && C > 1;
  if(useIneq)
    model->EnableDualIneq();

  ll rootObj = model->Optimize(ROOT);
  statusSol.rootObj = DD(rootObj);
  statusSol.integers = 1.0L * model->integers / ROUND(rootObj);
  statusSol.nonZero = 1.0L * model->nonZero / ROUND(rootObj);

  if(optimal){
    cout << "Optimal solution founded while the root was being resolved" << endl;
    statusSol.rootObj = statusSol.L2 = L2 = L1;
    statusSol.rootTime = model->UpdateTime();
    return true;
  }

  statusSol.L2 = L2 = ROUND(rootObj);
  statusSol.rootTime = model->UpdateTime();
  cout << "Solved Root" << endl;
  model->Print(true);

  if(PRUNE(rootObj))
    return optimal = true;

  /*if(model->PH(rootObj))
    return true;

  rootObj = model->Optimize(OPTIMIZE);
  model->Print(true);

  if(PRUNE(rootObj))
    return true;

  cout << "Cutting at Root" << endl;
  if(model->AddCuts(rootObj))
    return optimal = true;
  statusSol.integers = 1.0L * model->integers / ROUND(rootObj);
  statusSol.nonZero = 1.0L * model->nonZero / ROUND(rootObj);

  model->Print(true);*/

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
    int toRelax = 0;
    cout << "Begin Branch and Price" << endl;
    ll curObj = model->GetObj();

    while(!model->ReachTL()){
      assert(bm->GetLvl() >= 0);
      if(features.count("RF") && toRelax % 10 == 0 && relaxAndFix->Run(curObj))
        return;

      iterationsBP++;
      toRelax++;
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

void Solver::Solve(){
  try{
    if(!solveRoot()){
      branchAndPrice();
      statusSol.where = optimal ? 2 : -1;
    }
    else
      statusSol.where = 1;
  }catch(const char* msg){
    if(strcmp(msg, "TLE") != 0)
      throw msg;
    optimal = false;
    statusSol.rootTime = model->UpdateTime();
    statusSol.L2 = statusSol.rootObj = L1;
  }

  statusSol.status =  optimal ? "Optimal" : "Time_Limit";
  statusSol.Z = bestSol->Cost();
  statusSol.finalTime = model->UpdateTime();

  bestSol->PrintSolution(instance->pathOutSol, statusSol);
  if(!optimal)
    bm->Print(true);
  model->Print(true);

  cout << "IntegerSol: " << bestSol->Cost() << " - Status: " << statusSol.status;
  cout << " - Time: " << statusSol.finalTime <<  " - genCol: " << model->GetGenCols();
  cout << " - genCuts: " << model->GetGenCuts() << " - mxActCuts: " << model->GetMaxActCuts();
  cout << " - Int: " << statusSol.integers << " - NZ: " << statusSol.nonZero<< "\n";
  cout << " - MultiColors: " << model->wereThereMultiColorsItem() << "\n";
  
  cout << "Pricing Time: " << model->GetPricingTime();
  cout << " - LP Time: " << model->GetLPTime();
  cout << " - NumIssues: " << model->GetNumIssues() << endl;
  cout << endl << endl;

  ofstream logStream;
  logStream.open(instance->prefix + ".log");
  logStream << statusSol.Z << "," << statusSol.rootObj << ",";
  logStream << statusSol.finalTime << "," << (optimal ? "OPT" : "TLE") << "," << model->GetGenCols();
  logStream << "," << model->GetGenCuts() << "," << iterationsBP;
  logStream << "," << (optimal ? (L2 == statusSol.Z ? "AI" : "ANI") : "?");
  logStream << "," << model->GetMaxActCuts();
  logStream << "," << statusSol.integers << "," << statusSol.nonZero;
  logStream << "," << model->GetPricingTime() << "," << model->GetLPTime();
  logStream << "," << model->GetNumIssues() << "," << statusSol.where;
  logStream << "," << bm->GetNumBranches() << "," << model->GetNumPricingCalls().first;
  logStream << "," << model->GetNumPricingCalls().second;
  logStream << "," << model->wereThereMultiColorsItem();
  logStream << endl;
  logStream.close();

  delete bestSol;
  return;
}
