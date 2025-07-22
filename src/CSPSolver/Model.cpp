#include "Model.hpp"

Model::Model(Instance* instance): instance(instance),
  lastObj(INFLL), L0(instance->L0), L1(instance->L1), L2(0), N(instance->size),
  weights(instance->weights), freq(instance->freq), wToIt(binCapacity + 1, -1),
  dual(N, 2 * BASIS), compatible(N, vc(N, 1)), idToPat(Pattern::idToPat), itToIds(N){

  cout << fixed << setprecision(12);
  FOR(i, N)
    wToIt[weights[i]] = i;

  cpm = new CPM(model, vars, N, bounded, weights, freq, wToIt, ids, primal, dual, afim, itToIds, invIds);
  pricing = new Pricing(cpm, weights, freq, wToIt, conflictIts, dual, compatible);
  roundingHeuristic = new RoundingHeuristic(M, weights, freq, ids, wToIt, idToPat, primal);

  partial = new Solution(instance);
  initModel();
  begin = chrono::steady_clock::now();
}

Model::~Model(){
  delete model;
  delete pricing;
  delete cpm;
  delete roundingHeuristic;
  delete partial;

  for(auto sol : toAddGT)
    delete sol;

  if(vars != NULL)
    delete[] vars;

  if(constrs != NULL)
    delete[] constrs;

  FOR(it, SZ(weights)){
    for(int id : itToIds[it]){
      Pattern* pat = idToPat[id];
      for(int w2 : pat->its){
        int it2 = wToIt[w2];
        if(it2 != it)
          itToIds[it2].erase(id);
      }
      delete pat;
    }
    itToIds[it].clear();
  }
  idToPat.clear();
  actPatterns.clear();
  itToIds.clear();
}

void Model::initModel(){
  try{
    string logPath = instance->prefix + ".gb";
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", logPath);
    env.set(GRB_IntParam_Threads, 1);
    env.set(GRB_DoubleParam_TimeLimit, TIMELIMIT);
    env.set(GRB_IntParam_Presolve, 0);
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Method, 0);
    env.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    env.set(GRB_DoubleParam_OptimalityTol, 1e-9);
    env.set(GRB_IntParam_NumericFocus, 2);
    env.start();

    // Create an empty model
    model = new GRBModel(env);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_ScaleFlag, 0);

    // Add constraints of covering
    FOR(i, N)
      model->addConstr(GRBLinExpr(), GRB_GREATER_EQUAL, freq[i], "u" + to_string(i));
    model->update();

    useInitialHeuristic();
    RunPL();
    assert(bounded);
  } catch(GRBException e) {
    cout << "Error in Init Model\nError code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
}

void Model::getPrimal(){
  M = model->get(GRB_IntAttr_NumVars);
  assert(M == SZ(actPatterns) + numDualIneq);
  primal.resize(M);
  ids.resize(M);

  if(vars != NULL)
    delete[] vars;
  vars = model->getVars();
  FOR(i, M){
    primal[i] = bounded ? vars[i].get(GRB_DoubleAttr_X) : 0.0;
    string name = vars[i].get(GRB_StringAttr_VarName);
    ids[i] = name[0] == 'x' ? atoi(&name[1]) : -1;
    if(ids[i] != -1)
      assert(actPatterns.count(idToPat[ids[i]]));
  }
}

void Model::getDual(){
  obj = partial->Cost() * BASIS;

  updateConstrs();

  fill(all(dual), -1);
  cpm->ResetSR();
  FOR(i, NC){
    ll pi = (constrs[i].get(GRB_DoubleAttr_Pi) / MF) * 1.0L * BASIS;
    string name = constrs[i].get(GRB_StringAttr_ConstrName);
    int it = atoi(&name[1]);
    if(name[0] == 'u'){
      dual[it] = pi;
      obj += freq[it] * dual[it];
    }else if(name[0] == 'c'){
      //assert(pi <= 0);
      int rhs = cpm->SetDualSR(pi, it);
      obj += rhs * pi;
    } else if(name[0] == 'H'){
      int RHS = constrs[i].get(GRB_DoubleAttr_RHS);
      obj += RHS * pi;
    }else
      assert(false);
  }
  cpm->UpdMaxActCuts();

  safeObj = ((__int128) obj * BASIS) / (BASIS + epsM)  + (((__int128) obj * BASIS) % (BASIS + epsM) == 0 ? 0 : 1);
}

bool Model::colGen(){
  //cout << "Doing pricing: " << obj << endl;
  bool opt = false;

  ld bef = UpdateTime();
  vector<Pattern*> pats = pricing->Solve(opt, optPricing);
  pricingTime += UpdateTime() - bef;

  if(!opt){
    assert(!pats.empty() && pats[0] != NULL);
    for(Pattern* pat : pats)
      addPattern(pat);
  }
  else{
    if(bounded)
      safeObj = max(safeObj, FarleyDualBound());
    optPricing = true;
  }
  return opt;
}

void Model::addPattern(Pattern* pat){
  bounded = false;
  pat->update();
  assert(pat->id != -1);
  if(actPatterns.count(pat)){
    cout << setprecision(16) << DD(pat->val) << endl;
    assert(false);
  }
  assert(pricing->IsCompatible(pat));

  lastVal = pricing->UpdateVal(pat);

  GRBColumn col;
  FOR(i, SZ(pat->its)){
    int f = 1;
    while(i + 1 < SZ(pat->its) && pat->its[i] == pat->its[i + 1]){
      f++;
      i++;
    }
    int it = wToIt[pat->its[i]];
    assert(freq[it] >= f);
    itToIds[it].insert(pat->id);
    auto constr = model->getConstrByName("u" + to_string(it));
    col.addTerm(f, constr);
  }

  cpm->AddSRCoef(col, pat);
  if(insideGT)
    for(const auto& constr: constrsGT)
      if(constr.enable && constr.ids.count(pat->id))
          col.addTerm(1, constr.grbConstr);

  string varname = "x" + to_string(pat->id);
  model->addVar(0.0, GRB_INFINITY, masterCol ? 0.0 : MF, GRB_CONTINUOUS, col, varname);

  pricing->RemPool(pat);
  assert(actPatterns.insert(pat).second);
}

void Model::remPattern(Pattern* pat){
  bounded = false;

  assert(!pat->needUpdate);
  assert(actPatterns.erase(pat));
  auto var = model->getVarByName("x" + to_string(pat->id));
  model->remove(var);

  for(int w : pat->its){
    int it = wToIt[w];
    itToIds[it].erase(pat->id);
  }

  pricing->AddPool(pat);
}

ll Model::FarleyDualBound(){
  ll rcBest = !optPricing ? pricing->MinReducedCost() - 1 : 0;
  assert(rcBest <= 0);

  ll otherBound = partial->Cost() * BASIS;
  FOR(i, NC){
    ll pi = (constrs[i].get(GRB_DoubleAttr_Pi) / MF) * 1.0L * BASIS;
    string name = constrs[i].get(GRB_StringAttr_ConstrName);
    int rhs = constrs[i].get(GRB_DoubleAttr_RHS);
    otherBound += rhs * pi;
  }
  
  if (rcBest != 0)
    otherBound = ((__int128) otherBound  * BASIS) / (BASIS - rcBest) + (((__int128) otherBound * BASIS) % (BASIS - rcBest) == 0 ? 0 : 1);

  return otherBound;
}

void Model::useInitialHeuristic(){
  bestSol = Solution::GetBestInitial(instance);
  storeSolution(bestSol);
  bestInteger = bestSol->Cost();
}

void Model::storeSolution(Solution* sol){
  vi toAdd;
  for(Pattern* bin : sol->GetBins()){
    Pattern* pat = new Pattern(*bin);
    pat->update();
    Pattern::updateId(pat);
    toAdd.push_back(pat->id);
  }
  AddPatternList(toAdd);
}

void Model::rootRuleOptim(){
  if(numDualIneq == N){
    bool good = !masterCol && (safeObj > lastObj + eps10);
    FOR(i, N)
      good &= dualInq[i].get(GRB_DoubleAttr_X) < EPS9;
    if(good){
      lastObj = safeObj;
      ld multiW = safeObj / L0;
      FOR(i, N){
        ld RHS = ceil(multiW * weights[i] * 1e9) / binCapacity / 1e9;
        dualInq[i].set(GRB_DoubleAttr_Obj, MF * RHS);
      }
    }
    if(good){
      if(!PRUNE(GetObj()))
        PH(GetObj());
      ResidueOptim(ROOT);
      RunPL();
    }
  }else if(ResidueOptim(ROOT)){
    lastObj = min(lastObj, safeObj);
    RunPL();
  }
}

bool Model::costOptim(ll curObj){
  assert(bounded);
  if(PRUNE(curObj) || !features.count("Cost"))
    return false;

  vi toRem;
  for(auto pat : actPatterns){
    if(insideGT){
      bool any = false;
      for(const auto& constr : constrsGT)
        if(constr.ids.count(pat->id)){
          any = true;
          break;
        }
      if(any)
        continue;
    }
    if(GetObj() + pricing->UpdateVal(pat) - eps10 > max((bestInteger - 1) * BASIS, GetObj()))
      toRem.push_back(pat->id);
  }

  return RemPatternList(toRem);
}

void Model::computeInvIds(){
  invIds.clear();
  FOR(i, SZ(ids))
    invIds[ids[i]] = i;
}

void Model::updateConstrs(){
  if(constrs != NULL)
    delete[] constrs;
  constrs = model->getConstrs();
  NC = model->get(GRB_IntAttr_NumConstrs);
}

void Model::updateConflictIts(){
  conflictIts.clear();
  FOR(it, SZ(freq))
    if(freq[it] > 0)
      FOR(it2, SZ(freq))
        if(!compatible[it][it2]){
          conflictIts.push_back(it);
          break;
        }
}

void Model::updateIntegers(){
  integers = nonZero = partial->Cost();
  if(bounded)
    FOR(i, M){
      if(ids[i] == -1)
        continue;
      integers += floor(primal[i] + EPS9), nonZero += ceil(primal[i] - EPS9);
    }
}

void Model::checkModelValid(){
  static int cnt = 0;
  cnt++;
  bool firstTry = true;
  try_again:

  FOR(i, M){
    if(ids[i] == -1)
      continue;
    Pattern* pat = idToPat[ids[i]];
    assert(actPatterns.count(pat));
    assert(CanBeAdd(pat));
    if(!insideGT){
      if((primal[i] > EPS9 && !(pricing->UpdateVal(pat) < epsM)) ||
        (primal[i] <= EPS9 && !(pricing->UpdateVal(pat) > -epsM))){

        cout << "NumericError: " << primal[i] << " " << pricing->UpdateVal(pat) << endl;
        if(!firstTry)
          assert(false);
        firstTry = false;
        FastPL();
        goto try_again;
      }
    }
  }
}

ll Model::GetObj(){
  return safeObj;
}

ld Model::UpdateTime(){
  chrono::steady_clock::time_point cur = chrono::steady_clock::now();
  return timeCurrent = chrono::duration_cast<chrono::microseconds>(cur - begin).count() / 1e6L;
}

void Model::Print(bool force, vi st){
  UpdateTime();
  lastPrintBP += !st.empty();
  if(force || timeCurrent - lastPrint >= 5 || lastPrintBP >= 5){
    lastPrint = timeCurrent;
    if(!st.empty()){
      lastPrintBP = 0;
      cout << "Stack State: ";
      FOR(i, SZ(st))
        cout << st[i] << " "[i % 5 != 4];
    }else{
      cout << "ObjPD: " << DD(GetObj());
      cout << " - NumCols: " << model->get(GRB_IntAttr_NumVars);
      cout << " - NumCuts: " << cpm->GetActCuts();
    }
    cout << " - Time: " << setprecision(2) << timeCurrent << setprecision(12);
    cout << " - lastVal: " << DD(lastVal);
    cout << " - Incumbent: " << bestInteger;
    cout << endl;
  }
}

void Model::Update(bool computeInv){
  model->update();
  if(computeInv){
    getPrimal();
    computeInvIds();
  }
}

GRBConstr Model::AddConstr(vi patIds, char sense, int rhs, string name){
  GRBLinExpr linExpr;
  for(int id : patIds)
    if(invIds.count(id))
      linExpr += vars[invIds[id]];
  return model->addConstr(linExpr, sense, rhs, name);
}

void Model::RemConstr(GRBConstr constr){
  model->remove(constr);
}

int Model::GetSzMergeIts(){
  return SZ(mergedIts);
}

int Model::GetNumVars(){
  return M;
}

int Model::GetIntegers(){
  return integers;
}

int Model::GetNumCuts(){
  return cpm->GetActCuts();
}

int Model::GetGenCuts(){
  return cpm->GetCuts();
}

int Model::GetMaxActCuts(){
  return cpm->GetMaxActCuts();
}

int Model::GetRA(){
  return pricing->residueAccept;
}

int Model::GetGenCols(){
  return pricing->generatedCols;
}

ld Model::GetLPTime(){
  return LPTime;
}

ld Model::GetPricingTime(){
  return pricingTime;
}

bool Model::ReachTL(){
  return UpdateTime() >= TIMELIMIT;
}

void Model::ShrinkRA(){
  if (features.count("Waste"))
    pricing->ClearPool(min<int>(binCapacity - 1, floor(max<ld>(EPS10, (bestInteger - L0 - 1 + EPS10)) * binCapacity)));
}

vi Model::GetFreqInst(){
  vi freq2 = instance->freq;
  freq2.resize(SZ(freq));
  return freq2;
}

vi Model::GetFreqCur(){
  return freq;
}

bool Model::AddCuts(ll& curObj){
  assert(bounded);
  if(!features.count("Cut"))
    return false;
  int cnt = 0;
  while(cnt++ < 10){
    CalcAfim();

    if(!cpm->GenCuts())
      break;

    curObj = Optimize(OPTIMIZE);
    if(!optimal && PRUNE(curObj)){
      cout << "Pruned by SR Cuts\n";
      return true;
    }
    if(optimal || (PH(curObj) && PRUNE(curObj))){
      cout << "Improved while cutting" << endl;
      return true;
    }
  }
  return false;
}

ll Model::Optimize(WasteRule rule, ll oldObj){
  cpm->CheckCuts();
  static int cnt = 0;
  bool skip = false;

  if(optimal)
    return RunPL();

  if (oldObj == -1)
    oldObj = L0 * BASIS;

  if(rule == OPTIMIZE){
    ResidueOptim(rule);
    RunPL();
    if(insideRF && !PRUNE(oldObj) && BETTER_EQUAL_THAN(safeObj, oldObj))
      skip = true;
    if(!skip && bounded){
      checkModelValid();
      costOptim(GetObj());
    }
  }

  if(!skip)
    updateConflictIts();

  while(!skip){
    cnt++;
    if(ReachTL())
      throw "TLE";
    RunPL();
    if(insideRF && !PRUNE(oldObj) && BETTER_EQUAL_THAN(safeObj, oldObj)){
      skip = true;
      break;
    }
    if(!masterCol && !PRUNE(GetObj()) && PH(GetObj()) && optimal)
      break;
    if(rule == ROOT)
      rootRuleOptim();
    Print();
    if(colGen())
      break;
  }

  if(insideRF && masterCol){
    auto var = model->getVarByName("M");
    auto vars = model->getVars();
    int NN = model->get(GRB_IntAttr_NumVars);
    FOR(i, NN)
      if(vars[i].get(GRB_StringAttr_VarName)[0] == 'x')
        vars[i].set(GRB_DoubleAttr_Obj, MF);
    delete[] vars;

    model->remove(var);
    masterCol = false;
    Pattern::dfval = 1;
  }

  assert(!masterCol);

  if(bounded && rule == OPTIMIZE)
    checkModelValid();

  computeInvIds();
  if(!optimal && rule == ROOT)
    assert(BETTER_EQUAL_THAN(safeObj, lastObj));

  return GetObj();
}

ll Model::RunPL(){
  try {
    optPricing = false;
    repeat:
    model->update();
    ld bef = UpdateTime();
    model->optimize();
    LPTime += UpdateTime() - bef;
    bounded = model->get(GRB_IntAttr_Status) == GRB_OPTIMAL;

    if(bounded && model->get(GRB_IntAttr_NumVars) > 0 && model->get(GRB_DoubleAttr_DualVio) > EPS9){
      numIssues++;

      model->set(GRB_IntParam_Method, 2);
      model->reset();
      model->update();
      bef = UpdateTime();
      model->optimize();
      LPTime += UpdateTime() - bef;

      model->set(GRB_IntParam_Method, 0);
      model->update();
      bef = UpdateTime();
      model->optimize();
      LPTime += UpdateTime() - bef;
      assert(model->get(GRB_DoubleAttr_DualVio) < EPS9);
    }

    if(masterCol){
      assert(bounded);
      auto var = model->getVarByName("M");
      ld val = var.get(GRB_DoubleAttr_X);
      if(val < EPS9){
        auto vars = model->getVars();
        int NN = model->get(GRB_IntAttr_NumVars);
        FOR(i, NN)
          if(vars[i].get(GRB_StringAttr_VarName)[0] == 'x')
            vars[i].set(GRB_DoubleAttr_Obj, MF);
        delete[] vars;

        model->remove(var);
        masterCol = false;
        Pattern::dfval = 1;

        goto repeat;
      }
      // It can't be worse than using the master column
      assert(model->get(GRB_DoubleAttr_ObjVal) / MF <= 1 + EPS10);
      getDual();
      bounded = false;
      return obj = safeObj = INFLL;
    }

    if(!bounded){
      FOR(it, SZ(freq))
        dual[it] = freq[it] > 0 && itToIds[it].empty() ? 2*BASIS : 0;
      cpm->ResetSR();
      if(accumulate(all(dual), 0) == 0){
        auto vars = model->getVars();
        int NN = model->get(GRB_IntAttr_NumVars);
        FOR(i, NN)
          vars[i].set(GRB_DoubleAttr_Obj, 0.0);
        delete[] vars;

        masterCol = true;
        Pattern::dfval = 0;
        GRBColumn col;
        updateConstrs();

        FOR(i, NC){
          string name = constrs[i].get(GRB_StringAttr_ConstrName);
          int RHS = constrs[i].get(GRB_DoubleAttr_RHS);
          if(name[0] == 'u'){
            int it = atoi(&name[1]);
            assert(freq[it] > 0 && freq[it] == RHS);
            col.addTerm(freq[it], constrs[i]);
          }else if(name[0] == 'H')
            col.addTerm(RHS, constrs[i]);
          else if(name[0] == 's' || name[0] == 'G')
            col.addTerm(RHS, constrs[i]);
          else
            assert(name[0] == 'c' || name[0] == 'L');
        }

        model->addVar(0.0, GRB_INFINITY, MF, GRB_CONTINUOUS, col, "M");
        goto repeat;
      }
      return obj = safeObj = INFLL;
    }

    getPrimal();
    getDual();
    ld objGurobi = partial->Cost() +  model->get(GRB_DoubleAttr_ObjVal) / MF;
    assert(fabs(DD(obj) - objGurobi) < EPS10);
    assert(!optPricing || (DD(safeObj) < objGurobi + EPS10));
    updateIntegers();

    return GetObj();
  } catch(GRBException e) {
    cout << "Error in RunPL\nError code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  exit(1);
}

void Model::FastPL(){
  cpm->CheckCuts();
  model->set(GRB_IntParam_Method, 2);
  model->reset();
  RunPL();
  model->set(GRB_IntParam_Method, 0);
  RunPL();
}

void Model::CalcAfim(){
  assert(bounded);
  afim.assign(SZ(weights), vd(SZ(weights), 0.0));
  FOR(k, M)
    if(primal[k] > EPS9){
      vi& its = idToPat[ids[k]]->its;
      FOR(i, SZ(its)){
        int it1 = wToIt[its[i]];
        for(int j = i + 1; j < SZ(its); j++){
          int it2 = wToIt[its[j]];
          afim[it1][it2] += primal[k];
          if(it1 != it2)
            afim[it2][it1] += primal[k];
        }
      }
    }
}

vi Model::MergeIts(int it1, int it2, int it, vi idsToAdd){
  mergedIts.insert(make_pair(weights[it], ii(weights[it1], weights[it2])));
  AddIts({it}, idsToAdd);
  model->update();
  return RemoveIts({it1, it2});
}

vi Model::SplitIt(int it1, int it2, int it, vi idsToAdd){
  bool find = false;
  for(auto iter = mergedIts.begin(); iter !=  mergedIts.end(); iter++){
    if(iter->first == weights[it] && iter->second == ii(weights[it1], weights[it2])){
      mergedIts.erase(iter);
      find = true;
      break;
    }
  }
  assert(find);

  AddIts({it1, it2}, idsToAdd);
  model->update();
  return RemoveIts({it});
}

void Model::AddConstrGT(Solution* pSol){
  if(insideGT && pSol->Cost() <= constrsGT.front().rhs){
    cout << "Skip add Sol in GT\n";
    return;
  }
  Solution* sol = new Solution(*pSol);
  sol->Unmerge(mergedIts, false);
  toAddGT.push_back(sol);
  cout << "Add Sol in GT: " << sol->Cost() << "\n";
}

bool Model::PH(ll lb){
  int target = ROUND(lb);
  assert(bounded);

  roundingHeuristic->Run(bestSol, partial, mergedIts, target, integers);
  bool find = bestSol->Cost() < bestInteger;
  bestInteger = bestSol->Cost();
  optimal = bestInteger == L2 || bestInteger == L1;
  if(find)
    cout << "Find Integer Solution: " << bestInteger << " - Node LB: " << setprecision(12) << fixed << DD(lb) << endl;
  return target >= bestInteger;
}

void Model::AddIts(const vi& its, vi pats, bool force, bool improveIncumb, bool gt){
  bounded = false;

  vi its2 = its;
  Util::Unique(its2);

  for(int it : its2){
    GRBLinExpr linExp;
    int f = count(all(its), it);
    freq[it] += f;

    if(freq[it] == f){
      model->addConstr(linExp, GRB_GREATER_EQUAL, freq[it], "u" + to_string(it));
      N++;
    }else{
      auto constr = model->getConstrByName("u" + to_string(it));
      constr.set(GRB_DoubleAttr_RHS, freq[it]);
    }
  }
  model->update();
  AddPatternList(pats, force, improveIncumb, gt);
}

vi Model::RemoveIts(vi its){
  bounded = false;

  for(int it : its)
    freq[it]--;

  Util::Unique(its);

  vi idsToRem;
  for(int it : its)
    for(int id : itToIds[it])
      if(!pricing->IsFeasible(idToPat[id]))
        idsToRem.push_back(id);

  RemPatternList(idsToRem);

  for(int it : its){
    assert(freq[it] >= 0);
    auto constr = model->getConstrByName("u" + to_string(it));
    if(freq[it] == 0){
      model->remove(constr);
      N--;
    }else
      constr.set(GRB_DoubleAttr_RHS, freq[it]);
  }

  model->update();
  return idsToRem;
}

void Model::AddPatternList(vi& l, bool force, bool improveIncumb, bool gt){
  Util::Unique(l);

  for(auto id : l){
    Pattern* pat = idToPat[id];
    if(actPatterns.count(pat))
      continue;
    if(CanBeAdd(pat))
      addPattern(pat);
    else{
      assert(!force || gt || improveIncumb);
      pricing->AddPool(pat);
    }
  }
}

bool Model::RemPatternList(vi& l){
  Util::Unique(l);
  for(int id : l)
    remPattern(idToPat[id]);
  return !l.empty();
}

int Model::AppendIt(int wn){
  int it = wToIt[wn] = SZ(weights);
  FOR(i, SZ(weights))
    compatible[i].push_back(true);
  freq.push_back(0);
  weights.push_back(wn);
  compatible.push_back(vc (SZ(weights), true));
  dual.push_back(0);
  itToIds.push_back(si());
  cpm->AddIt();
  assert(freq[it] == 0);
  return it;
}

void Model::AddPoolInactive(Pattern* pat){
  if(!actPatterns.count(pat))
    pricing->AddPool(pat);
}

bool Model::ResidueOptim(WasteRule rule){
  if(optimal || !features.count("Waste"))
    return false;
  if(rule == OPTIMIZE)
    pricing->residueAccept = min<int>(binCapacity - 1, max<int>(0, floor(binCapacity * (bestInteger - L0 - 1 + EPS10))));
  else{
    if(rule == ROOT){
      if(!bounded)
        return false;
      pricing->residueAccept = min<int>(pricing->residueAccept, floor(binCapacity * (max(L1, DD(obj)) - L0 + EPS10)));
    }
  }
  assert(pricing->residueAccept >= 0);

  vi toRem;
  for(auto pat : actPatterns)
    if(pat->residue() > pricing->residueAccept)
      toRem.push_back(pat->id);

  return RemPatternList(toRem);
}

void Model::EnableDualIneq(){
  assert(bounded);

  FOR(i, N)
    assert(constrs[i].get(GRB_StringAttr_ConstrName) == "u" + to_string(i));

  lastObj = bestInteger * BASIS;
  ld multiW = bestInteger / L0;
  FOR(i, N){
    GRBColumn col;
    col.addTerm(1, constrs[i]);
    ld RHS = ceil(multiW * weights[i] * 1e9) / binCapacity / 1e9;
    dualInq.push_back(model->addVar(0.0, GRB_INFINITY, MF * RHS, GRB_CONTINUOUS, col, "d" + to_string(i)));
  }
  numDualIneq = SZ(dualInq);

  int numItemsPricing = 0;
  FOR(it, SZ(freq))
    numItemsPricing += min(binCapacity / weights[it], freq[it]);

  if(numItemsPricing > SZ(freq) * 6 / 5){
    cout << "Initially Binary Pricing - Num Items: " << SZ(freq) << " - Total: " << numItemsPricing << endl;
    pricing->unique = true;
    Optimize(ROOT);
    pricing->unique = false;
    cout << "Finish Binary Pricing" << endl;
  }else{
    Optimize(ROOT);
  }

  numDualIneq = 0;
  for(auto var : dualInq)
    model->remove(var);
  dualInq.clear();

  Print(true);
  Optimize(ROOT);
}

ii Model::GetNumPricingCalls(){
  return ii(pricing->cnt, pricing->cntDP);
}

bool Model::CanBeAdd(Pattern* pat){
  return pricing->CanBeAdd(pat);
}
