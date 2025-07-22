#include "RFHeuristic.hpp"

void RFHeuristic::fixVars(vector<NodeRF>& st, vi idsToRem){
  vi freq2 = model->GetFreqCur();
  vi its;
  st.push_back(NodeRF());
  int remToConstrGT = 0;

  assert(!insideGT || model->constrsGT[0].enable);

  FOR(h, SZ(idsToRem)){
    int id = idsToRem[h];
    bool any = false, alls = true;
    const vi& patIts = idToPat[id]->its;

    FOR(i, SZ(patIts)){
      int f = 1;
      while(i + 1 < SZ(patIts) && patIts[i] == patIts[i + 1]){
        f++;
        i++;
      }
      int it = wToIt[patIts[i]];
      any |= freq2[it] > 0;
      alls &= freq2[it] >= f;
    }

    assert(any);
    if(!its.empty() && !alls)
      continue;

    Pattern* pat = new Pattern();
    for(int w : patIts){
      int it = wToIt[w];
      if(freq2[it] > 0){
        freq2[it]--;
        its.push_back(it);
        pat->add(w);
      }
    }
    partial->InsertBin(pat);
    
    if(insideGT && model->constrsGT[0].ids.count(id)){
      auto& constr = model->constrsGT[0];
      int RHS = constr.grbConstr.get(GRB_DoubleAttr_RHS);
      constr.grbConstr.set(GRB_DoubleAttr_RHS, RHS - 1);
      model->Update();
      remToConstrGT++;
    }

    st.back().partialBins.push_back(pat);
    st.back().sz++;
  }
  assert(!its.empty());

  st.back().its = its;
  st.back().patIds = model->RemoveIts(its);
}

void RFHeuristic::fixPatsVar(vector<NodeRF>& st){
  static int cnt = 0;
  cnt++;

  assert(!insideGT || model->constrsGT[0].enable);

  vi idsToRem;
  vector<tuple<ld, int, int>> idsFrat;

  FOR(i, M){
    ld val = primal[i];
    int id = ids[i];
    assert(id != -1);
    while(val > 1 - EPS9)
      idsToRem.push_back(id), val--;
    val = round(100 * val) / 100;
    if(val > EPS9)
      idsFrat.push_back({val, idToPat[id]->its.front(), id});
  }

  if(idsToRem.empty()){
    ld gap = (DD(model->GetObj()) - bestInteger) - 1 - EPS10;
    if(gap <= 0)
      idsToRem = {get<2>(*max_element(all(idsFrat)))};
    else{
      sort(rall(idsFrat));
      vi freq2 = model->GetFreqCur();
      for(auto pr : idsFrat){
        ld val = get<0>(pr);
        if((1 - val > gap || val < 0.51) && !idsToRem.empty())
          break;

        int id = get<2>(pr);
        bool alls = true;
        for(int w : idToPat[id]->its)
          alls &= freq2[wToIt[w]] > 0, freq2[wToIt[w]]--;

        if(!alls && !idsToRem.empty()){
          for(int w : idToPat[id]->its)
            freq2[wToIt[w]]++;
          continue;
        }

        gap -= 1 - val;
        idsToRem.push_back(id);
      }
    }
  }
  assert(!idsToRem.empty());
  fixVars(st, idsToRem);
}

void RFHeuristic::remBackNodeRF(vector<NodeRF>& stRF, int h, bool improveIncumb, bool gt){
  NodeRF& node = stRF[h];

  model->AddIts(node.its, node.patIds, true, improveIncumb, gt);
  for(auto bin : node.partialBins){
    partial->RemoveBin(bin);
    if(insideGT){
      Pattern* pat = Pattern::getTracked(bin);
      auto& constr = model->constrsGT[0];
      assert(constr.enable);
      if(constr.ids.count(pat->id)){
        int RHS = constr.grbConstr.get(GRB_DoubleAttr_RHS);
        constr.grbConstr.set(GRB_DoubleAttr_RHS, RHS + 1);
        model->Update();
      }
    }
    delete bin;
  }
  stRF.erase(stRF.begin() + h);
}

bool RFHeuristic::Run(ll& curObj){
  static int cnt = 0;
  insideRF = true;

  cout << "Running Relax and Fix" << endl;
  int iterations = 0, maxIntegers = 0;
  int runsGT = 2 * insideGT;
  int prevBestInteger = bestInteger;
  vector<NodeRF> stRF;

  resetGT:
  int runsRF = 3;
  curObj = model->Optimize(OPTIMIZE);
  maxIntegers = max(maxIntegers, integers);

  resetRF:
  runsRF--;
  
  //if(insideGT)
  //  for(int id : model->constrsGT[0].ids)
  //    assert(model->actPatterns.count(idToPat[id]));

  while(!optimal){
    assert(bounded);
    iterations++;
    cnt++;

    fixPatsVar(stRF);

    curObj = model->Optimize(OPTIMIZE, curObj);
    
    if(!bounded || PRUNE(curObj))
      break;
      
    model->PH(curObj);
    maxIntegers = max(maxIntegers, integers);
  }

  model->ResidueOptim(OPTIMIZE);
  bool last = true, improveIncum = bestInteger > prevBestInteger;

  int stSZ = SZ(stRF);
  REV(h, SZ(stRF)){
    if(!optimal && last && stSZ == 1)
      model->AddConstrGT(partial), last = false;
    if(improveIncum || runsRF == 0 || h >= 3 * stSZ / 4)
      remBackNodeRF(stRF, h, improveIncum, true);
    if(!optimal && last)
      model->AddConstrGT(partial), last = false;
  }
  if(!optimal && runsRF > 0){
    curObj = model->Optimize(OPTIMIZE);
    goto resetRF;
  }

  assert(stRF.empty() && partial->Cost() == 0);

  if(!optimal){
    if(improveIncum || --runsGT > 0){
      prevBestInteger = bestInteger;
      if(insideGT){
        auto& constr = model->constrsGT[0];
        constr.grbConstr.set(GRB_DoubleAttr_RHS, constr.rhs - 12);
      }
      goto resetGT;
    }
    insideRF = false;
    model->FastPL();

    curObj = model->Optimize(OPTIMIZE);
    model->PH(curObj);
  }else
    insideRF = false;

  if(optimal)
    cout << "Find Optimal Solution in Relax and Fix" << endl;
  else{
    assert(model->bounded);
    cout << "Max integers variables: " << maxIntegers << endl;
  }
  
  assert(partial->Cost() == 0);
  return optimal;
}
