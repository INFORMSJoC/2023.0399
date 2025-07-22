#include "BranchManager.hpp"

void BranchManager::BranchToLeft(ll& curObj){
  assert(lvl == SZ(st));
  pair<ii, int> ans = getItsToMerge();
  ii pr = ans.first;
  mergeIts(pr.first, pr.second, ans.second);
  curObj = model->Optimize(OPTIMIZE, curObj);
  numBranches++;
}

pair<ii, int> BranchManager::getItsToMerge(){
  model->CalcAfim();
  
  map<ii, int> prior;
  if (features.count("Historic"))
    FOR(i, SZ(globalPrior))
      prior[globalPrior[i]] = SZ(globalPrior) - i;

  vector<tuple<int, int, ii>> vans;
  FOR(it1, SZ(weights))
    FOR(it2, SZ(weights)){
      ld val = afim[it1][it2], valF = FAC(val);
      if(valF < 1e-3)
        continue;
      if(freq[it1] > 0 && freq[it2] > (it1 == it2) && weights[it1] >= weights[it2]){
        int dp = prior.count(ii(it1, it2)) ? prior[ii(it1, it2)] : 0;
        int w = weights[it1] + weights[it2];
        vans.push_back({dp, w, ii(it1, it2)});
      }
    }

  if(vans.empty()){
    FOR(i, SZ(model->primal)){
      ld val = model->primal[i], valF = FAC(val);
      if(valF < 1e-3)
        continue;
      Pattern* pat = idToPat[model->ids[i]];
      FOR(j, SZ(pat->its))
        REV(k, j)
          vans.push_back({max(pat->its[j], pat->its[k]), i, ii(wToIt[pat->its[j]], wToIt[pat->its[k]])});
    }
    assert(!vans.empty());
  }
  ii ret = get<ii>(*max_element(all(vans)));
  int qt = min({freq[ret.first], freq[ret.second], (int)ceil(afim[ret.first][ret.second] - 1e-7)});
  assert(qt > 0);
  return pair<ii, int>(ret, ret.first == ret.second ? 1 : qt);
}

void BranchManager::mergeIts(int it1, int it2, int qt){
  FOR(q, qt){
    assert(compatible[it1][it2]);
    st.push_back(Node(it1, it2));
    Node& node = st.back();
    lvl++;

    int w1 = weights[it1], w2 = weights[it2];

    assert(freq[it1] >= 1 + (it1 == it2));
    assert(freq[it2] >= 1);

    try{
      int it, wn = w1 + w2;
      if(wToIt[wn] != -1){
        it = wToIt[wn];
      }else{
        it = model->AppendIt(wn);
      }

      node.it = it;

      idsMergedPats(it1, it2, node.join);
      vi remIds = model->MergeIts(it1, it2, it, node.join);
      findSplitAndBT(remIds, node.split, node.backtrack, w1, w2);

      auto diffConfl = conflictDAG.AddMerge(it1, it2, it);
      Util::Append(node.split, conflictsPropagation(diffConfl));
      if(q + 1 != qt)
        model->Update();
    } catch(GRBException e) {
      cout << "Error in merge itens" << endl;
      cout << "Error code = " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch(...){
      cout << "Error in merge itens" << endl;
    };
  }
}

Pattern* BranchManager::mergeItPat(Pattern* pat, int w1, int w2){
  Pattern* nw = new Pattern(w1 + w2);
  bool fw1 = false, fw2 = false;
  for(int w : pat->its){
    if(!fw1 && w == w1){
      fw1 = true;
      continue;
    }
    if(!fw2 && w == w2){
      fw2 = true;
      continue;
    }
    nw->add(w);
  }
  return nw;
}

void BranchManager::idsMergedPats(int it1, int it2, vi& idsPats){
  vi interc;
  set_intersection(all(itToIds[it1]), all(itToIds[it2]), back_inserter(interc));
  int w1 = weights[it1], w2 = weights[it2];
  for(int id : interc){
    if(it1 == it2 && count(all(idToPat[id]->its), w1) < 2)
      continue;
    Pattern* pat = new Pattern(*idToPat[id]);
    pat->rem(w1);
    pat->rem(w2);
    pat->add(w1 + w2);
    Pattern::updateId(pat);
    idsPats.push_back(pat->id);
  }
}

void BranchManager::idsSplitPats(int it1, int it2, vi& idsPats){
  int w1 = weights[it1], w2 = weights[it2], it = wToIt[w1 + w2];
  for(int id : itToIds[it]){
    Pattern* pat = new Pattern(*idToPat[id]);
    pat->rem(w1 + w2);
    pat->add(w1);
    pat->add(w2);
    Pattern::updateId(pat);
    idsPats.push_back(pat->id);
  }
}

void BranchManager::findSplitAndBT(const vi& from, vi& split, vi& bt, int w1, int w2){
  for(int id : from)
    if(count(all(idToPat[id]->its), w1) == 0 || count(all(idToPat[id]->its), w2) < (1 + (w1 == w2)))
      split.push_back(id);
    else
      bt.push_back(id);
}

vi BranchManager::reCalcAndPropag(){
  return conflictsPropagation(conflictDAG.ReCalc());
}

vi BranchManager::conflictsPropagation(const vector<ii>& diff){
  vi idsToRem;
  for(ii pr : diff)
    getConflictList(idsToRem, pr.first, pr.second);

  model->RemPatternList(idsToRem);
  return idsToRem;
}

void BranchManager::getConflictList(vi& co, int it1, int it2){
  if(it1 != it2){
    set_intersection(all(itToIds[it1]), all(itToIds[it2]), back_inserter(co));
  }else{
    int w = weights[it1];
    for(int id : itToIds[it1])
      if(count(all(idToPat[id]->its), w) >= 2)
        co.push_back(id);
  }
}

void BranchManager::swapState(){
  bounded = false;

  Node& node = st.back();
  node.state ^= true;
  int it1 = node.it1, it2 = node.it2, it = node.it;

  compatible[it1][it2] = compatible[it2][it1] = !node.state;
  
  if(node.state){
    conflictDAG.AddConflict(it1, it2);
    node.join = model->SplitIt(it1, it2, it, node.split);
  }
  else{
    node.split = model->MergeIts(it1, it2, it, node.join);
  }
  Util::Append(node.state? node.join : node.split, reCalcAndPropag());
}

int BranchManager::backtracking(){
  static int cnt = 0;
  cnt++;
  swapState();
  model->Update();

  vi rem(SZ(st), 0);
  if (features.count("Splay"))
    REV(i, SZ(st) - 1){
      if(st[i].fixed || st[i].state)
        break;
      rem[i] = true;
    }

  vector<vi> m(SZ(st), freq);

  REV(i, SZ(st) - 1){
    m[i] = m[i + 1];
    if(!st[i].state){
      m[i][st[i].it1]++;
      m[i][st[i].it2]++;
      m[i][st[i].it]--;
    }
  }

  REV(i, SZ(st)){
    if(rem[i]){
      assert(!st[i].state);
      for(int j = i + 1; j < SZ(st); j++){
        assert(!st[j].state);
        if(!rem[j]){
          m[j][st[i].it1]++;
          m[j][st[i].it2]++;
          m[j][st[i].it]--;
          if(m[j][st[i].it] < 0){
            rem[i] = false;
            for(; j >= i + 1; j--)
              if(!rem[j]){
                m[j][st[i].it1]--;
                m[j][st[i].it2]--;
                m[j][st[i].it]++;
              }
            break;
          }
        }
      }
    }
  }
  
  int eliminated = accumulate(all(rem), 0);
  if(eliminated == 1 && rem[SZ(rem) - 2])
    eliminated = 0;

  {
    int it1 = st.back().it1, it2 = st.back().it2;
    auto iter = find(all(globalPrior), ii(it1, it2));
    if(iter == globalPrior.end())
      globalPrior.push_back(ii(it1, it2));
    else if(iter != globalPrior.begin()){
      auto iter2 = prev(iter);
      swap(*iter, *iter2);
    }
  }
  if(eliminated){
    st.back().fixed = true;

    REV(i, SZ(rem))
      if(rem[i])
        remPosStack(i);
  }else{
    int szSt = SZ(st);
    remPosStack(SZ(st) - 1);
    while(!st.empty() && st.back().state)
      remPosStack(SZ(st) - 1);
    if(st.empty())
      return (optimal = true, szSt);
    swapState();
    eliminated = szSt - SZ(st);
  }
  Util::Append(st.back().backtrack, reCalcAndPropag());

  vi freq2 = model->GetFreqInst();
  for(auto s : st)
    if(!s.state){
      freq2[s.it1]--;
      freq2[s.it2]--;
      assert(freq2[s.it1] >= 0);
      assert(freq2[s.it2] >= 0);
      freq2[s.it]++;
    }
  assert(freq2 == freq);

  return eliminated;
}

void BranchManager::remPosStack(int p){
  bounded = false;

  Node& node = st[p];
  int it1 = node.it1, it2 = node.it2, it = node.it;

  if(node.state)
    compatible[node.it1][node.it2] = compatible[node.it2][node.it1] = true;

  if(!node.state){
    idsSplitPats(it1, it2, node.backtrack);
    Util::Append(node.backtrack, node.split);
    model->SplitIt(it1, it2, it, node.backtrack);
  }else
    model->AddPatternList(node.backtrack);

  st.erase(st.begin() + p);
  lvl--;
}

bool BranchManager::FixST(ll& curObj, bool print){
  assert(!st.empty());
  static int cnt = 0;
  cnt++;

  while(true){
    assert(lvl == SZ(st));
    if(!PRUNE(curObj))
      break;

    assert(model->optPricing);
    print = true;
    if(!st.back().state){
      numBranches++;
      cout << "Swaping: " << DD(curObj) << " - Its " << st.back().it1 << " and " << st.back().it2 << " - ST: " << st.back().state << " - lvl " << lvl << endl;
      swapState();
    }else{
      cout << "Fail: " << DD(curObj) << " - Its " << st.back().it1 << " and " << st.back().it2 << " - ST: " << st.back().state;
      cout << " - Subiu " << backtracking() << " niveis" << " - cur lvl " << lvl << endl;
      if(optimal)
        return true;
    }
    curObj = model->Optimize(OPTIMIZE);
    if(optimal)
      return true;
  }
  if(!st.back().state){
    int it1 = st.back().it1, it2 = st.back().it2;
    auto iter = find(all(globalPrior), ii(it1, it2));
    if(iter != globalPrior.end()){
      auto iter2 = next(iter);
      if(iter2 != globalPrior.end())
        swap(*iter, *iter2);
    }
  }
  //if(!st.back().state)
  //  depthRem.erase(ii(st.back().it1, st.back().it2));

  if(print){
    cout << "Sucess: " << DD(curObj) << " - Its " << st.back().it1 << " and " << st.back().it2;
    cout << " - ST: " << st.back().state << " - Integers: " << model->GetIntegers();
    cout << " - NumCols: " << model->GetNumVars() << " - NumCuts: " << model->GetNumCuts();
    cout << " - lvl " << lvl << endl;
    Print();
  }
  return false;
}

ll BranchManager::CleanStack(bool hardReset, bool run){
  oldST.clear();
  if(st.empty())
    return run? model->Optimize(OPTIMIZE) : model->GetObj();

  int ini = 0;
  while(ini < SZ(st) && st[ini].state)
    ini++;

  if(!hardReset){
    if(ini < SZ(st) && !st[ini].state){
      ini++;
      bool any = false;
      while(ini < SZ(st) && st[ini].state)
        ini++, any = true;
      if(!any)
        ini--;
    }
    //lvlIts[st[ini].it1] = 0;
  }else if(run)
    oldST = st;

  for(int i = SZ(st) - 1; i >= ini; i--)
    remPosStack(i);
  reCalcAndPropag();

  if(hardReset){
    vi freq2 = model->GetFreqInst();
    assert(freq2 == freq);
  }

  if(run){
    model->FastPL();
    return model->Optimize(OPTIMIZE);
  }

  return model->GetObj();
}

void BranchManager::RecoverST(){
  for(int i = SZ(st); i < SZ(oldST); i++){
    const auto& s = oldST[i];
    mergeIts(s.it1, s.it2, 1);
    if(s.state)
      swapState();
    st[i].fixed = oldST[i].fixed;
  }
  oldST.clear();
}

void BranchManager::Print(bool force){
  vi bST;
  for(auto node: st)
    bST.push_back(node.state);
  model->Print(force, bST);
}

bool BranchManager::AllRight(){
  for(auto& node : st)
    if(!node.state)
      return false;
  return true;
}

int BranchManager::GetLvl(){
  return lvl;
}

int BranchManager::GetNumBranches(){
  return numBranches;
}
