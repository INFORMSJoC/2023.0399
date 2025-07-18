#include "CutPlaneManager.hpp"

CPM::CPM(GRBModel*& model, GRBVar*& vars, const int& N, bool& bounded,
  const vi& freq, const vi& ids, const vd& primal,
  const vll& dual, const vector<vd>& afim, const vector<set<int>>& itToIds,
  map<int, int>& invIds):
    model(model), vars(vars), N(N), bounded(bounded), items(Pattern::items), freq(freq),
    ids(ids), primal(primal), dual(dual), afim(afim), itToIds(itToIds),
    invIds(invIds), itsToSRAll(N, vector<IM>()), itsToSR(N, vector<IM>()){}

void CPM::updateItemList(){
  if(SZ(itsToSR) < SZ(items)){
    itsToSR.resize({});
    itsToSRAll.resize({});
  }
}

vector<IM> CPM::getIdsCut(Cut cut){
  map<int, int> acumIdx;
  for(IM im : cut)
    for(int id : itToIds[im.i])
      acumIdx[id] += im.m;

  vector<IM> ret;
  for(auto pr : acumIdx)
    if(pr.second >= cut.d)
      ret.push_back({pr.first, pr.second / cut.d});
  return ret;
}

pair<ld, vector<IM>> CPM::calcCut(Cut cut){
  ld acum = 0;
  vector<IM> idMultCuts = getIdsCut(cut);
  for(IM im : idMultCuts)
    acum += im.m * primal[invIds[im.i]];
  return {acum, idMultCuts};
}

multimap<ld, pair<Cut, vector<IM>>, greater<ld>> CPM::calcCuts(){
  vector<vi> rel(SZ(items));
  FOR(it1, SZ(items))
    if(freq[it1] == 1)
      for(int it2 = it1 + 1; it2 < SZ(items); it2++)
        if(freq[it2] == 1 && afim[it1][it2] > EPS8)
          rel[it1].push_back(it2);

  int cnt = 0;
  multimap<ld, pair<Cut, vector<IM>>, greater<ld>> cuts;
  FOR(it1, SZ(items)){
    for(int it2 : rel[it1]){
      vi ls;
      cnt += SZ(rel[it1]) + SZ(rel[it2]);
      std::set_intersection(all(rel[it1]), all(rel[it2]), back_inserter(ls));
      for(int it3 : ls){
        ld val = afim[it1][it2] + afim[it2][it3] + afim[it3][it1];

        if(val > 1 + EPS8){
          cnt += SZ(itToIds[it1]) + SZ(itToIds[it2]) + SZ(itToIds[it3]);
          vector<IM> ims = {{it1, 1}, {it2, 1}, {it3, 1}};
          sort(all(ims));

          Cut cut(ims, 1, 2);
          pair<ld, vector<IM>> valIneq = calcCut(cut);

          if(valIneq.first > cut.RHS + EPS8)
            cuts.insert({valIneq.first - cut.RHS, {cut, valIneq.second}});
        }
      }
    }
  }
  return cuts;
}

multimap<ld, pair<Cut, vector<IM>>, greater<ld>> CPM::calcCutsUp(){
  int MX = N * binCapacity / 20;

  vector<vi> rel(SZ(items));
  FOR(it1, SZ(items))
    if(freq[it1] == 1)
      FOR(it2, SZ(items))
        if(freq[it2] == 1 && afim[it1][it2] > EPS8)
          rel[it1].push_back(it2);

  int cnt = 0;
  vi itsCut;
  set<vi> cutsVI;
  multimap<ld, pair<Cut, vector<IM>>, greater<ld>> cuts;
  FOR(it1, SZ(items)){
    if(cnt > MX)
      break;

    for(int it2 : rel[it1]){
      if(cnt > MX)
        break;

      itsCut = {it1, it2};
      vi ls1;
      cnt += SZ(rel[it1]) + SZ(rel[it2]);
      std::set_intersection(all(rel[it1]), all(rel[it2]), back_inserter(ls1));

      for(int it3 : ls1){
        if(cnt > MX)
          break;

        if(count(all(itsCut), it3))
          continue;

        itsCut = {it1, it2, it3};
        vi ls2;
        cnt += SZ(rel[it2]) + SZ(rel[it3]);
        std::set_intersection(all(rel[it2]), all(rel[it3]), back_inserter(ls2));

        for(int it4 : ls2){
          if(cnt > MX)
            break;
          if(count(all(itsCut), it4))
            continue;
          itsCut = {it1, it2, it3, it4};

          sort(all(itsCut));
          if(cutsVI.count(itsCut))
            continue;

          rank4(cuts, itsCut, cutsVI, cnt);

          vi ls3;
          cnt += SZ(rel[it3]) + SZ(rel[it4]);
          std::set_intersection(all(rel[it3]), all(rel[it4]), back_inserter(ls3));

          for(int it5 : ls3){
            if(cnt > MX)
              break;

            itsCut = {it1, it2, it3, it4, it5};
            sort(all(itsCut));
            if(cutsVI.count(itsCut))
              continue;
            rank5(cuts, itsCut, cutsVI, cnt);
          }
        }
      }
    }
  }
  return cuts;
}

void CPM::rank4(multimap<ld, pair<Cut, vector<IM>>, greater<ld>>& cuts, vi& itsCut, set<vi>& cutsVI, int& cnt){
  ld affim = 0;
  FOR(i, SZ(itsCut))
    for(int j = i + 1; j < SZ(itsCut); j++)
      affim += afim[itsCut[i]][itsCut[j]];
  if(affim < 1 + EPS8)
    return;

  int szIds = 0;
  for(int it : itsCut)
    szIds += SZ(itToIds[it]);

  vi perm = {1, 1, 1, 2};
  pair<ld, pair<Cut, vector<IM>>> mxCut = {0, {Cut({}, 0, 0), {}}};
  do{
    cnt += szIds;

    vector<IM> ims;
    FOR(i, SZ(itsCut))
      ims.push_back({itsCut[i], perm[i]});
    Cut cut(ims, 1, 3);

    pair<ld, vector<IM>> valIneq = calcCut(cut);
    if(mxCut.first < valIneq.first - cut.RHS - EPS8)
      mxCut = {valIneq.first - cut.RHS, {cut, valIneq.second}};
  }while(next_permutation(all(perm)));

  if(mxCut.first > EPS8){
    cutsVI.insert(itsCut);
    cuts.insert(mxCut);
  }
}

void CPM::rank5(multimap<ld, pair<Cut, vector<IM>>, greater<ld>>& cuts, vi& itsCut, set<vi>& cutsVI, int& cnt){
  int szIds = 0;
  for(int it : itsCut)
    szIds += SZ(itToIds[it]);

  ld affim = 0;
  FOR(i, SZ(itsCut))
    for(int j = i + 1; j < SZ(itsCut); j++)
      affim += afim[itsCut[i]][itsCut[j]];
  if(affim < 2 + EPS8)
    return;

  vvi perms = {{1, 1, 1, 2, 2}, {1, 1, 1, 1, 3}, {1, 1, 2, 2, 3}, {1, 1, 2, 2, 2}, {1, 2, 2, 3, 3}};
  vi dem = {4, 4, 5, 3, 4};
  vi RHSs = {1, 1, 1, 2, 2};
  //vvi perms = {{1, 1, 1, 1, 1}, {2, 2, 2, 2, 2}, {1, 1, 2, 2, 2}};
  //vi dem = {2, 3, 3};
  //vi RHSs = {2, 3, 2};
  pair<ld, pair<Cut, vector<IM>>> mxCut = {0, {Cut({}, 0, 0), {}}};
  FOR(i, SZ(perms)){
    vi perm = perms[i];
    int d = dem[i];
    int RHS = RHSs[i];
    do{
      cnt += szIds;

      vector<IM> ims;
      FOR(i, SZ(itsCut))
        ims.push_back({itsCut[i], perm[i]});

      Cut cut(ims, RHS, d);
      pair<ld, vector<IM>> valIneq = calcCut(cut);

      if(mxCut.first < valIneq.first - RHS - EPS8)
        mxCut = {valIneq.first - RHS, {cut, valIneq.second}};
    }while(next_permutation(all(perm)));
  }
  if(mxCut.first > EPS8){
    cutsVI.insert(itsCut);
    cuts.insert(mxCut);
  }
}

int CPM::addCutPool(Cut cut){
  if(cutsPool.count(cut)){
    cut = *cutsPool.find(cut);
    countSR[cut.id] = 0;
    assert(mapSR[cut.id] == -2);
    mapSR[cut.id] = -1;
  }else{
    int id = cut.id = SZ(mapSR);
    mapSR.push_back(-1);
    countSR.push_back(0);
    setsSRAll.push_back(cut);
    for(IM im : setsSRAll.back())
      itsToSRAll[im.i].push_back({id, im.m});
    cutsPool.insert(cut);
  }
  return cut.id;
}

void CPM::AddSRCoef(GRBColumn& col, Pattern* pat){
  vi cnt(SZ(setsSRAll), 0);

  for(int it : pat->its)
    for(IM im : itsToSRAll[it])
      cnt[im.i] += im.m;

  FOR(id, SZ(cnt))
    if(cnt[id] >= setsSRAll[id].d && mapSR[id] >= -1){
      auto constr = model->getConstrByName("c" + to_string(id));
      col.addTerm(cnt[id] / setsSRAll[id].d, constr);
    }
}

void CPM::AddIt(){
  itsToSR.push_back({});
  itsToSRAll.push_back({});
}

bool CPM::GenCuts(){
  for(auto pr : invIds){
    int id = stoi(&vars[pr.second].get(GRB_StringAttr_VarName)[1]);
    assert(pr.first == id);
  }

  updateItemList();
  auto cuts = calcCuts();

  //if(cuts.empty())
  //  cuts = calcCutsUp();

  int allowed = 20; //max(1, min(20, MAX_CUT_ALLOWED - GetActCuts()));
  if(SZ(cuts) > allowed)
    cuts.erase(next(cuts.begin(), allowed), cuts.end());

  if(cuts.empty())
    return false;

  for(auto& pr : cuts){
    const Cut& cut = pr.second.first;
    int id = addCutPool(cut);

    GRBLinExpr constr;
    for(IM im : pr.second.second)
      constr += im.m * vars[invIds[im.i]];
    model->addConstr(constr, GRB_LESS_EQUAL, cut.RHS, "c" + to_string(id));
  }

  return true;
}

int CPM::SetDualSR(ll val, int it){
  if(-val > 0){
    int id = SZ(dualSR);
    mapSR[it] = id;
    countSR[it] = 0;
    dualSR.push_back(val);
    seenSR.push_back(0);
    setsSR.push_back(setsSRAll[it]);
    for(IM im : setsSR.back())
      itsToSR[im.i].push_back({id, im.m});
    countSR[it] += -val < eps13;
  }else{
    countSR[it]++;
    mapSR[it] = -1;
  }
  actCuts++;
  return setsSRAll[it].RHS;
}

void CPM::ResetSR(){
  actCuts = 0;
  dualSR.clear();
  setsSR.clear();
  seenSR.clear();
  FOR(it, SZ(freq))
    itsToSR[it].clear();
}

void CPM::CheckCuts(){
  FOR(id, SZ(setsSRAll)){
    if(mapSR[id] >= -1){
      bool good = true;
      for(IM im : setsSRAll[id])
        good &= freq[im.i] == 1;
      if(!good || countSR[id] > 100){
        bounded = false;
        mapSR[id] = -2;
        auto constr = model->getConstrByName("c" + to_string(id));
        model->remove(constr);
      }
    }
  }
}

int CPM::GetActCuts(){
  return actCuts;
}

int CPM::GetCuts(){
  return SZ(mapSR);
}

int CPM::GetMaxActCuts(){
  return mxActCuts;
}

void CPM::UpdMaxActCuts(){
  mxActCuts = max(mxActCuts, actCuts);
}
