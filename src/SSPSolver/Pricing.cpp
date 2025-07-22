#include "Pricing.hpp"

#define BEST_THAN(val1, val2) ((val1) > (val2))
#define WORST_EQUAL_THAN(val1, val2) ((val1) <= (val2))
#define IMPROVE(val) BEST_THAN(val, epsM - 1)

Pricing::Pricing(CPM* cpm, vi& weights, vi& freq, vi& wToIt, vi& conflictIts, vll& v, vector<vc>& compatible, bool& insideGT):
    W(binCapacity), weights(weights), freq(freq), wToIt(wToIt), conflictIts(conflictIts),
    v(v), vSR(cpm->dualSR), itsToSR(cpm->itsToSR), setsSR(cpm->setsSR), compatible(compatible),
    seenSR(cpm->seenSR), insideGT(insideGT), capAccepted(2 * binCapacity - 1){}

Pricing::~Pricing(){
  for(Pattern* pat : pool)
    delete pat;
}

vector<Pattern*> Pricing::getFromPool(){
  vector<Pattern*> bests;

  for(auto pat : pool){
    if(pat->getCap() > capAccepted)
      continue;
    bool good = IsFeasible(pat);
    if(good && IMPROVE(pat->val))
      if(isValid(pat) && IMPROVE(pat->val))
        bests.push_back(pat);
  }
  dryBests(bests, false);
  return bests;
};

void Pricing::dryBests(vector<Pattern*>& bests, bool del){
  if(!bests.empty()){
    struct compVal{
      bool operator () (Pattern* p1, Pattern* p2) const{
        return BEST_THAN(p1->val, p2->val);
      }
    };

    vector<Pattern*> ret;
    sort(all(bests), compVal());
    freqRet.assign(SZ(freq), 0);
    assert(!del);

    for(Pattern* pat : bests){
      assert(pat->id != -1);
      assert(isValid(pat, pat->val));
      bool can = true;
      for(int w : pat->its)
        can &= freqRet[wToIt[w]] < mxFreqRet;
      if(can){
        for(int w : pat->its)
          freqRet[wToIt[w]]++;
        ret.push_back(pat);
      }else{
        pool.insert(pat);
      } 
    }
    ret.swap(bests);
  }
}

bool Pricing::isValid(Pattern* pat, ll val2){
  bool isGood = true;
  FOR(i, pat->sz())
    for(int j = i + 1; j < pat->sz(); j++)
      isGood &= compatible[wToIt[pat->its[i]]][wToIt[pat->its[j]]];

  if(!isGood)
    return false;

  int cap = 0;
  UpdateVal(pat);
  for(int w : pat->its)
    cap += weights[wToIt[w]];

  assert(val2 == INFLL || pat->val == val2);
  assert(cap == pat->getCap() && cap >= binCapacity);

  return true;
}

void Pricing::recoverBest(int i, Pattern* pat){
  if(i <= 0 || pat->getCap() >= binCapacity){
    assert(pat->getCap() <= capAccepted);
    assert(pat->getCap() >= binCapacity && IMPROVE(pat->val));
    Pattern* temp = new Pattern(*pat);
    Pattern::updateId(temp);

    if(insideGT && gtSet.count(temp->id))
      return;
    if(incumb->id == -1)
      delete incumb;

    incumb = temp;
    return;
  }

  trysDP--;
  if(trysDP < 0)
    return;

  int it = its[i], iu = i - 1, w2 = pat->getCap() + weights[it];
  while(iu > 0 && weights[its[iu]] == weights[it])
    iu--;

  bool swapOrder = false;
  const ll dpRight = pat->val - dp[iu][pat->getCap()];
  if(w2 <= capAccepted){
    const ll dpLeft = pat->val - dp[i - 1][w2] - v[it];
    if(BEST_THAN(dpLeft, incumb->val)){

      if(BEST_THAN(dpRight, dpLeft)){
        recoverBest(iu, pat);
        swapOrder = true;
        if(trysDP < 0)
          return;
      }

      if(canPack(pat, it)){
        addIt(pat, it);
        if(BEST_THAN(pat->val - dp[i - 1][w2], incumb->val)){
          recoverBest(i - 1, pat);
          if(trysDP < 0)
            return;
        }
        remIt(pat, it);
      }
    }
  }

  if(!swapOrder && BEST_THAN(dpRight, incumb->val))
    recoverBest(iu, pat);
}

void Pricing::recoverDPLexic(int i, Pattern* pat){
  if(i <= 0 || pat->getCap() >= capAccepted){
    assert(pat->getCap() <= capAccepted);
    Pattern* bin = new Pattern(*pat);
    Pattern::updateId(bin);
    if(insideGT && gtSet.count(bin->id))
      return void(goodAll = false);

    retDP.push_back(bin);
    for(int w : pat->its)
      freqRet[wToIt[w]]++;
    return;
  }

  trysDP--;
  if(trysDP < 0 && !retDP.empty())
    return;

  int it = its[i], w2 = pat->getCap() + weights[it];

  if(w2 <= capAccepted){
    ll dpTo = pat->val - dp[i - 1][w2] - v[it];
    if(IMPROVE(dpTo)){
      if(canPack(pat, it)){
        addIt(pat, it);
        if(IMPROVE(pat->val - dp[i - 1][w2]))
          recoverDPLexic(i - 1, pat);
        else
          goodAll = false;
        remIt(pat, it);
      }else
        goodAll = false;
    }
  }

  for(int w : pat->its)
    if(freqRet[wToIt[w]] > 2 * mxFreqRet)
      return;

  while(i > 0 && weights[its[i]] == weights[it])
    i--;

  if(IMPROVE(pat->val - dp[i][pat->getCap()]))
    recoverDPLexic(i, pat);
}

bool Pricing::updateDPsize(int N, ll INFMAX){
  if (N > SZ(dp) || SZ(dp[0]) < capAccepted + 1){
    ll dpSize = 8LL * N  * (capAccepted + 1);
    if(dpSize >= 60'000'000'000LL){
      cerr << "N: " << N << " - Cap: " << capAccepted << " - DP size: " << dpSize << endl;
      throw runtime_error("Memory limit reached in DP");
    }
    dp.clear();
    dp.assign(N, vll(capAccepted + 1, INFMAX));
    return true;
  }
  return false;
}

vector<Pattern*> Pricing::DP(bool& opt, bool& feasible){
  cntDP++;
  int N = its.size(); 
  const ll INFMAX = LONG_LONG_MAX / 2;//2 * accumulate(all(freq), 0) * BASIS;

  long double acumAll = 0;
  FOR(i, SZ(its))
    acumAll += v[its[i]];
  assert(2 * acumAll - 1 < INFMAX);

  its.insert(its.begin(), SZ(freq));

  if(!updateDPsize(SZ(its), INFMAX))
    FOR(i, N + 1)
      fill(all(dp[i]), INFMAX);

  for(int w = binCapacity; w <= capAccepted; w++)
    dp[0][w] = 0;

  FORI(i, N){
    int it = its[i], wg = weights[it];
    FOR(w, capAccepted + 1)
      dp[i][w] = min(dp[i - 1][w], w + wg <= capAccepted ? dp[i - 1][w + wg] + v[it] : INFMAX);
  }

  opt = !IMPROVE(Pattern::GetDefaultValue() - dp[N][0]);
  if(opt){
    feasible = Pattern::GetDefaultValue() - dp[N][0] <= 0;
    return {};
  }

  
  Pattern* temp = new Pattern();
  trysDP = N * 1LL * W / 20;
  goodAll = true;
  retDP.clear();
  fill(all(seenSR), 0);
  freqRet.assign(SZ(freq), 0);

  if(!features.count("Best")){
    recoverDPLexic(N, temp);

    if(!retDP.empty()){
      FOR(i, SZ(retDP))
        if(BEST_THAN(retDP[i]->val, retDP[0]->val))
          swap(retDP[i], retDP[0]);
      incumb = new Pattern(*retDP[0]);
      assert(IMPROVE(incumb->val));

      trysDP = N * 1LL * W / 50;
      recoverBest(N, temp);
      if(BEST_THAN(incumb->val, retDP[0]->val))
        retDP.push_back(incumb);
      else 
        delete incumb;
    }
  }else{
    incumb = new Pattern(*temp);
    incumb->val = epsM - 1;
    recoverBest(N, temp);
    goodAll = false;
    if(IMPROVE(incumb->val))
      retDP.push_back(incumb);
    else
      delete incumb;
  }
  delete temp;

  if(retDP.empty()){
    if(goodAll){
      cout << setprecision(16) << dp[N][0] <<  " - N: " << N << " - cnt: "  << cntDP << endl;
      int i = N, w = 0;
      ll val = Pattern::GetDefaultValue();
      while(i > 0 && w < binCapacity){
        int it = its[i];
        if(weights[it] <= w && dp[i][w] == (dp[i - 1][w + weights[it]] - v[it])){
          val -= v[it];
          w += weights[it];
          cout << it << " - v: " << v[it] << " - w: " << weights[it] << " - V: "  <<  val << " - W: " << w  << " - i: " << i << endl;
        }
        i--;
      }
      throw runtime_error("Nao achou a coluna boa");
    }
    opt = true;
  }
  for(auto pat : retDP)
    assert(IMPROVE(pat->val));
  dryBests(retDP, false);
  opt |= retDP.empty();
  return retDP;
}

bool Pricing::canPack(Pattern* bin, int it){
  for(int w2 : bin->its)
    if(!compatible[it][wToIt[w2]])
      return false;
  return true;
}

void Pricing::addItVal(Pattern* pat, int it){
  pat->val -= v[it];

  for(IM im: itsToSR[it]){
    int id = im.i;
    seenSR[id] += im.m;
    if(seenSR[id] >= setsSR[id].d){
      pat->val -= vSR[id];
      seenSR[id] -= setsSR[id].d;
    }
  }
}

void Pricing::addIt(Pattern* pat, int it){
  pat->add(weights[it]);
  addItVal(pat, it);
}

void Pricing::remIt(Pattern* pat, int it){
  pat->rem(weights[it]);
  pat->val += v[it];

  for(IM im: itsToSR[it]){
    int id = im.i;
    seenSR[id] -= im.m;
    if(seenSR[id] < 0){
      pat->val += vSR[id];
      seenSR[id] += setsSR[id].d;
    }
  }
}

void Pricing::AddPool(Pattern* pat, bool force){
  auto pr = pool.insert(pat);
  assert((*pr.first)->id == pat->id);
  assert(!force || pr.second);
}

void Pricing::RemPool(Pattern* pat, bool force){
  assert(pool.erase(pat) || !force);
}

void Pricing::ClearPool(int cap){
  for(auto iter = pool.begin(); iter != pool.end();){
    Pattern* pat = *iter;
    if(pat->getCap() > cap){
      iter = pool.erase(iter);
      Pattern::allPats.erase(pat);
      Pattern::idToPat[pat->id] = NULL;
      delete pat;
    }else
      iter++;
  }
}

bool Pricing::IsCompatible(Pattern* pat){
  FOR(i, pat->sz()){
    int it1 = wToIt[pat->its[i]];
    for(int j = i + 1; j < pat->sz(); j++)
      if(!compatible[it1][wToIt[pat->its[j]]])
        return false;
  }
  return true;
};

bool Pricing::IsFeasible(Pattern* pat){
  assert(!pat->needUpdate);
  if(pat->getCap() > capAccepted)
    return false;
  pat->val = Pattern::GetDefaultValue();
  FOR(i, SZ(pat->its)){
    int f = 1;
    while(i + 1 < SZ(pat->its) && pat->its[i] == pat->its[i + 1]){
      f++;
      i++;
    }
    int it = wToIt[pat->its[i]];
    pat->val -= v[it] * f;
    if(freq[it] < f)
      return false;
  }
  return true;
}

bool Pricing::CanBeAdd(Pattern* pat){
  assert(!pat->needUpdate);
  return IsFeasible(pat) && IsCompatible(pat);
}

/* must only be called after calling updateInstance at least once */
vector<Pattern*> Pricing::Solve(bool& opt, bool& feasible) {
  cnt++;
  assert(SZ(freq) == SZ(weights));
  assert(SZ(seenSR) == SZ(vSR));

  vector<Pattern*> ret;

  ret = getFromPool();
  if(!ret.empty())
    return ret;

  its.clear();
  FOR(it, SZ(freq))
    if(freq[it] > 0){
      int qt = min(W / weights[it], (!unique && compatible[it][it]) ? freq[it] : 1);
      FOR(q, qt)
        its.push_back(it);
    }

  vi prior(SZ(freq));
  FOR(i, SZ(setsSR))
    for(IM im : setsSR[i])
      prior[im.i] |= (-vSR[i] > eps13) ? 2 : 0;

  for(int it : conflictIts)
    prior[it] |= 1;

  stable_sort(all(its), [&](int i1, int i2){
    return prior[i1] != prior[i2] ? prior[i1] < prior[i2]: weights[i1] < weights[i2];
  });
  ret = DP(opt, feasible);

  if(opt)
    return {};

  assert(!ret.empty());
  generatedCols += SZ(ret);
  return ret;
}


ll Pricing::UpdateVal(Pattern* pat){
  fill(all(seenSR), 0);
  pat->val = Pattern::GetDefaultValue();
  for(int w : pat->its)
    addItVal(pat, wToIt[w]);
  return pat->val;
}
