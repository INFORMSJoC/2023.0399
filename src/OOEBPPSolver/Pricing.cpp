#include "Pricing.hpp"

#define BEST_THAN(val1, val2) ((val1) < (val2))
#define WORST_EQUAL_THAN(val1, val2) ((val1) >= (val2))
#define IMPROVE(val) BEST_THAN(val, -epsM + 1)

Pricing::Pricing(CPM* cpm, vi& freq, vi& canMid, vi& canEnd, vi& conflictIts, vll& v, vector<vc>& compatible):
    W(binCapacity), items(Pattern::items), freq(freq), canMid(canMid), canEnd(canEnd), 
    conflictIts(conflictIts), v(v), vSR(cpm->dualSR), itsToSR(cpm->itsToSR),
    setsSR(cpm->setsSR), compatible(compatible), seenSR(cpm->seenSR),
    residueAccept(binCapacity - 1){}

Pricing::~Pricing(){
  for(Pattern* pat : pool)
    delete pat;
}

vector<Pattern*> Pricing::getFromPool(){
  vector<Pattern*> bests;

  for(auto pat : pool){
    if(pat->residue() > residueAccept)
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

    for(Pattern* pat : bests){
      assert(isValid(pat, pat->val));
      bool can = true;
      for(int it : pat->its)
        can &= freqRet[it] < mxFreqRet;
      if(can){
        for(int it : pat->its)
          freqRet[it]++;
        Pattern::updateId(pat);
        ret.push_back(pat);
      }else if(del)
        delete pat;
    }
    ret.swap(bests);
  }
}

bool Pricing::isValid(Pattern* pat, ll val2){
  bool isGood = true;
  pat->update();
  FOR(i, pat->sz()){
    if(i == 0)
      isGood &= canEnd[pat->its[i]];
    else
      isGood &= canMid[pat->its[i]];
    for(int j = i + 1; j < pat->sz(); j++)
      isGood &= compatible[pat->its[i]][pat->its[j]];
  }

  if(!isGood)
    return false;

  int residue = binCapacity;
  UpdateVal(pat);
  for(int it : pat->its)
    residue -= items[it].w;

  pat->update();
  residue += items[pat->its[0]].w - items[pat->its[0]].wi;

  assert(val2 == INFLL || pat->val == val2);
  assert(residue == pat->residue() && residue >= 0);

  return true;
}

void Pricing::recoverBest(int i, Pattern* pat){
  if(i <= 0 || pat->residue() == 0){
    assert(pat->residue() <= residueAccept);
    assert(BEST_THAN(pat->val, incumb->val));
    delete incumb;
    incumb = new Pattern(*pat);
    incumb->update();
    return;
  }

  trysDP--;
  if(trysDP < 0)
    return;

  int it = its[i], iu = i - 1, w2 = pat->residue() - items[it].w;

  bool swapOrder = false;
  const ll dpRight = dp[iu][pat->residue()] + pat->val;
  if(w2 >= 0 && canMid[it]){
    const ll dpLeft = dp[i - 1][w2] + pat->val - v[it];
    if(BEST_THAN(dpLeft, incumb->val)){

      if(BEST_THAN(dpRight, dpLeft)){
        recoverBest(iu, pat);
        swapOrder = true;
        if(trysDP < 0)
          return;
      }

      if(canPack(pat, it)){
        addIt(pat, it);
        if(BEST_THAN(dp[i - 1][w2] + pat->val, incumb->val)){
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
  if(i <= 0 || pat->residue() == 0){
    assert(pat->residue() <= residueAccept);
    Pattern* bin = new Pattern(*pat);
    bin->update();
    retDP.push_back(bin);
    for(int it : pat->its)
      freqRet[it]++;
    return;
  }

  trysDP--;
  if(trysDP < 0 && !retDP.empty())
    return;

  int it = its[i], w2 = pat->residue() - items[it].w;

  if(w2 >= 0 && canMid[it]){
    ll dpTo = dp[i - 1][w2] + pat->val - v[it];
    if(IMPROVE(dpTo)){
      if(canPack(pat, it)){
        addIt(pat, it);
        if(IMPROVE(dp[i - 1][w2] + pat->val))
          recoverDPLexic(i - 1, pat);
        else
          goodAll = false;
        remIt(pat, it);
      }else
        goodAll = false;
    }
  }

  for(int it : pat->its)
    if(freqRet[it] > 2 * mxFreqRet)
      return;

  i--;
  if(IMPROVE(dp[i][pat->residue()] + pat->val))
    recoverDPLexic(i, pat);
}

bool Pricing::updateDPsize(int N, ll INFMAX){
  if (N > SZ(dp)){
    ll dpSize = N  * (W + 1);
    if(dpSize * 8 >= 4'500'000'000){
      cerr << "N: " << N << " - W: " << W << " - DP size: " << dpSize << endl;
      throw runtime_error("Memory limit reached in DP");
    }
    dp.assign(N, vll(W + 1, INFMAX));
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

  FOR(w, residueAccept + 1)
    dp[0][w] = 0;

  FORI(i, N){
    int it = its[i], wg = items[it].w;
    REV(w, W + 1)
      dp[i][w] = min(dp[i - 1][w], w - wg >= 0 && canMid[it] ? dp[i - 1][w - wg] - v[it] : INFMAX);
  }

  opt = true;
  FORI(i, N)
    opt &= !canEnd[its[i]] || !IMPROVE(dp[i - 1][W - items[its[i]].wi] - v[its[i]] + Pattern::GetDefaultValue());
  if(opt){
    feasible = true;
    FORI(i, N)
      feasible &= !canEnd[its[i]] || (dp[i - 1][W - items[its[i]].wi] - v[its[i]] + Pattern::GetDefaultValue() >= 0);
    return {};
  }

  vector<Pattern*> retDPGlobal;
  goodAll = true;
  FORI(i, N){
    if(canEnd[its[i]] && IMPROVE(dp[i - 1][W - items[its[i]].wi] - v[its[i]] + Pattern::GetDefaultValue())){
      Pattern* temp = new Pattern();
      trysDP = N * W / 20;
      retDP.clear();
      fill(all(seenSR), 0);
      freqRet.assign(SZ(freq), 0);
      addIt(temp, its[i]);
      recoverDPLexic(i - 1, temp);
      if(i - 1 > 0 && !retDP.empty()){
        FOR(i, SZ(retDP))
          if(BEST_THAN(retDP[i]->val, retDP[0]->val))
            swap(retDP[i], retDP[0]);
        incumb = new Pattern(*retDP[0]);
        assert(IMPROVE(incumb->val));

        trysDP = (i - 1) * W / 50;
        recoverBest(i - 1, temp);
        if(BEST_THAN(incumb->val, retDP[0]->val))
          retDP.push_back(new Pattern(*incumb));
        delete incumb;
      }
      delete temp;
      retDPGlobal.insert(retDPGlobal.end(), all(retDP));
    }
  }
  retDP = retDPGlobal;


  if(retDP.empty()){
    if(goodAll){
      // cout << setprecision(16) << DD(dp[N][W]) <<  " - N: " << N << " - cnt: "  << cntDP << endl;
      // int i = N, w = W;
      // ll val = Pattern::GetDefaultValue();
      // while(i > 0 && w >= 0){
      //   int it = its[i];
      //   if(weights[it] <= w && dp[i][w] == (dp[i - 1][w - weights[it]] - v[it]))){
      //     val -= v[it];
      //     w -= weights[it];
      //     cout << it << " - v: " << DD(v[it]) << " - w: " << weights[it] << " - V: "  <<  DD(val) << " - W: " << w  << " - i: " << i << endl;
      //   }
      //   i--;
      // }
      throw runtime_error("Nao achou a coluna boa");
    }
    opt = true;
  }
  for(auto pat : retDP)
    assert(IMPROVE(pat->val));
  dryBests(retDP, true); // can delete patterns with propagation error
  opt |= retDP.empty();
  return retDP;
}

bool Pricing::canPack(Pattern* bin, int it){
  if(!bin->canPack(it))
    return false;
  for(int it2 : bin->its)
    if(!compatible[it][it2])
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
  pat->add(it);
  addItVal(pat, it);
}

void Pricing::remIt(Pattern* pat, int it){
  pat->rem(it);
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

void Pricing::ClearPool(int residue){
  for(auto iter = pool.begin(); iter != pool.end();){
    Pattern* pat = *iter;
    if(pat->residue() > residue){
      iter = pool.erase(iter);
      Pattern::allPats.erase(pat);
      Pattern::idToPat[pat->id] = NULL;
      delete pat;
    }else
      iter++;
  }
}

bool Pricing::IsCompatible(Pattern* pat){
  pat->update();
  FOR(i, pat->sz()){
    if(i == 0 && !canEnd[pat->its[i]])
      return false;
    if(i > 0 && !canMid[pat->its[i]])
      return false;
    int it1 = pat->its[i];
    for(int j = i + 1; j < pat->sz(); j++)
      if(!compatible[it1][pat->its[j]])
        return false;
  }
  return true;
};

bool Pricing::IsFeasible(Pattern* pat){
  assert(!pat->needUpdate);
  if(pat->residue() > residueAccept)
    return false;
  pat->val = Pattern::GetDefaultValue();
  FOR(i, SZ(pat->its)){
    int f = 1, it = pat->its[i];
    while(i + 1 < SZ(pat->its) && it == pat->its[i + 1]){
      f++;
      i++;
    }
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
  assert(SZ(freq) == SZ(items));
  assert(SZ(seenSR) == SZ(vSR));

  vector<Pattern*> ret;

  ret = getFromPool();
  if(!ret.empty())
    return ret;

  its.clear();
  FOR(it, SZ(freq))
    if(freq[it] > 0){
      assert(freq[it] == 1);
      its.push_back(it);
    }

  stable_sort(all(its), [&](int i1, int i2){
    return items[i1].p < items[i2].p;
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
  for(int it : pat->its)
    addItVal(pat, it);
  return pat->val;
}
