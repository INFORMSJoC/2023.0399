#include "Pricing.hpp"

#define BEST_THAN(val1, val2) ((val1) < (val2))
#define WORST_EQUAL_THAN(val1, val2) ((val1) >= (val2))
#define IMPROVE(val) BEST_THAN(val, -epsM + 1)

Pricing::Pricing(CPM* cpm, vi& weights, vi& freq, vi& wToIt, vi& conflictIts, vll& v, vector<vc>& compatible):
    W(binCapacity), weights(weights), freq(freq), wToIt(wToIt), conflictIts(conflictIts),
    v(v), vSR(cpm->dualSR), itsToSR(cpm->itsToSR), setsSR(cpm->setsSR), compatible(compatible),
    seenSR(cpm->seenSR), residueAccept(binCapacity - 1){}

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
      for(int w : pat->its)
        can &= freqRet[wToIt[w]] < mxFreqRet;
      if(can){
        for(int w : pat->its)
          freqRet[wToIt[w]]++;
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
  FOR(i, pat->sz())
    for(int j = i + 1; j < pat->sz(); j++)
      isGood &= compatible[wToIt[pat->its[i]]][wToIt[pat->its[j]]];

  if(!isGood)
    return false;

  int residue = binCapacity;
  UpdateVal(pat);
  for(int w : pat->its)
    residue -= weights[wToIt[w]];

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

  int it = its[i], iu = i - 1, w2 = pat->residue() - weights[it];
  while(iu > 0 && weights[its[iu]] == weights[it])
    iu--;

  bool swapOrder = false;
  const ll dpRight = dp[iu][pat->residue()] + pat->val;
  if(w2 >= 0){
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
    for(int w : pat->its)
      freqRet[wToIt[w]]++;
    return;
  }

  trysDP--;
  if(trysDP < 0 && !retDP.empty())
    return;

  int it = its[i], w2 = pat->residue() - weights[it];

  if(w2 >= 0){
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

  for(int w : pat->its)
    if(freqRet[wToIt[w]] > 2 * mxFreqRet)
      return;

  while(i > 0 && weights[its[i]] == weights[it])
    i--;

  if(IMPROVE(dp[i][pat->residue()] + pat->val))
    recoverDPLexic(i, pat);
}

bool Pricing::updateDPsize(int N, ll INFMAX){
  if (N > SZ(dp)){
    ll dpSize = N  * 1LL * (W + 1);
    if(dpSize * 8 >= 25'000'000'000LL){
      cerr << "N: " << N << " - W: " << W << " - DP size: " << dpSize << endl;
      throw runtime_error("Memory limit reached in DP");
    }
    dp.assign(N, vll(W + 1, INFMAX));
    return true;
  }
  return false;
}

vector<Pattern*> Pricing::DP(bool& opt, bool& feasible){
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
    int it = its[i], wg = weights[it];
    REV(w, W + 1)
      dp[i][w] = min(dp[i - 1][w], w - wg >= 0 ? dp[i - 1][w - wg] - v[it] : INFMAX);
  }

  opt = !IMPROVE(dp[N][W] + Pattern::GetDefaultValue());
  if(opt){
    feasible = dp[N][W] + Pattern::GetDefaultValue() >= 0;
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
        retDP.push_back(new Pattern(*incumb));
      delete incumb;
    }
  }else{
    incumb = new Pattern(*temp);
    incumb->val = -epsM + 1;
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
      cout << setprecision(16) << DD(dp[N][W]) <<  " - N: " << N << " - cnt: "  << cntDP << endl;
      int i = N, w = W;
      ll val = Pattern::GetDefaultValue();
      while(i > 0 && w >= 0){
        int it = its[i];
        if(weights[it] <= w && dp[i][w] == (dp[i - 1][w - weights[it]] - v[it])){
          val -= v[it];
          w -= weights[it];
          cout << it << " - v: " << DD(v[it]) << " - w: " << weights[it] << " - V: "  <<  DD(val) << " - W: " << w  << " - i: " << i << endl;
        }
        i--;
      }
      throw runtime_error("Nao achou a coluna boa");
    }
    opt = true;
  }
  for(auto pat : retDP)
    assert(IMPROVE(pat->val));
  dryBests(retDP, true); // can delete patterns with propagation error
  opt |= retDP.empty();
  cntDP += !retDP.empty();
  return retDP;
}

bool Pricing::canPack(Pattern* bin, int it){
  if(!bin->canPack(weights[it]))
    return false;
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
  if(pat->residue() > residueAccept)
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

ll Pricing::MinReducedCost(){
  Pos::dp = &dp;

  if(SZ(its) <= 1)
    return -INFLL;

  assert(its[0] == SZ(freq));
  int n = its.size() - 1;
  priority_queue<Pos> q;
  q.push({n, new Pattern(), vi(SZ(seenSR))});

  Pattern* best = new Pattern();
  best->val = -eps13;
  int cntTrys = 0, mxCnt = SZ(its) * binCapacity / (SZ(seenSR) + 1) / 20;
  ll rcBest = best->val;

  while(!q.empty()){
    cntTrys++;
    Pos pos = q.top(); q.pop();
    if(cntTrys >= mxCnt || WORST_EQUAL_THAN(pos.val(), best->val)){
      rcBest = min(rcBest, pos.val());
      delete pos.p;
      continue;
    }

    Pattern* pat = pos.p;
    int i = pos.i;
    if(i <= 0 || pat->residue() == 0){
      assert(pat->residue() <= residueAccept);
      assert(i >= 0 && BEST_THAN(pos.val(), best->val));
      pat->update();
      swap(best, pat);
      delete pat;
      continue;
    }
    int it = its[i], w2 = pat->residue() - weights[it];
    if(w2 >= 0){
      ll val2 = pat->val - v[it];
      if(BEST_THAN(dp[i - 1][w2] + val2, best->val)){
        if(canPack(pat, it)){
          seenSR = pos.seen;
          addIt(pat, it);
          if(BEST_THAN(dp[i - 1][w2] + pat->val, best->val))
            q.push({i - 1, new Pattern(*pat), seenSR});
          remIt(pat, it);
        }
      }
    }

    while(i > 0 && weights[it] == weights[its[i]])
      i--;
    if(BEST_THAN(dp[i][pat->residue()] + pat->val, best->val))
      q.push({i, pat, pos.seen});
    else
      delete pat;
  }
  if(best->sz())
    isValid(best, best->val);
  rcBest = min(rcBest, best->val);
  delete best;

  return rcBest;
}

ll Pricing::UpdateVal(Pattern* pat){
  fill(all(seenSR), 0);
  pat->val = Pattern::GetDefaultValue();
  for(int w : pat->its)
    addItVal(pat, wToIt[w]);
  return pat->val;
}
