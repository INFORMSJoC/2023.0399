#include "Pricing.hpp"

#define BEST_THAN(val1, val2) ((val1) < (val2))
#define WORST_EQUAL_THAN(val1, val2) ((val1) >= (val2))
#define IMPROVE(val) BEST_THAN(val, -epsM + 1)
#define VALID(val) (item.w <= b ? (val) : INFMAX)

Pricing::Pricing(CPM* cpm, vi& freq, vi& conflictIts, vll& v, vector<vc>& compatible, bool& multiColor):
    W(binCapacity), items(Pattern::items), freq(freq), conflictIts(conflictIts), v(v),
    vSR(cpm->dualSR), itsToSR(cpm->itsToSR), setsSR(cpm->setsSR), compatible(compatible),
    seenSR(cpm->seenSR), seenColors(Q, 0), residueAccept(binCapacity - 1), multiColor(multiColor){}

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
  FOR(i, pat->sz())
    for(int j = i + 1; j < pat->sz(); j++)
      isGood &= compatible[pat->its[i]][pat->its[j]];

  if(!isGood)
    return false;

  int residue = binCapacity;
  UpdateVal(pat);
  for(int it : pat->its)
    residue -= items[it].w;

  assert(val2 == INFLL || pat->val == val2);
  assert(residue == pat->residue() && residue >= 0);

  return true;
}

void Pricing::recoverBest(int i, int b, int c, bool u, Pattern* pat){
  assert(dp[c][i][b][u] != INFMAX);
  if(i == 0 || b == 0 || (c == 0 && !u)){
    assert(BEST_THAN(pat->val, incumb->val));
    assert(pat->residue() <= residueAccept);
    if(!multiColor || IsColorFeasible(pat)){
      delete incumb;
      incumb = new Pattern(*pat);
      incumb->update();
    }
    return;
  }

  trysDP--;
  if(trysDP < 0)
    return;
  
  int it = its[i];
  Item item = items[it];
  int b2 = b - item.w;
  bool eq = item.c[0] == items[its[i - 1]].c[0];

  if(!u && multiColor){
    for(int id : pat->its){
      for(int i = 1; i < SZ(items[id].c); i++){
        if(items[id].c[i] == item.c[0]){
          c--;
          u = true;
          assert(c >= 0);
          RecursiveDP(i, b, c, u);
          goodAll = false;
          break;
        }
      }
      if(u)
        break;
    }
  }

  bool swapOrder = false;
  int iu = i - 1;
  ll dpRight = dp[c][iu][b][eq && u] + pat->val;
  if(b2 >= 0){
    ll dpLeft = dp[c - !u][i - 1][b2][eq] + pat->val - v[it];
    if(BEST_THAN(dpLeft, incumb->val)){
      if(BEST_THAN(dpRight, dpLeft)){
        recoverBest(iu, b, c, eq && u, pat);
        swapOrder = true;
        if(trysDP < 0)
          return;
      }

      if(canPack(pat, it)){
        addIt(pat, it);
        if(BEST_THAN(dp[c - !u][i - 1][b2][eq] + pat->val, incumb->val)){
          recoverBest(i - 1, b2, c - !u, eq, pat);
          if(trysDP < 0)
            return;
        }
        remIt(pat, it);
      }
    }
  }

  if(!swapOrder and BEST_THAN(dpRight, incumb->val))
    recoverBest(iu, b, c, eq && u, pat);
}

void Pricing::recoverDPLexic(int i, int b, int c, bool u, Pattern* pat){
  assert(dp[c][i][b][u] != INFMAX);

  if(i == 0 || b == 0 || (c == 0 && !u)){
    assert(IMPROVE(pat->val));
    // assert(pat->residue() <= residueAccept);
    if(IsColorFeasible(pat)){
      Pattern* bin = new Pattern(*pat);
      bin->update();
      retDP.push_back(bin);
      for(int it : pat->its)
        freqRet[it]++;
    }
    return;
  }

  trysDP--;
  if(trysDP < 0 && !retDP.empty())
    return;
  
  int it = its[i];
  Item item = items[it];
  int b2 = b - item.w;
  bool eq = item.c[0] == items[its[i - 1]].c[0];

  if(!u && multiColor){
    for(int id : pat->its){
      for(int i = 1; i < SZ(items[id].c); i++){
        if(items[id].c[i] == item.c[0]){
          c--;
          u = true;
          assert(c >= 0);
          RecursiveDP(i, b, c, u);
          goodAll = false;
          break;
        }
      }
      if(u)
        break;
    }
  }

  if(b2 >= 0 && IMPROVE(dp[c - !u][i - 1][b2][eq] + pat->val - v[it])){
    if(canPack(pat, it)){
      addIt(pat, it);
      if(IMPROVE(dp[c - !u][i - 1][b2][eq] + pat->val))
        recoverDPLexic(i - 1, b2, c - !u, eq, pat);
      else
        goodAll = false;
      remIt(pat, it);
    }else
      goodAll = false;
  }
  for(int it : pat->its)
    if(freqRet[it] > 2 * mxFreqRet)
      return;
  if(IMPROVE(dp[c][i - 1][b][eq && u] + pat->val))
    recoverDPLexic(i - 1, b, c, eq && u, pat);
}

bool Pricing::updateDPsize(int N, ll INFMAX){
  if (N > SZ(dp)){
    ll dpSize = N * (W + 1) * (C + 1) * 2;
    if(dpSize * 8 >= 4'500'000'000){
      cerr << "N: " << N << " - W: " << W << " - DP size: " << dpSize << endl;
      throw runtime_error("Memory limit reached in DP");
    }
    dp.assign(C + 1, vector<vector<pll>>(N, vector<pll>(W + 1, pll(INFMAX, INFMAX))));
    return true;
  }
  return false;
}

ll Pricing::RecursiveDP(int i, int b, int c, bool u){
  if(dp[c][i][b][u] != INFMAX)
    return dp[c][i][b][u];

  //if(i == 0 || b == 0 || (c == 0 && !u))
  //  return dp[c][i][b][u] = 0;

  int it = its[i];
  Item item = items[it];
  bool eq = item.c[0] == items[its[i - 1]].c[0]; 
  return dp[c][i][b][u] = min(RecursiveDP(i - 1, b, c, eq && u),
                              VALID(RecursiveDP(i - 1, b - item.w, c - !u, eq) - v[it]));
}

vector<Pattern*> Pricing::DP(bool& opt, bool& feasible){
  cntDP++;
  int N = its.size();

  long double acumAll = 0;
  FOR(i, SZ(its))
    acumAll += v[its[i]];
  assert(2 * acumAll - 1 < INFMAX);

  its.insert(its.begin(), SZ(freq));
  items.push_back(Item(-1, -1));

  if(!updateDPsize(N + 1, INFMAX))
    FOR(c, C + 1)
      FOR(i, N + 1)
        FOR(b, W + 1)
          dp[c][i][b] = pll(INFMAX, INFMAX);

  FOR(c, C + 1)
    FOR(b, W + 1)
        dp[c][0][b] = pll(0, 0);

  FOR(c, C + 1)
    FOR(i, N + 1)
      dp[c][i][0] = pll(0, 0);

  FORI(i, N){
    int it = its[i];
    Item item = items[it];
    int wit = item.w;
    bool eq = item.c[0] == items[its[i - 1]].c[0]; 
    ll vit = v[it];
    FOR(b, W + 1){
      dp[0][i][b] = pll(0, min(dp[0][i - 1][b][eq],
                                VALID(dp[0][i - 1][b - wit][eq] - vit)));
    }
  }
  
  FORI(c, C)
    FORI(i, N){
      int it = its[i];
      Item item = items[it];
      int wit = item.w;
      bool eq = item.c[0] == items[its[i - 1]].c[0]; 
      ll vit = v[it];
      FOR(b, W + 1){
        FOR(u, 2)
          dp[c][i][b][u] = min(dp[c][i - 1][b][eq && u],
                                    VALID(dp[c - !u][i - 1][b - wit][eq] - vit));
      }
    }
  //RecursiveDP(N, W, C, false);

  opt = !IMPROVE(dp[C][N][W][false] + Pattern::GetDefaultValue());
  if(opt){
    feasible = dp[C][N][W][false] + Pattern::GetDefaultValue() >= 0;
    items.pop_back();
    return {};
  }

  Pattern* temp = new Pattern();
  trysDP = N * W * C/ 20;
  goodAll = true;
  retDP.clear();
  fill(all(seenSR), 0);
  freqRet.assign(SZ(freq), 0);

  recoverDPLexic(N, W, C, false, temp);

  if(!retDP.empty()){
    FOR(i, SZ(retDP))
      if(BEST_THAN(retDP[i]->val, retDP[0]->val))
        swap(retDP[i], retDP[0]);
    incumb = new Pattern(*retDP[0]);
    assert(IMPROVE(incumb->val));

    trysDP = N * W / 50;
    recoverBest(N, W, C, false, temp);
    if(BEST_THAN(incumb->val, retDP[0]->val))
      retDP.push_back(new Pattern(*incumb));
    delete incumb;
  }
  delete temp;

  if(retDP.empty()){
    if(goodAll){
      cout << dp[C][N][W][false] + Pattern::GetDefaultValue() << endl;
      // cout << setprecision(16) << dp[N][W] <<  " - N: " << N << " - cnt: "  << cntDP << endl;
      // int i = N, w = W;
      // ll val = Pattern::GetDefaultValue();
      // while(i > 0 && w >= 0){
      //   int it = its[i];
      //   if(weights[it] <= w && fabs(dp[i][w] - (dp[i - 1][w - weights[it]] - v[it])) < EPSP){
      //     val -= v[it];
      //     w -= weights[it];
      //     cout << it << " - v: " << v[it] << " - w: " << weights[it] << " - V: "  <<  val << " - W: " << w  << " - i: " << i << endl;
      //   }
      //   i--;
      // }
      throw runtime_error("Nao achou a coluna boa");
    }
    opt = true;
  }
  items.pop_back();
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

bool Pricing::IsColorFeasible(Pattern* pat){
  int q = 0;
  for(int it : pat->its){
    Item item = items[it];
    for(int c : item.c)
      if(!seenColors[c])
        q++, seenColors[c]++;
  }
  for(int it : pat->its)
    for(int c : items[it].c)
      seenColors[c] = 0;
  return q <= C;
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
    return items[i1] < items[i2];
  });
  FOR(i, SZ(its) - 1)
    assert(items[its[i]].c[0] <= items[its[i + 1]].c[0]);

  ret = DP(opt, feasible);

  if(opt)
    return {};

  assert(!ret.empty());
  generatedCols += SZ(ret);
  return ret;
}

// ll Pricing::MinReducedCost(){
//   Pos::dp = &dp;

//   if(SZ(its) <= 1)
//     return -INF;

//   assert(its[0] == SZ(freq));
//   int n = its.size() - 1;
//   priority_queue<Pos> q;
//   q.push({n, new Pattern(), vi(SZ(seenSR))});

//   Pattern* best = new Pattern();
//   best->val = -EPS13;
//   int cntTrys = 0, mxCnt = SZ(its) * binCapacity / (SZ(seenSR) + 1) / 20;
//   ll rcBest = best->val;

//   while(!q.empty()){
//     cntTrys++;
//     Pos pos = q.top(); q.pop();
//     if(cntTrys >= mxCnt || WORST_EQUAL_THAN(pos.val(), best->val)){
//       rcBest = min(rcBest, pos.val());
//       delete pos.p;
//       continue;
//     }

//     Pattern* pat = pos.p;
//     int i = pos.i;
//     if(i <= 0 || pat->residue() == 0){
//       assert(pat->residue() <= residueAccept);
//       assert(i >= 0 && BEST_THAN(pos.val(), best->val));
//       pat->update();
//       swap(best, pat);
//       delete pat;
//       continue;
//     }
//     int it = its[i], w2 = pat->residue() - weights[it];
//     if(w2 >= 0){
//       ll val2 = pat->val - v[it];
//       if(BEST_THAN(dp[i - 1][w2] + val2, best->val)){
//         if(canPack(pat, it)){
//           seenSR = pos.seen;
//           ll oldVal = pat->val;
//           addIt(pat, it);
//           if(BEST_THAN(dp[i - 1][w2] + pat->val, best->val))
//             q.push({i - 1, new Pattern(*pat), seenSR});
//           remIt(pat, it, oldVal);
//         }
//       }
//     }

//     while(i > 0 && weights[it] == weights[its[i]])
//       i--;
//     if(BEST_THAN(dp[i][pat->residue()] + pat->val, best->val))
//       q.push({i, pat, pos.seen});
//     else
//       delete pat;
//   }
//   if(best->sz())
//     isValid(best, best->val);
//   rcBest = min(rcBest, best->val);
//   delete best;

//   return rcBest;
// }

ll Pricing::UpdateVal(Pattern* pat){
  fill(all(seenSR), 0);
  pat->val = Pattern::GetDefaultValue();
  for(int it : pat->its)
    addItVal(pat, it);
  return pat->val;
}
