#include "RoundingHeuristic.hpp"

RoundingHeuristic::RoundingHeuristic(int& M, vi& weights, vi& freq,
  vi& ids, vi& wToIt, vector<Pattern*>& idToPat, vd& primal): M(M),
  weights(weights), freq(freq), ids(ids), wToIt(wToIt), idToPat(idToPat),
  primal(primal){}

RoundingHeuristic::~RoundingHeuristic(){}

const ld accepted = 0.6;

vector<pair<ld, vi>> RoundingHeuristic::recover(){
  vector<pair<ld, vi>> bins;
  REV(i, M){
    if(primal[i] > accepted && ids[i] >= 0){
      int f = floor(primal[i] + EPS9);
      Pattern* pat = idToPat[ids[i]];
      vi bin;
      for(int w : pat->its)
        bin.push_back(wToIt[w]);

      FOR(j, f)
        bins.push_back({1.0, bin});
      if(primal[i] > f + accepted)
        bins.push_back({primal[i] - f, bin});
    }
  }
  sort(rall(bins));
  return bins;
}

Solution* RoundingHeuristic::roundBFD(vector<pair<ld, vi>>& fratBins, Solution* partial, Instance* instance, int target){
  Solution* sol = new Solution(*partial);

  vi freq2 = freq;
  int SA = round((target - instance->L0) * binCapacity);

  for(auto& pr : fratBins){
    vi& bin = pr.second;
    int RATemp = binCapacity, unused = 0;

    FOR(i, SZ(bin)){
      int f = 1;
      int it = bin[i];
      while(i + 1 < SZ(bin) && it == bin[i + 1]){
        f++;
        i++;
      }
      unused += freq2[it];
      RATemp -= min(freq2[it], f) * weights[it];
    }

    if(unused >= 2 && RATemp <= SA){
      SA -= RATemp;
      Pattern* pat = new Pattern();
      for(int it : bin)
        if(freq2[it] > 0){
          pat->add(weights[it]);
          freq2[it]--;
        }
      sol->InsertBin(pat);
    }
  }

  vi rem;
  FOR(it, SZ(freq2))
    while(freq2[it] > 0){
      rem.push_back(it);
      freq2[it]--;
    }

  sort(all(rem), [weights = weights](int i1, int i2){
    return weights[i1] > weights[i2];
  });
  for(int it : rem)
    sol->PackingBFD(weights[it]);
    
  return sol;
}

bool RoundingHeuristic::Run(Solution*& bestSol, Solution* partial, const multimap<int, ii>& mergedIts, int target, int integers){
  vector<pair<ld, vi>> fratBins = recover();
  Solution* sol = roundBFD(fratBins, partial, bestSol->instance, target);

  sol->Unmerge(mergedIts, true);
  assert(target - integers > 0 || integers == sol->Cost());
  if(bestSol->Cost() > sol->Cost())
    swap(bestSol, sol);
  delete sol;
  return bestSol->Cost() <= target;
}
