#include "Solution.hpp"

Solution::Solution(Instance* instance){
  this->instance = instance;
  this->weights = instance->weights.data();
  this->freq = instance->freq.data();
}

Solution::Solution(const Solution& sol){
  this->instance = sol.instance;
  this->weights = sol.weights;
  this->freq = sol.freq;
  copyBins(sol);
}

void Solution::copyBins(const Solution& sol){
  for(Pattern* bin : sol.sol)
    InsertBin(new Pattern(*bin));
}

Solution::~Solution(){
  for(Pattern* bin : sol)
    delete bin;
}

void Solution::Check() const{
  if(instance->L1 == 0)
    return;
  int n = instance->size;
  vi freqCnt(n, 0);

  for(Pattern* bin : sol){
    int cap = 0;
    for(int w : bin->its){
      cap += w;
      int it = lower_bound(weights, weights + n, w) - weights;
      assert(w == weights[it]);
      freqCnt[it]++;
    }
    assert(cap >= binCapacity);
  }
  FOR(i, n)
    assert(freqCnt[i] == freq[i]);
}

void Solution::PrintSolution(string pathOutSol, const StatusSolution& statusSol){
  Check();
  ofstream outSolStream;
  outSolStream.open(pathOutSol);
  ll sumSize = 0;
  outSolStream << setprecision(15);
  outSolStream << "N: " << instance->NT << endl;
  outSolStream << "W: " << binCapacity << endl;
  outSolStream << "L0: " << instance->L0 << endl;
  outSolStream << "L1: " << instance->L1 << endl;
  outSolStream << "Root_Obj: " << statusSol.rootObj << endl;
  outSolStream << "L2: " << statusSol.L2 << endl;
  outSolStream << "Z: " << statusSol.Z << endl;
  outSolStream << "Status: " << statusSol.status << endl;
  outSolStream << "Root_Time: " << statusSol.rootTime << endl;
  outSolStream << "Final_Time: " << statusSol.finalTime << endl;

  outSolStream << endl;

  int i = 0;
  for(Pattern* bin: sol){
    outSolStream << "Bin " << ++i << ": " << bin->getCap() << " of_capacity"<< endl;
    outSolStream << "Amount_of_items: " << bin->sz() << endl;
    outSolStream << "Items: " << endl;
    for(int w : bin->its){
      outSolStream << "Size: " << w << endl;
      sumSize += w;
    }
    outSolStream << endl;
  }
  outSolStream << "Total_size: " << sumSize << endl;
  outSolStream.close();

  string commandSystem = instance->pathCheckSol + " " + instance->pathOutSol + " " + instance->pathInstance;
  instance->resultCheckSol |= system(commandSystem.c_str());
  //assert(!instance->resultCheckSol);
}

vector<Pattern*> Solution::GetBins(){
  return sol;
}

Solution* Solution::GetBestInitial(Instance* instance){
  Solution* sol = new Solution(instance);
  sol->FFD();
  return sol;
}

void Solution::FFD(){
  if(instance->L1 == 0)
    return;
  vi ws;
  REV(i, instance->size)
    FOR(j, freq[i])
      ws.push_back(weights[i]);
  VectorPackingBFD(ws);
  Check();
}

void Solution::Unmerge(multimap<int, ii> mergedIts, bool Check){
  for(auto bin : sol){
    bool original = false;
    while(!original){
      original = true;
      vi pIts = bin->its;
      for(int w : pIts){
        auto iter = mergedIts.find(w);
        if(iter != mergedIts.end()){
          ii pr = iter->second;
          mergedIts.erase(iter);
          bin->rem(w);
          bin->add(pr.first);
          bin->add(pr.second);
          original = false;
        }
      }
    }
  }

  int n = instance->size;
  for(auto bin : sol){
    for(int w : bin->its){
      int it = lower_bound(weights, weights + n, w) - weights;
      assert(w == weights[it]);
      assert(!mergedIts.count(w));
    }
  }

  vector<Pattern*> toAdd(all(sol));
  sol.clear();
  for(auto bin : toAdd)
    InsertBin(bin);

  if(Check)
    this->Check();
}


void Solution::VectorPackingBFD(vector<int>& ws){
  sort(all(ws));
  if(accumulate(all(ws), 0LL) < binCapacity){
    auto bin = sol.back();
    sol.pop_back();
    for(int w : ws)
      bin->add(w);
    InsertBin(bin);
    return;
  }
  
  Pattern* bin = new Pattern();
  while(!ws.empty()){
    int it = 0;
    if(ws[0] + bin->getCap() < binCapacity){
      auto iter = lower_bound(all(ws), binCapacity - bin->getCap());
      if(iter == ws.end())
        iter--;
      it = iter - ws.begin();
    }
    bin->add(ws[it]);
    ws.erase(next(ws.begin(), it));
    
    if(bin->getCap() >= binCapacity){
      ll rem = accumulate(all(ws), 0LL);
      if(rem < binCapacity){
        for(int w : ws)
          bin->add(w);
        ws.clear();
        InsertBin(bin);
      }else{
        InsertBin(bin);
        bin = new Pattern();
      }
    }
  }
  assert(ws.empty());
}