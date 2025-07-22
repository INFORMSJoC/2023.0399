#include "Solution.hpp"

Solution::Solution(Instance* instance){
  this->instance = instance;
  this->freq = instance->freq.data();
}

Solution::Solution(const Solution& sol){
  this->instance = sol.instance;
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
  int n = instance->NT;
  vi freqCnt(n, 0);

  for(Pattern* bin : sol){
    int cap = 0;
    assert(!bin->needUpdate);
    for(int id : bin->its){
      assert(id < instance->NT);
      cap += Pattern::items[id].w;
      freqCnt[id]++;
    }
    Item first = Pattern::items[bin->its[0]];
    assert(cap - first.w + first.wi <= binCapacity);
  }
  FOR(i, n)
    assert(freqCnt[i] == freq[i]);
}

void Solution::PrintSolution(string pathOutSol, const StatusSolution& statusSol){
  Check();
  ofstream outSolStream;
  outSolStream.open(pathOutSol);
  ll sumSize = 0, totalFreeSpace = 0;
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
    totalFreeSpace += bin->residue();
    outSolStream << "Bin " << ++i << ": " << bin->residue() << " of_free_capacity"<< endl;
    outSolStream << "Amount_of_items: " << bin->sz() << endl;
    outSolStream << "Items: " << endl;
    sort(all(bin->its), [&] (int it1, int it2){return instance->inv[Pattern::items[it1].p] < instance->inv[Pattern::items[it2].p];});
    for(int id : bin->its){
      ii pr = instance->inv[Pattern::items[id].p];
      outSolStream << "Size: " << Pattern::items[id].w << " Prior: " << pr.first << " Id: " << pr.second << endl;
      sumSize += Pattern::items[id].w;
    }
    outSolStream << endl;
  }
  outSolStream << "Total_size: " << sumSize << endl;
  outSolStream << "Total_free_space: " << totalFreeSpace << endl << endl;
  outSolStream.close();

  string commandSystem = instance->pathCheckSol + " " + instance->pathOutSol + " " + instance->pathInstance;
  instance->resultCheckSol |= system(commandSystem.c_str());
  assert(!instance->resultCheckSol);
}

vector<Pattern*> Solution::GetBins(){
  return sol;
}

Solution* Solution::GetBestInitial(Instance* instance){
  Solution* sol = new Solution(instance);
  sol->BFD();
  return sol;
}

void Solution::BFD(){
  REV(i, instance->NT)
    FOR(j, freq[i])
      PackingBF(i);
  Check();
}

void Solution::Unmerge(multimap<int, ii> mergedIts, bool Check){
  for(auto bin : sol){
    bool original = false;
    while(!original){
      original = true;
      vi pIts = bin->its;
      for(int id : pIts){
        auto iter = mergedIts.find(id);
        if(iter != mergedIts.end()){
          ii pr = iter->second;
          mergedIts.erase(iter);
          bin->rem(id);
          bin->add(pr.first);
          bin->add(pr.second);
          original = false;
        }
      }
    }
  }

  for(auto bin : sol)
    for(int id : bin->its)
      assert(!mergedIts.count(id) && id < instance->NT);

  vector<Pattern*> toAdd(all(sol));
  sol.clear();
  for(auto bin : toAdd)
    InsertBin(bin);

  if(Check)
    this->Check();
}
