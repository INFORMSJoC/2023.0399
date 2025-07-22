#ifndef SOL_H
#define SOL_H

#include "Instance.hpp"

struct StatusSolution{
  ld rootObj, rootTime, finalTime, integers = -1, nonZero = -1;
  int L2, Z, where = -2;
  string status;
};

class Solution{
  private:
    vector<Pattern*> sol;
    int *freq;

    void copyBins(const Solution& sol);

    inline void addItem(Pattern* bin, int id){
      RemoveBin(bin);
      bin->add(id);
      InsertBin(bin);
    }

    inline void vectorPackingBF(vi& its){
      for(int id : its)
        PackingBF(id);
    }

  public:
    Instance* instance;

    Solution(Instance* instance);
    Solution(const Solution& sol);
    ~Solution();
    
    vector<Pattern*> GetBins();
    void Check() const;
    static Solution* GetBestInitial(Instance* instance);
    void BFD();
    void Unmerge(multimap<int, ii> mergedIts, bool check);

    void PrintSolution(string pathOutSol, const StatusSolution& statusSol);
    void PrintSolution(const char* s, const StatusSolution& statusSol){
      PrintSolution(string(s), statusSol);
    }

    inline void RemoveBin(Pattern* bin){
      bin->update();
      sol.erase(lower_bound(all(sol), bin, lessPatternID()));
    }

    inline void InsertBin(Pattern* bin){
      assert(bin->id == -1);
      bin->update();
      sol.insert(upper_bound(all(sol), bin, lessPatternID()), bin);
    }

    inline void PackingBF(int id){
      ii best(-1, -1);
      Item item = Pattern::items[id];

      FOR(i, SZ(sol)){
        auto pat = sol[i];
        if(pat->canPack(id) && best.first < pat->cap + item.w)
          best = ii(pat->cap + item.w, i);
      }
      if(best.first != -1){
        addItem(sol[best.second], id);
      }
      else{
        Pattern* pat = new Pattern();
        pat->add(id);
        InsertBin(pat);
      }
    }

    int Cost() const{
      return SZ(sol);
    }

    bool operator < (const Solution& solB) const{
      if(Cost() != solB.Cost())
        return Cost() < solB.Cost();

      FOR(i, SZ(sol))
        if(sol[i]->cap != solB.sol[i]->cap)
          return sol[i]->cap > solB.sol[i]->cap;

      return false;
    }
};

struct GreaterSol{
  bool operator () (const Solution* sol1, const Solution* sol2) const{
    return *sol2 < *sol1;
  }
};

#endif
