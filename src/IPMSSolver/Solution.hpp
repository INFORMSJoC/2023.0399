#ifndef SOL_H
#define SOL_H

#include "Instance.hpp"

struct StatusSolution{
  ld finalTime;
  int Z;
  string status;
};

class Solution{
  private:
    vector<Pattern*> sol;
    int* weights, *freq;

    void copyBins(const Solution& sol);

    inline void addItem(Pattern* bin, int w){
      RemoveBin(bin);
      bin->add(w);
      InsertBin(bin);
    }

    inline void vectorPackingBFD(vector<int>& items){
      for(int it : items)
        PackingBFD(it);
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

    inline void PackingBFD(int weight){
      auto it = lower_bound(all(sol), weight, [](Pattern* pat, int val){
        return pat->residue() < val;
      });

      if(it != sol.end())
        addItem(*it, weight);
      else
        InsertBin(new Pattern(weight));
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
