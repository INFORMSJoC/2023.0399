#ifndef ROUNDINGHEURISTIC_H
#define ROUNDINGHEURISTIC_H

#include "Solution.hpp"

class RoundingHeuristic{
  private:
  const int& M;
  const vw &items;
  const vi &freq, &ids;
  const vector<Pattern*>& idToPat;
  const vd& primal;

  vector<pair<ld, vi>> recover();
  Solution* roundBFD(vector<pair<ld, vi>>& fratBins, Solution* partial, Instance* instance, int target);

  public:
  RoundingHeuristic(int& M, vi& freq, vi& ids, vector<Pattern*>& idToPat, vd& primal);
  ~RoundingHeuristic();
  
  bool Run(Solution*& bestSol, Solution* partial, const multimap<int, ii>& mergedIts, int target, int integers);
};

#endif
