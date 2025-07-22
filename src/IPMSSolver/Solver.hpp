#ifndef SOLVER_H
#define SOLVER_H

#include "GTHeuristic.hpp"

class Solver{
  Instance* instance;
  Model* model;
  BranchManager* bm;
  RFHeuristic* relaxAndFix;
  GTHeuristic* gtHeuristic;
  Solution*& bestSol;
  
  const ld &L0;
  const int &bestInteger, &L1;
  int &L2;
  bool &optimal, &bounded;
  int iterationsBP = 0;

  bool solveRoot();
  bool relax(ll& curObj, bool print);
  Solution* Graham();
  void branchAndPrice();

  public:
  Solver(Instance* instance);
  ~Solver();
  void Solve();
};

#endif
