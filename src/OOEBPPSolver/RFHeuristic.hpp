#ifndef RF_HEURISTIC
#define RF_HEURISTIC

#include "Model.hpp"

struct NodeRF {
  vi its, patIds;
  vector<Pattern*> partialBins;
  int sz;
  NodeRF(): sz(0){}
};

class RFHeuristic {
  private:
    Model* model;
    Solution*& partial;
    const bool& bounded, &optimal, &insideGT;
    bool& insideRF;
    const int& M, &bestInteger, &integers;
    vi& ids;
    vd& primal;
    const vector<Pattern*>& idToPat;

    void fixPatsVar(vector<NodeRF>& st);
    void fixVars(vector<NodeRF>& st, vi idsToRem);
    void remBackNodeRF(vector<NodeRF>& stRF, int h, bool improveIncumb, bool gt);

  public:
    RFHeuristic(Model* model): model(model), partial(model->partial), bounded(model->bounded),
    optimal(model->optimal), insideGT(model->insideGT), insideRF(model->insideRF), M(model->M),
    bestInteger(model->bestInteger), integers(model->integers), ids(model->ids),
    primal(model->primal), idToPat(model->idToPat){}
    ~RFHeuristic(){}

    bool Run(ll& curObj);
};

#endif
