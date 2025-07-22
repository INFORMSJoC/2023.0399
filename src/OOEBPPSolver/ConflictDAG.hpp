#ifndef CONFLICTDAG_H
#define CONFLICTDAG_H

#include "Pattern.hpp"

template<class Node>
struct ConflictDAG{
  struct Vertex{
    set<int> childs, confl;
  };

  //map<int, Vertex> adj;
  set<ii> baseConfl, allConfl;
  vector<vector<char>>& compatible;
  vector<Node>& st;
  vector<ii> partialDiff;

  ConflictDAG(vector<vector<char>>& compatible, vector<Node>& st):
     compatible(compatible), st(st){}

  vector<ii> AddMerge(int i, int j, int k){
    //adj[i].childs.insert(k);
    //adj[j].childs.insert(k);
    partialDiff.clear();
    FOR(h, SZ(compatible[i]))
      if(!compatible[h][i] && compatible[h][k])
        setConflict(h, k);

    FOR(h, SZ(compatible[j]))
      if(!compatible[h][j] && compatible[h][k])
        setConflict(h, k);
    //propagate(i);
    //propagate(j);
    return partialDiff;
  }

  vector<ii> AddConflict(int i, int j){
    partialDiff.clear();
    baseConfl.insert(minmax(i, j));
    setConflict(i, j);
    //propagate(i);
    //propagate(j);
    return partialDiff;
  }

  vector<ii> ReCalc(){
    clearList();

    for(auto nd : st){
      if(nd.onEnd)
        continue;
      if(nd.state)
        AddConflict(nd.it1, nd.it2);
      else
        AddMerge(nd.it1, nd.it2, nd.it);
    }

    return vector<ii>(all(allConfl));
  }

  private:
  void setConflict(int i, int j){
    //adj[i].confl.insert(j);
    //adj[j].confl.insert(i);
    allConfl.insert(minmax(i, j));
    partialDiff.push_back(minmax(i, j));
    compatible[i][j] = compatible[j][i] = false;
  }

  /*
  void propagate(int i){
    for(int k : adj[i].childs){
      bool any = false;
      for(int h : adj[i].confl)
        if(!adj[k].confl.count(h)){
          setConflict(k, h);
          any = true;
        }
      if(any)
        propagate(k);
    }
  }*/

  void clearList(){
    /*for(auto pr : adj)
      for(int j : pr.second.confl)
        compatible[pr.first][j] = true;
    adj.clear();*/
    baseConfl.clear();
    allConfl.clear();
    FOR(i, SZ(compatible))
      FOR(j, SZ(compatible))
        compatible[i][j] = true;
    //assert(compatible[i][j]);
  }
};

#endif
