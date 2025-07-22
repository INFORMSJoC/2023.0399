#ifndef DYNAMICDAG_H
#define DYNAMICDAG_H

#include "Pattern.hpp"

class DynamicDAG{
  private:
  int n;
  vvi toUp;

  void extend(int v){
    while(v >= n){
      n++;
      toUp.push_back(vi());
    }
  }
  
  public:
  DynamicDAG(): n(0){}

  bool AddEdge(int v, int p){
    extend(v);
    extend(p);
    if(!count(all(toUp[v]), p)){
      toUp[v].push_back(p);
      return true;
    }
    return false;
  }

  void RemEdge(int v, int p){
    auto it = find(all(toUp[v]), p);
    if(it != toUp[v].end())
      toUp[v].erase(it);
  }

  bool IsAncestor(int v, int p){
    extend(v);
    extend(p);
    set<int> seen;
    stack<int> st;
    st.push(v);
    seen.insert(v);
    while(!st.empty()){
      int v = st.top(); st.pop();
      if(v == p)
        return true;
      for(int u : toUp[v]){
        if(!seen.count(u)){
          seen.insert(u);
          st.push(u);
        }
      }
    }
    return false;
  }
};

#endif
