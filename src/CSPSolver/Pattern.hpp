#ifndef PATTERN_H
#define PATTERN_H

#include "Util.hpp"

class Pattern{
  public:
  
  struct lessPatternEQ{
    bool operator() (Pattern* const& b1, Pattern* const& b2) const{
      if(b1->cap != b2->cap)
        return b1->cap > b2->cap;
      assert(!b1->needUpdate);
      assert(!b2->needUpdate);
      return b1->its > b2->its;
    }
  };

  vi its;
  bool needUpdate;
  int cap, id;
  ll val, binId;
  inline static ll lastId = -1;
  inline static vi dels;
  inline static set<Pattern*, lessPatternEQ> allPats;
  inline static vector<Pattern*> idToPat;
  inline static ll dfval = 1;

  static ll GetDefaultValue(){
    return dfval * BASIS;
  }

  inline Pattern(const Pattern& b): its(b.its), needUpdate(b.needUpdate),
    cap(b.cap), id(-1), val(b.val){

    this->binId = ++lastId;
    dels.push_back(false);
  }

  inline Pattern(int cur = -1): cap(0), id(-1), val(GetDefaultValue()){
    this->binId = ++lastId;
    dels.push_back(false);
    if (cur != -1)
      add(cur);
    needUpdate = false;
  }

  inline void add(int weight){
    this->its.push_back(weight);
    this->cap += weight;
    needUpdate = true;
  }

  inline void rem(int weight){
    auto iter = find(all(its), weight);
    assert(iter != its.end());
    its.erase(iter);
    this->cap -= weight;
    needUpdate = true;
  }

  inline void update(){
    if(needUpdate){
      sort(rall(its));
      needUpdate = false;
    }
  }

  inline int residue(){
    return binCapacity - cap;
  }

  inline int sz(){
    return its.size();
  }

  inline bool canPack(int weight){
    return cap + weight <= binCapacity;
  }

  static void updateId(Pattern*& pat){
    pat->update();
    auto iter = Pattern::allPats.find(pat);
    if(iter != Pattern::allPats.end()){
      if((*iter)->id == pat->id)
        return;
      (*iter)->val = pat->val;
      delete pat;
      pat = *iter;
    }else{
      pat->id = Pattern::idToPat.size();
      Pattern::idToPat.push_back(pat);
      Pattern::allPats.insert(pat);
    }
  }

  static Pattern* getTracked(Pattern* unTracked){
    Pattern* pat = new Pattern(*unTracked);
    Pattern::updateId(pat);
    return pat;
  }

  ~Pattern(){
    dels[binId] = true;
  }
};

struct lessPatternID{
  bool operator() (Pattern* const& b1, Pattern* const& b2) const{
    if(b1->cap != b2->cap)
      return b1->cap > b2->cap;
    assert(!b1->needUpdate);
    assert(!b2->needUpdate);
    return b1->its > b2->its || (b1->binId < b2->binId && b1->its == b2->its);
  }
};

#endif
