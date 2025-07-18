#ifndef PATTERN_H
#define PATTERN_H

#include "Util.hpp"

struct Item{
  int w;
  vi c;
  Item (int w, vi c): w(w), c(c){}
  Item (int w, int c): w(w), c({c}){}

  bool operator < (const Item& it) const{
    return iiv(c[0], w, c) < iiv(it.c[0], it.w, it.c);
  }

  bool operator > (const Item& it) const{
    return iiv(c[0], w, c) > iiv(it.c[0], it.w, it.c);
  }

  bool operator == (const Item& it) const{
    return iiv(c[0], w, c) == iiv(it.c[0], it.w, it.c);
  }
};

using vw = vector<Item>;

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

  vi its, colors;
  bool needUpdate;
  int cap, id, nc;
  ll val, binId;
  inline static ll lastId = -1;
  inline static vi dels;
  inline static set<Pattern*, lessPatternEQ> allPats;
  inline static vector<Pattern*> idToPat;
  inline static vw items;
  inline static ll dfval = 1;

  static ll GetDefaultValue(){
    return dfval * BASIS;
  }

  inline Pattern(const Pattern& b): its(b.its), colors(b.colors), 
    needUpdate(b.needUpdate), cap(b.cap), id(-1), nc(b.nc), val(b.val){

    this->binId = ++lastId;
    dels.push_back(false);
  }

  inline Pattern(): colors(Q), cap(0), id(-1), nc(0), val(GetDefaultValue()){
    this->binId = ++lastId;
    dels.push_back(false);
    needUpdate = false;
  }

  inline void add(int id){
    this->cap += items[id].w;
    this->its.push_back(id);
    for(int cl : items[id].c)
      colors[cl]++, nc += colors[cl] == 1;
    needUpdate = true;
  }

  inline void rem(int id){
    auto iter = find(all(its), id);
    assert(iter != its.end());
    its.erase(iter);
    this->cap -= items[id].w;
    for(int cl : items[id].c)
      colors[cl]--, nc -= colors[cl] == 0;
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

  inline bool canPack(int id){
    Item item = items[id];
    if(cap + item.w > binCapacity)
      return false;
    int mc = nc;
    for(int cl : items[id].c)
      mc += colors[cl] == 0;
    return mc <= C;
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
