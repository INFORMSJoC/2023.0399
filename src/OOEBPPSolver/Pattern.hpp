#ifndef PATTERN_H
#define PATTERN_H

#include "Util.hpp"

struct Item{
  int w, wi, p;
  Item (int w, int wi, int p): w(w), wi(wi), p(p){}

  bool operator < (const Item& it) const{
    return iii(w, wi, p) < iii(it.w, it.wi, it.p);
  }

  bool operator > (const Item& it) const{
    return iii(w, wi, p) > iii(it.w, it.wi, it.p);
  }

  bool operator == (const Item& it) const{
    return iii(w, wi, p) == iii(it.w, it.wi, it.p);
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

  vi its;
  bool needUpdate;
  int cap, id;
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

  inline Pattern(const Pattern& b): its(b.its), needUpdate(b.needUpdate),
    cap(b.cap), id(-1), val(b.val){

    this->binId = ++lastId;
    dels.push_back(false);
  }

  inline Pattern(): cap(0), id(-1), val(GetDefaultValue()){
    this->binId = ++lastId;
    dels.push_back(false);
    needUpdate = false;
  }

  inline void add(int id){
    this->cap += items[id].w;
    this->its.push_back(id);
    needUpdate = true;
  }

  inline void rem(int id){
    auto iter = find(all(its), id);
    assert(iter != its.end());
    its.erase(iter);
    this->cap -= items[id].w;
    needUpdate = true;
  }

  inline void update(){
    if(needUpdate){
      if(!its.empty()){
        FOR(i, SZ(its))
          if(items[its[0]].p < items[its[i]].p)
            swap(its[0], its[i]);
        sort(next(its.begin()), its.end());
      }
      needUpdate = false;
    }
  }

  inline int residue(){
    assert(!its.empty());
    update();
    return binCapacity - (cap - items[its[0]].w + items[its[0]].wi);
  }

  inline int sz(){
    return its.size();
  }

  inline bool canPack(int id){
    Item item = items[id];
    if(its.empty() || cap + items[id].w <= binCapacity)
      return true;
    update();
    Item first = items[its[0]];
    return first.p < item.p ? (cap + item.wi <= binCapacity) :
      ((cap + item.w - first.w + first.wi) <= binCapacity);
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
