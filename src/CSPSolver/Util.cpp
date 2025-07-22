#include "Util.hpp"

unordered_set<string> features;

vi Util::mask(const vi& v, const vi& M){
  vi ret;
  for(int e : v)
    ret.push_back(M[e]);
  return ret;
}
