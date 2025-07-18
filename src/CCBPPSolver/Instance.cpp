#include "Instance.hpp"

int binCapacity, C, Q;

Instance::Instance(): L0(-1), NT(0), L1(-1) {}

void compute(vw items, ld& L0, int& L1){
  ll sumSize = 0;
  for(auto it : items)
    sumSize += it.w;
  L0 = 1.0L * sumSize / binCapacity;
  L1 = (sumSize / binCapacity) + (sumSize % binCapacity != 0);
}

Instance* Instance::ReadInstance(string pathInstance){
  Instance* instance = new Instance();
  int &numItems = instance->NT;
  ifstream instanceStream;

  instanceStream.open(pathInstance);

  if (!instanceStream.good()){
    throw new invalid_argument("Not found instance file");
  }

  instanceStream >> numItems >> Q >> binCapacity >> C;
  vw& items = Pattern::items;
  instance->freq.assign(numItems, 1);

  FOR(i, numItems){
    int c, w;
    instanceStream >> w >> c;
    items.push_back({w, c});
  }
  instanceStream.close();

  sort(all(items));
  if(C > 1)
    compute(items, instance->L0, instance->L1);
  else{
    vector<vw> byColor(Q);
    for(auto it : items)
      byColor[it.c[0]].push_back(it);
    instance->L1 = 0;
    FOR(i, Q){
      ld L0;
      int L1;
      compute(byColor[i], L0, L1);
      instance->L1 += L1;
    }
    instance->L0 = instance->L1;
  }

  return instance;
}
