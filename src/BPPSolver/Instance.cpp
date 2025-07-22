#include "Instance.hpp"

int binCapacity;

Instance::Instance(): L0(-1), NT(0), L1(-1) {}


Instance* Instance::ReadInstance(string pathInstance){
  Instance* instance = new Instance();
  int &numItems = instance->NT;
  ifstream instanceStream;

  instanceStream.open(pathInstance);

  if (!instanceStream.good()){
    throw new invalid_argument("Not found instance file");
  }

  instanceStream >> numItems >> binCapacity;
  vw& items = Pattern::items;
  instance->freq.assign(numItems, 1);

  FOR(i, numItems){
    int w;
    instanceStream >> w;
    items.push_back({w});
  }
  instanceStream.close();

  sort(all(items));
  ll sumSize = 0;
  for(auto it : items)
    sumSize += it.w;
  instance->L0 = 1.0L * sumSize / binCapacity;
  instance->L1 = (sumSize / binCapacity) + (sumSize % binCapacity != 0);

  return instance;
}
