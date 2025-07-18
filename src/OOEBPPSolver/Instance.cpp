#include "Instance.hpp"

int binCapacity;

Instance::Instance(): L0(-1), NT(0), L1(-1) {}

Instance* Instance::ReadInstance(string pathInstance){
  Instance* instance = new Instance();
  int &numItems = instance->NT, lx;
  ifstream instanceStream;

  instanceStream.open(pathInstance);

  if (!instanceStream.good()){
    throw new invalid_argument("Not found instance file");
  }

  instanceStream >> numItems >> binCapacity >> lx;
  vw& items = Pattern::items;
  instance->freq.assign(numItems, 1);

  map<ii, int> mp;
  FOR(i, numItems){
    int c, w, p;
    instanceStream >> c >> w >> p;
    mp[ii(p, c)] = w;
    assert(w <= binCapacity);
  }
  assert(SZ(mp) == numItems);
  instanceStream.close();

  int id = 0;
  for(auto e : mp){
    instance->inv[id] = e.first;
    items.push_back(Item(e.second, 1, id)), id++;
  }

  { 
    vi linWeighths;
    for(Item it : items)
      linWeighths.push_back(it.w);
    ll sumSize = 0;
    std::sort(all(linWeighths));
    int BC1 = binCapacity - 1;
    bool find = false;

    FOR(i, numItems){
      sumSize += linWeighths[i];
      ll curL1 = sumSize / BC1 + (sumSize % BC1 != 0);
      if(curL1 >= (numItems - i - 1)){
        instance->L0 = 1.0L * sumSize / BC1;
        instance->L1 = curL1;
        find = true;
        break;
      }
    }

    assert(find);
  }

  return instance;
}