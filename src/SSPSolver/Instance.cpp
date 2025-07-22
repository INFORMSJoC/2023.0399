#include "Instance.hpp"

int binCapacity;

Instance::Instance(): L0(-1), NT(0), size(0), L1(-1) {}

void Instance::add(int weight){
  this->weights.push_back(weight);
  this->size++;
  this->NT++;
}

void Instance::sort(){
  std::sort(all(this->weights));
}

Instance* Instance::ReadInstance(string pathInstance){
  Instance* instance;
  int size, numItems;
  ifstream instanceStream;

  instanceStream.open(pathInstance);

  if (!instanceStream.good()){
    throw new invalid_argument("Not found instance file");
  }

  instanceStream >> numItems >> binCapacity;
  instance = new Instance();

  FOR(i, numItems){
    instanceStream >> size;
    instance->add(size);
  }

  instanceStream.close();
  instance->sort();
  vi& freq = instance->freq;
  vi& weights = instance->weights;

  ll sumSize = 0;
  for(int w : weights)
    sumSize += w;
  instance->L0 = 1.0L * sumSize / binCapacity;
  instance->L1 = sumSize / binCapacity;

  freq.resize(instance->size);
  freq[0] = 1;

  int id = 0, last = 0;

  for(int i = 1; i < instance->size; i++){
    if(weights[last] != weights[i]){
      id++;
      last = i;
    }
    freq[id]++;
  }

  instance->size = id + 1;
  freq.resize(id + 1);
  weights.resize(unique(all(weights)) - weights.begin());
  assert(instance->size == SZ(weights));

  return instance;
}