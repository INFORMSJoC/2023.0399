#include "Instance.hpp"

int binCapacity, numMachines;

Instance::Instance(): L0(-1), NT(0), size(0), L1(-1) {}

void Instance::add(int weight){
  this->weights.push_back(weight);
  this->size++;
  this->NT++;
}

void Instance::sort(){
  std::sort(all(this->weights));
}

void Instance::computeLBs(){
  L0 = 1.0L * sumSize / binCapacity;
  L1 = (sumSize / binCapacity) + (sumSize % binCapacity != 0);
}

Instance* Instance::ReadInstance(string pathInstance){
  Instance* instance;
  int size, numItems;
  ifstream instanceStream;

  instanceStream.open(pathInstance);

  if (!instanceStream.good()){
    throw new invalid_argument("Not found instance file");
  }

  instanceStream >> numItems >> numMachines;
  instance = new Instance();

  FOR(i, numItems){
    instanceStream >> size;
    instance->add(size);
  }

  instanceStream.close();
  instance->sort();
  vi& freq = instance->freq;
  vi& weights = instance->weights;

  ll& sumSize = instance->sumSize;
  for(int w : weights)
    sumSize += w;

  instance->B0 = 1.0L * sumSize / numMachines;
  instance->B1 = (sumSize / numMachines) + (sumSize % numMachines != 0);
  assert(weights.front() <= weights.back());
  int B2 = weights.back();
  instance->B0 = max(B2, instance->B0);
  instance->B1 = max(B2, instance->B1);
  /*if(numItems > numMachines){
    int B3 = weights[numItems - numMachines] + weights[numItems - numMachines - 1];
    instance->B0 = max(B3, instance->B0);
    instance->B1 = max(B3, instance->B1);
  }*/

  binCapacity = weights.back();
  instance->computeLBs();
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