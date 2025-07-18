#ifndef INSTANCE_H
#define INSTANCE_H

#include "Pattern.hpp"

class Instance{
  private:
    void add(int weight);
    void sort();

  public:
    string pathInstance, pathOutSol, pathCheckSol, prefix;
    bool resultCheckSol = false;
    ld L0, B0;
    int NT, size, L1, B1;
    ll sumSize = 0;
    vi weights, freq;
    Instance();
    static Instance* ReadInstance(string pathInstance);
    void computeLBs();
};

#endif
