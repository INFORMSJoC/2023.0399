#ifndef INSTANCE_H
#define INSTANCE_H

#include "Pattern.hpp"

class Instance{
  public:
    string pathInstance, pathOutSol, pathCheckSol, prefix;
    bool resultCheckSol = false;
    ld L0;
    int NT, L1;
    vi freq;
    map<int, ii> inv;
    Instance();
    static Instance* ReadInstance(string pathInstance);
};

#endif
