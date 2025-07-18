#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using ld = long double;

int main(int argc, char* argv[]){
  ifstream solIn, instIn;
  solIn.open(argv[1]);
  instIn.open(argv[2]);
  int WSol, capInst, nSol, nInst, L1, L2, Z;
  ld L0, L2frat, rootTime, finalTime;
  string solStatus, lx;

  instIn >> nInst >> capInst;
  
  map<int, int> sizesInst, sizesSol;
  for(int i = 0; i < nInst; i++){
    int sz;
    instIn >> sz;
    sizesInst[sz]++;
  }

  solIn >> lx >> nSol;
  solIn >> lx >> WSol;
  solIn >> lx >> L0;
  solIn >> lx >> L1;
  solIn >> lx >> L2frat;
  solIn >> lx >> L2;
  solIn >> lx >> Z;
  solIn >> lx >> solStatus;
  solIn >> lx >> rootTime;
  solIn >> lx >> finalTime;

  if(solIn.eof()){
    cout << "Problema na leitura da instancia" << endl;
    exit(1);
  }

  if(L1 == 0)
    return 0;

  if (L0 < L1 || L2 < Z || L1 < Z || rootTime > finalTime){
    cout << "Incoerencia na solucao" << endl;
    exit(1);
  }
  
  for(int i = 0; i < Z; i++){
    int capUsed, sizeBin, cap = 0, size;
    solIn >> lx >> lx >> capUsed >> lx;
    solIn >> lx >> sizeBin;
    solIn >> lx;
    for(int j = 0; j < sizeBin; j++){
      solIn >> lx >> size;
      if(size < 0){
        cout << "Item com cap negativa na bin: " << i << endl;
        exit(1);
      }
      cap += size;
      sizesSol[size]++;
    }
    if(cap < WSol || cap != capUsed){
      cout << "Cap incorreta na bin: " << i << endl;
      exit(1);
    }
  }
  {
    string bin;
    solIn >> bin >> lx >> lx >> lx >> lx;
    if (bin == "Bin"){
      cout << "Solucao mentiu ao dizer o numero de bins" << endl;
      exit(1);
    }
  }

  if (sizesInst.size() > sizesSol.size()){
    cout << "Ta faltando " << sizesInst.size()- sizesSol.size() << " item(s) na solucao" << endl;
    exit(1);
  }

  if (sizesInst.size() < sizesSol.size()){
    cout << "Ta sobrando " << sizesSol.size()- sizesInst.size() << " item(s) na solucao" << endl;
    exit(1);
  }

  for(auto pr : sizesInst){
    int size = pr.first, freq = pr.second;
    if(freq < sizesSol[size]){
      cout << "Frequencia do item " << size << " deveria ser " << freq << " mas eh " <<  sizesSol[size] << endl;
      exit(1);
    }
  }

  cout << "Solution ok" << endl;
  solIn.close();
  instIn.close();
  return 0;
}
