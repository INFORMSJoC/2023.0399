#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using ld = long double;

int main(int argc, char* argv[]){
  ifstream solIn, instIn;
  solIn.open(argv[1]);
  instIn.open(argv[2]);
  int MSol, MInst, nSol, nInst, B1, Z;
  ld B0, finalTime;
  string solStatus, lx;

  instIn >> nInst >> MInst;
  
  map<int, int> sizesInst, sizesSol;
  for(int i = 0; i < nInst; i++){
    int sz;
    instIn >> sz;
    sizesInst[sz]++;
  }

  solIn >> lx >> nSol;
  solIn >> lx >> MSol;
  solIn >> lx >> B0;
  solIn >> lx >> B1;
  solIn >> lx >> Z;
  solIn >> lx >> solStatus;
  solIn >> lx >> finalTime;

  if(solIn.eof()){
    cout << "Problema na leitura da instancia" << endl;
    exit(1);
  }

  if (B0 > B1 || B1 > Z || MSol > MInst){
    cout << "Incoerencia na solucao" << endl;
    exit(1);
  }
  
  for(int i = 0; i < MSol; i++){
    int capFree, sizeBin, cap = 0, size;
    string binId;
    solIn >> lx >> binId >> capFree >> lx;
    
    binId.pop_back();
    assert(stoi(binId) == i + 1);

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
    if(cap > Z || Z - cap != capFree){
      cout << cap << " " << Z << " " << capFree << endl;
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