#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using ld = long double;
using ii = pair<int, int>;

struct Item{
  int w, c;
  bool operator < (Item it2) const{
    return ii(w, c) < ii(it2.w, it2.c);
  }
};

int main(int argc, char* argv[]){
  ifstream solIn, instIn;
  solIn.open(argv[1]);
  instIn.open(argv[2]);
  int WSol, CInst, QInst, capInst, nSol, nInst, L1, L2, Z;
  ld L0, L2frat, rootTime, finalTime;
  string solStatus, lx;

  instIn >> nInst >> QInst >> capInst >> CInst;
  
  map<Item, int> sizesInst, sizesSol;
  for(int i = 0; i < nInst; i++){
    int sz, cl;
    instIn >> sz >> cl;
    sizesInst[{sz, cl}]++;
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

  if (L0 > L1 || L2 > Z || L1 > Z || rootTime > finalTime){
    cout << "Incoerencia na solucao" << endl;
    exit(1);
  }
  
  for(int i = 0; i < Z; i++){
    int capSol, sizeBin, capIn = 0, size, cl;
    solIn >> lx >> lx >> capSol >> lx;
    solIn >> lx >> sizeBin;
    solIn >> lx;
    map<int, int> colors;
    for(int j = 0; j < sizeBin; j++){
      solIn >> lx >> size >> lx >> cl;
      if(size < 0){
        cout << "Item com cap negativa na bin: " << i << endl;
        exit(1);
      }
      colors[cl]++;
      capIn += j + 1 == sizeBin ? 1 : size;
      sizesSol[{size, cl}]++;
    }
    if(capIn > WSol){
      cout << "Cap incorreta na bin: " << i << endl;
      exit(1);
    }
    if(colors.size() > CInst){
      cout << "Number incorrect of colors in bin: " << i << endl;
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
    Item item = pr.first;
    int freq = pr.second;
    if(freq < sizesSol[item]){
      cout << "Frequencia do item " << item.c << " deveria ser " << freq << " mas eh " <<  sizesSol[item] << endl;
      exit(1);
    }
  }

  cout << "Solution ok" << endl;
  solIn.close();
  instIn.close();
  return 0;
}