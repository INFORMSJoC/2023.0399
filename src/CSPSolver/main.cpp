#include "fileHandler.hpp"
#include "Solver.hpp"

int TIMELIMIT = 3600;

int main(int argc, char* argv[]){
  if(argc < 3){
    cerr << "Numero de argumentos invalido" << endl;
    exit(1);
  }

  string outdir = outdirName(argv[0]);
  Instance* instance = Instance::ReadInstance(argv[1]);
  instance->pathInstance = argv[1];
  instance->prefix = outdir + "/" + getFileName(instance->pathInstance);
  instance->pathOutSol = instance->prefix + ".sol";
  instance->pathCheckSol = argv[2];

  stringstream ss;

  for(int i = 3; i < argc; i++)
    ss << argv[i] << " ";

  unordered_set<string> valid_features = {"RF", "CRF", "Waste", "Cost", "Cut", "Best", "Ineq", "Splay", "Historic"};
  string tp, val;
  while(ss >> tp){
    if(tp[0] == '-'){
      tp.erase(remove(tp.begin(), tp.end(), '-'), tp.end());
      if(tp.empty() && (ss >> tp))
        throw invalid_argument("Invalid input arguments");
      if(tp == "TL"){
        if(ss >> val)
          TIMELIMIT = stoi(val);
        else
          throw invalid_argument("Invalid input arguments");
      }else{
        if(!valid_features.count(tp))
          throw invalid_argument("Invalid input arguments");
        features.insert(tp);
      }
    }
    else
      throw invalid_argument("Invalid input arguments");
  }
  if (SZ(features) == 0){
    features = valid_features;
    features.erase("Best");
  }
    
  for (auto feat : features)
    cout << feat << " ";
  cout << "\n";

  createDir(outdir);

  chrono::steady_clock::time_point begin, now;
  begin = chrono::steady_clock::now();

  Solver* solver =  new Solver(instance);
  solver->Solve();
  delete solver;

  FOR(i, SZ(Pattern::dels))
    assert(Pattern::dels[i]);
  delete instance;

  return 0;
}
