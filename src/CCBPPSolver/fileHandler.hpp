#include "bits/stdc++.h"
#include <sys/stat.h>

using namespace std;

string getFileName(string path){
  int begin = 0;
  size_t found = path.find("/");
  while(found != string::npos){
    begin = found + 1;
    found = path.find("/", begin);
  }
  return path.substr(begin);
}

string outdirName(string currentDir){
  string outdir = currentDir;

  while(outdir.back() != '/')
    outdir.pop_back();

  if(outdir == (string)"./")
    outdir.pop_back(), outdir.pop_back();

  outdir += "out";
  return outdir;
}

void createDir(string& dir){
  struct stat buffer;
  if (stat(dir.c_str(), &buffer) != 0){
    string command = "mkdir " + dir;
    assert(system(command.c_str()) != -1);
  }
}
