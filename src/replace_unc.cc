#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

#include "root_safe_get.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

template<typename T> using pp = std::pair<T,T>;
  
int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "usage: " << argv[0]
         << " weighted.root unweighted.root" << endl;
    return 0;
  }
  
  TFile *fw = new TFile(argv[1],"read");
  if (fw->IsZombie()) return 1;

  vector<pair<string,pp<pp<Double_t>>>> stats;
  
  TTree *tw = get<TTree>(fw,"stats");
  auto *branches = tw->GetListOfBranches();
  stats.reserve(branches->GetEntries());
  for (auto* x : *branches) {
    //cout << "Branch: " << x->GetName() << endl;
    stats.emplace_back(x->GetName(),pp<pp<Double_t>>());
    tw->SetBranchAddress(x->GetName(), (void*)&stats.back().second.first);
  }
  tw->GetEntry(0);
  
  cout << endl;
  for (auto& x : stats)
    cout <<setw(25)<< x.first <<" = "<< x.second.first.first <<"\t"<< x.second.first.second << endl;
  cout << endl;
  
  TFile *fu = new TFile(argv[2],"read");
  if (fu->IsZombie()) return 1;
  
  TTree *tu = get<TTree>(fu,"stats");
  for (auto& x : stats) {
    //cout << "Branch: " << x.first << endl;
    tu->SetBranchAddress(x.first.c_str(), (void*)&x.second.second);
  }
  tu->GetEntry(0);
  
  cout << endl;
  for (auto& x : stats)
    cout <<setw(25)<< x.first <<" = "<< x.second.second.first <<"\t"<< x.second.second.second << endl;
  cout << endl;
  
  for (auto& x : stats) {
    if (x.first!="sigma_68")
      x.second.first.second = x.second.second.second;
  }
  
  cout << endl;
  for (auto& x : stats)
    cout <<setw(25)<< x.first <<" = "<< x.second.first.first <<"\t"<< x.second.first.second << endl;
  cout << endl;


  string fname(argv[1]);
  TFile *fo = new TFile(
    (fname.substr(0,fname.rfind('.'))+"_fixed.root").c_str(),"recreate");
  if (fo->IsZombie()) return 1;
  
  fo->cd();
  get<TGraph>(fw,"hist")->Write("hist");
  get<TGraph>(fw,"fcn")->Write("fcn");
  
  TTree *to = new TTree(tw->GetName(),tu->GetTitle());
  for (auto& x : stats)
    to->Branch(x.first.c_str(), (void*)&x.second.first,
               (x.first+"[2]/D").c_str());
  to->Fill();
  
  fo->Write(0,TObject::kOverwrite);
  
  delete fo;
  delete fu;
  delete fw;

  return 0;
}
