#include <iostream>

#include <TFile.h>
#include <TTree.h>

using namespace std;

int main(int argc, char** argv)
{
  if (argc!=2) return 1;

  TFile *file = new TFile(argv[1],"read");
  cout << file->GetName() << endl;
  TTree *tree = (TTree*)file->Get("CollectionTree");

  Float_t m_yy = 0.;
  tree->SetBranchAddress("HGamEventInfoAuxDyn.m_yy",&m_yy);
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("HGamEventInfoAuxDyn.m_yy",1);

  for (Long64_t nent=tree->GetEntries(), ent=0; ent<nent; ++ent) {
    tree->GetEntry(ent);
    if (m_yy!=0.) cout << m_yy << endl; // Problem: never prints!
  }

  return 0;
}
