#include <iostream>
#include <string>
#include <utility>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>

#include "root_safe_get.hh"
#include "assert_ext.hh"
#include "catstr.hh"
#include "branches.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "usage: " << argv[0] << " mc.root output.pdf" << endl;
    return 1;
  }
  assert_ext(argv[1],{"root"});
  assert_ext(argv[2],{"pdf"});

  TFile *file = new TFile(argv[1],"read");
  TTree *tree = get<TTree>(file,"mc");

  Double_t myy, weight;
  branches(tree,
    "m_yy", &myy,
    "weight", &weight
  );

  TH1D hist("m_yy","",1000,105,160);

  for (Long64_t i=0, n=tree->GetEntries(); i<n; ++i) {
    tree->GetEntry(i);
    hist.Fill(myy,weight);
  }
  
  TCanvas canv;
  //canv.SetLogy();
  hist.Fit("gaus","","",120,130);
  hist.SetAxisRange(124,126);
  hist.Draw();
  canv.SaveAs(argv[2]);

  delete file;
  return 0;
}
