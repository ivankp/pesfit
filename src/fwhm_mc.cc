#include <iostream>
#include <string>
#include <utility>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
// #include <TCanvas.h>

#include "root_safe_get.hh"
#include "assert_ext.hh"
#include "catstr.hh"
#include "branches.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout << "usage: " << argv[0] << " mc.root" << endl;
    return 1;
  }
  assert_ext(argv[1],{"root"});

  TFile *file = new TFile(argv[1],"read");
  TTree *tree = get<TTree>(file,"mc");

  Double_t myy, weight;
  branches(tree,
    "m_yy", &myy,
    "weight", &weight
  );

  TH1D hist("m_yy","",500,105,160);
  const double bw = hist.GetBinWidth(1);

  for (Long64_t i=0, n=tree->GetEntries(); i<n; ++i) {
    tree->GetEntry(i);
    hist.Fill(myy,weight);
  }

  // TCanvas canv;
  // canv.SetLogy();
  // hist.Draw();
  // canv.SaveAs("mc_hist.pdf");

  int k = 1;
  for (int i=(120-105)/bw, b=(130-105)/bw; i<b; ++i)
    if (hist[i] > hist[k]) k = i;

  int a = k;
  for (int i=k; i>0; --i)
    if (hist[i] < hist[k]/2) {
      a = i + 1;
      break;
    }

  int b = k;
  for (int i=k, end=hist.GetNbinsX(); i<end; ++i)
    if (hist[i] < hist[k]/2) {
      b = i;
      break;
    }

  double left  = hist.GetBinLowEdge(a);
  double right = hist.GetBinLowEdge(b);

  test(left)
  test(right)
  cout << "width = " << right-left << endl;

  delete file;
  return 0;
}
