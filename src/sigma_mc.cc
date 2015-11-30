#include <iostream>
#include <string>
#include <utility>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

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
    cout << "usage: " << argv[0] << " frac mc.root" << endl;
    return 1;
  }
  double frac = atof(argv[1]);
  double tail = (1.-frac)/2.;
  assert_ext(argv[2],{"root"});

  TFile *file = new TFile(argv[2],"read");
  TTree *tree = get<TTree>(file,"mc");

  Double_t myy, weight;
  branches(tree,
    "m_yy", &myy,
    "weight", &weight
  );

  TH1D hist("m_yy","",55000,105,160);

  for (Long64_t i=0, n=tree->GetEntries(); i<n; ++i) {
    tree->GetEntry(i);
    hist.Fill(myy,weight);
  }

  Double_t total = 0.;
  for (Int_t i=1, n=hist.GetSize()-1; i<n; ++i)
    total += hist[i];
  test(total)

  frac *= total;
  tail *= total;

  cout << "\n" << frac << " interval from tails" << endl;

  Double_t left = 0.;
  for (Int_t i=1, n=hist.GetSize()-1; i<n; ++i) {
    left += hist[i];
    if (left > tail) {
      left = hist.GetBinLowEdge(i);
      break;
    }
  }
  test(left)

  Double_t right = 0.;
  for (Int_t i=hist.GetSize()-2; i>0; --i) {
    right += hist[i];
    if (right > tail) {
      right = hist.GetBinLowEdge(i+1);
      break;
    }
  }
  test(right)
  cout << "width = " << right-left << endl;

  cout << "\n" << frac << " interval from center" << endl;

  const Int_t bin125 = hist.FindFixBin(125.);
  Int_t a=bin125, b=bin125;
  bool l = true;
  Double_t center = hist[bin125];
  while (center < frac) {
    center += ( (l = !l) ? hist[--a] : hist[++b] );
  }
  left  = hist.GetBinLowEdge(l ? a+1 : a);
  right = hist.GetBinLowEdge(l ? b+1 : b);
  test(left)
  test(right)
  cout << "width = " << right-left << endl;

  delete file;
  return 0;
}
