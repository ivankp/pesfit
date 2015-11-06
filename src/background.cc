#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooFitResult.h>
#include <RooCurve.h>

#include "catstr.hh"
#include "seqmap.hh"
#include "structmap.hh"
#include "root_safe_get.hh"
#include "workspace.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

structmap(double,hist_t,
  (nominal)(scale_down)(scale_up)(res_down)(res_up));

Double_t exp2(const Double_t* x, const Double_t* p) noexcept {
  const Double_t m = (*x-100.)/100.;
  return exp(p[0]*m - p[1]*m*m);
}

int main(int argc, char** argv)
{
  TFile *fin = new TFile(argv[1],"read");
  if (fin->IsZombie()) return 1;

  seqmap<hist_t> stats;
  {
    TTree *tree = get<TTree>(fin,"stats");
    auto *branches = tree->GetListOfBranches();
    stats.reserve(branches->GetEntries());
    for (auto x : *branches) {
      cout << "Branch: " << x->GetName() << endl;
      tree->SetBranchAddress(x->GetName(),&stats[x->GetName()]);
    }
    tree->GetEntry(0);
  }
  cout << endl;

  TH1 *h_nominal = get<TH1>(fin,"nominal");
  const double Nsig = 124.84;
  h_nominal->Scale(Nsig/h_nominal->Integral());

  Double_t params[2] = { -3.4472e+00, 4.6032e-01 };
  TF1 *fexp2 = new TF1("exp2",exp2,105,140,2);
  fexp2->SetParameters(params);
  
  TH1 *hexp2 = new TH1D("hexp2",";m_{#gamma#gamma} [GeV];N Events",
    h_nominal->GetNbinsX(),
    h_nominal->GetXaxis()->GetXmin(),
    h_nominal->GetXaxis()->GetXmax());
  hexp2->Eval(fexp2);
  hexp2->SetStats(0);
  hexp2->SetLineColor(43);
  hexp2->SetMarkerColor(43);
  hexp2->SetLineWidth(2);
  const double Nbg = 24740.;
  hexp2->Scale(Nbg/hexp2->Integral());
  hexp2->Add(h_nominal);
  
  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  
  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(17);
  lbl.SetNDC();

  auto draw_stat = [&lbl](double y, const string& str){
    lbl.DrawLatex(0.75,y,str.c_str());
  };
  
  //canv.SetLogy();
  hexp2->Draw();
  // h_nominal->Draw("same");
  // hexp2->SetMinimum(0.5);
  // hexp2->SetMaximum(500);
  
  draw_stat(0.85,cat("a = ",params[0]));
  draw_stat(0.80,cat("b = ",params[1]));
  draw_stat(0.75,cat("N_{bg} = ",Nbg));
  draw_stat(0.70,cat("N_{sig} = ",Nsig));
  
  canv.SaveAs("exp2.pdf");
  
  delete fexp2;
  delete hexp2;

  return 0;
}
