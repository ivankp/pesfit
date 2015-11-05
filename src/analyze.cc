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
#include <TH2.h>
#include <TGraph.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
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

int main(int argc, char** argv)
{
  TFile *fin = new TFile(argv[2],"read");
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

  TH1
    // *h_nominal    = get<TH1>(fin,"nominal"),
    *h_scale_down = get<TH1>(fin,"scale_down"),
    *h_scale_up   = get<TH1>(fin,"scale_up"),
    *h_res_down   = get<TH1>(fin,"res_down"),
    *h_res_up     = get<TH1>(fin,"res_up");

  TGraph
    // *f_nominal    = get<TGraph>(fin,"nominal_fit"),
    *f_scale_down = get<TGraph>(fin,"scale_down_fit"),
    *f_scale_up   = get<TGraph>(fin,"scale_up_fit"),
    *f_res_down   = get<TGraph>(fin,"res_down_fit"),
    *f_res_up     = get<TGraph>(fin,"res_up_fit");

  workspace ws(argv[1]);

  auto fix = [&ws](const char* name, double x){
    auto *var = ws->var(name);
    test(x)
    var->setRange(x,x);
    var->setVal(x);
  };

  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  canv.SetLogy();
  canv.SaveAs(cat(argv[3],'[').c_str());

  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(20);
  lbl.SetNDC();

  auto draw_stat = [&lbl](TH1* hist, double y, double x){
    auto lblp = lbl.DrawLatex(0.14,y,hist->GetName());
    lblp->SetTextColor(hist->GetLineColor());
    lblp->DrawLatex(0.31,y,cat(x).c_str());
  };

  // SCALE ************************************************

  h_scale_down->SetTitle("Unconstrained fit");
  h_scale_down->Draw();
  f_scale_down->Draw("same");
  h_scale_up->Draw("same");
  f_scale_up->Draw("same");

  draw_stat(h_scale_down,0.84,stats["mean_offset_bin0"].scale_down);
  draw_stat(h_scale_up  ,0.80,stats["mean_offset_bin0"].scale_up  );

  canv.SaveAs(argv[3]);

  h_scale_down->SetTitle("Mean offset set to window mean - 125");
  fix("mean_offset_bin0", stats["hist_window_mean"].scale_down - 125.);
  auto fit_res = ws.fit(h_scale_down);
  h_scale_down->Draw();
  fit_res.second->Draw("same");
  draw_stat(h_scale_down,0.84,
    static_cast<RooRealVar*>(
      fit_res.first->floatParsFinal().find("mean_offset_bin0")
    )->getVal());

  fix("mean_offset_bin0", stats["hist_window_mean"].scale_up - 125.);
  fit_res = ws.fit(h_scale_up);
  h_scale_up->Draw("same");
  fit_res.second->Draw("same");
  draw_stat(h_scale_up,0.80,
    static_cast<RooRealVar*>(
      fit_res.first->floatParsFinal().find("mean_offset_bin0")
    )->getVal());

  canv.SaveAs(argv[3]);

  // RESOLUTION *******************************************

  h_res_down->SetTitle("Unconstrained fit");
  h_res_down->Draw();
  f_res_down->Draw("same");
  h_res_up->Draw("same");
  f_res_up->Draw("same");

  draw_stat(h_res_down,0.84,stats["sigma_offset_bin0"].res_down);
  draw_stat(h_res_up  ,0.80,stats["sigma_offset_bin0"].res_up  );

  canv.SaveAs(argv[3]);

  h_res_down->SetTitle("Sigma offset set to HWHM");
  fix("sigma_offset_bin0", stats["FWHM"].res_down/2.);
  fit_res = ws.fit(h_res_down);
  h_res_down->Draw();
  fit_res.second->Draw("same");
  draw_stat(h_res_down,0.84,
    static_cast<RooRealVar*>(
      fit_res.first->floatParsFinal().find("sigma_offset_bin0")
    )->getVal());

  fix("sigma_offset_bin0", stats["FWHM"].res_up/2.);
  fit_res = ws.fit(h_res_up);
  h_res_up->Draw("same");
  fit_res.second->Draw("same");
  draw_stat(h_res_up,0.80,
    static_cast<RooRealVar*>(
      fit_res.first->floatParsFinal().find("sigma_offset_bin0")
    )->getVal());

  canv.SaveAs(argv[3]);

  canv.SaveAs(cat(argv[3],']').c_str());

  return 0;
}
