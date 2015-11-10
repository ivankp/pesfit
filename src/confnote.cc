#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>

#include "catstr.hh"
#include "structmap.hh"
#include "seqmap.hh"
#include "root_safe_get.hh"
#include "val_err.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

structmap(val_err<double>,hist_t,
  (nominal)(scale_down)(scale_up)(res_down)(res_up));

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "usage: " << argv[0] << " in.root out.pdf" << endl;
    return 0;
  }

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

  TH1
    *h_nominal    = get<TH1>(fin,"nominal");
    // *h_scale_down = get<TH1>(fin,"scale_down"),
    // *h_scale_up   = get<TH1>(fin,"scale_up"),
    // *h_res_down   = get<TH1>(fin,"res_down"),
    // *h_res_up     = get<TH1>(fin,"res_up");

  TGraph
    *f_nominal    = get<TGraph>(fin,"nominal_fit");
    // *f_scale_down = get<TGraph>(fin,"scale_down_fit"),
    // *f_scale_up   = get<TGraph>(fin,"scale_up_fit"),
    // *f_res_down   = get<TGraph>(fin,"res_down_fit"),
    // *f_res_up     = get<TGraph>(fin,"res_up_fit");

  // DRAW ************************************************

  TCanvas canv;
  canv.SetMargin(0.07,0.04,0.1,0.02);
  // canv.SetLogy();
  canv.SetTicks();
  canv.SaveAs(cat(argv[2],'[').c_str());

  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(18);
  lbl.SetNDC();

  h_nominal->SetLineWidth(1);
  h_nominal->SetLineColor(1);
  h_nominal->SetMarkerColor(1);
  h_nominal->GetYaxis()->SetTitleOffset(0.7);
  h_nominal->SetMinimum(0);
  
  f_nominal->SetLineColorAlpha(2,0.65);
  f_nominal->SetMarkerColor(2);
  f_nominal->SetFillColor(0);

  h_nominal->SetTitle("");
  h_nominal->Draw("E1");
  f_nominal->Draw("same");
  
  const double lxmin = 0.67;
  
  lbl.DrawLatex(lxmin,0.9,"ATLAS")->SetTextFont(73);
  lbl.DrawLatex(lxmin+0.095,0.9,"Internal");
  lbl.DrawLatex(lxmin,0.85,"#it{#sqrt{s}} = 8 TeV");
  lbl.DrawLatex(lxmin,0.80,"#it{H#rightarrow#gamma#gamma}, #it{m_{H}} = 125 GeV");
  
  TLegend leg(lxmin,0.73,lxmin+0.21,0.785);
  leg.AddEntry(h_nominal,"MC","le");
  leg.AddEntry(f_nominal,"Model","l");
  leg.SetBorderSize(0);
  leg.SetNColumns(2);
  leg.SetFillStyle(0);
  leg.Draw();
  
  lbl.DrawLatex(lxmin,0.65,"Fit parameters:");
  lbl.DrawLatex(lxmin,0.60,cat(
    "#mu_{CB} = ",
    make_pair(stats["mean_offset_bin0"].nominal+125.," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,0.55,cat(
    "#sigma_{CB} = ",
    make_pair(stats["sigma_offset_bin0"].nominal," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,0.50,cat(
    "#alpha_{CB} = ",
    make_pair(stats["crys_alpha_bin0"].nominal," #pm ")
  ).c_str());
  lbl.DrawLatex(lxmin,0.45,cat(
    "#it{n}_{CB} = ",
    setprecision(3),
    stats["crys_norm_bin0"].nominal.val,
    " (fixed)"
  ).c_str());
  lbl.DrawLatex(lxmin,0.40,cat(
    "#mu_{GA} = ",
    make_pair(stats["gaus_mean_offset_bin0"].nominal+125.," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,0.35,cat(
    "#kappa_{GA} = ",
    make_pair(stats["gaus_kappa_bin0"].nominal," #pm ")
  ).c_str());
  lbl.DrawLatex(lxmin,0.30,cat(
    "#phi_{CB} = ",
    make_pair(stats["fcb_bin0"].nominal," #pm ")
  ).c_str());
  lbl.DrawLatex(lxmin,0.24,cat(
    "FWHM_{ fit} = ",
    setprecision(3), stats["FWHM"].nominal.val,
    " [GeV]"
  ).c_str());

  canv.SaveAs(argv[2]);
  
  canv.SaveAs(cat(argv[2],']').c_str());

  delete fin;
  return 0;
}
