#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>

#include "regex.hh"
#include "catstr.hh"
#include "root_safe_get.hh"
#include "TGraph_fcns.hh"
#include "binned.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout << "usage: " << argv[0] << " in.root out.pdf" << endl;
    return 0;
  }

  binned<TH1*> hmap({0,5,10,15,20,25,30});

  for (size_t i=1, n=hmap.nbins(); i<=n; ++i)
    hmap.at(i) = new TH1D(
      cat("m_yy_nvert[",hmap.left_edge(i),',',hmap.right_edge(i),')').c_str(),
      ";m_{#gamma#gamma} [GeV];",100,105,140);

  for (TH1 *h : hmap) cout << h->GetName() << endl;

  TFile *fin = new TFile(argv[1],"read");
  if (fin->IsZombie()) return 1;

/*
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

  TH1    *h_nominal = get<TH1>   (fin,"nominal");
  TGraph *f_nominal = get<TGraph>(fin,"nominal_fit");

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

  double lxmin = 0.67;
  double ly = 0.9;

  lbl.DrawLatex(lxmin,ly,"ATLAS")->SetTextFont(73);
  lbl.DrawLatex(lxmin+0.095,ly,"Internal");
  lbl.DrawLatex(lxmin,ly-=0.06,"#it{#sqrt{s}} = 8 TeV");
  lbl.DrawLatex(lxmin,ly-=0.06,"#it{H#rightarrow#gamma#gamma}, #it{m_{H}} = 125 GeV");

  ly-=0.075;
  TLegend leg(lxmin,ly,lxmin+0.21,ly+0.055);
  leg.AddEntry(h_nominal,"MC","le");
  leg.AddEntry(f_nominal,"Model","l");
  leg.SetBorderSize(0);
  leg.SetNColumns(2);
  leg.SetFillStyle(0);
  leg.Draw();

  lxmin = 0.12;
  ly = 0.9;

  lbl.DrawLatex(lxmin,ly,"Fit parameters:");
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#mu_{CB} = ",
    make_pair(stats["mean_offset_bin0"].nominal+125.," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#sigma_{CB} = ",
    make_pair(stats["sigma_offset_bin0"].nominal," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#alpha_{CB} = ",
    make_pair(stats["crys_alpha_bin0"].nominal," #pm ")
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#it{n}_{CB} = ",
    fixed, setprecision(1),
    stats["crys_norm_bin0"].nominal.val,
    " (fixed)"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#mu_{GA} = ",
    make_pair(stats["gaus_mean_offset_bin0"].nominal+125.," #pm "),
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#kappa_{GA} = ",
    make_pair(stats["gaus_kappa_bin0"].nominal," #pm ")
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "#phi_{CB} = ",
    make_pair(stats["fcb_bin0"].nominal," #pm ")
  ).c_str());

  lbl.DrawLatex(lxmin,ly-=0.10,cat(
    "FWHM_{ fit} = ",
    setprecision(3), stats["FWHM"].nominal.val,
    " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "68% : (",
    fixed, setprecision(1),
    ltailx(f_nominal,0.16), ',', rtailx(f_nominal,0.16),
    ") [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,cat(
    "90% : (",
    fixed, setprecision(1),
    ltailx(f_nominal,0.05), ',', rtailx(f_nominal,0.05),
    ") [GeV]"
  ).c_str());

  canv.SaveAs(argv[2]);

  canv.SaveAs(cat(argv[2],']').c_str());
  */

  delete fin;
  return 0;
}
