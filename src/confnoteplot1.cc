#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
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
#include "TGraph_fcns.hh"
#include "golden_min.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

structmap(val_err<double>,hist_t,
  (nominal)(scale_down)(scale_up)(res_down)(res_up));
  
int main(int argc, char** argv)
{
  if (argc!=3 && argc!=4) {
    cout << "usage: " << argv[0] << " in.root out.pdf [minsig]" << endl;
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
      //cout << "Branch: " << x->GetName() << endl;
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
  
  const Double_t integral = integrate(f_nominal);
  
  if (argc==4 && string("minsig")==argv[3]) {
    golden_min gm;
    
    for (auto frac : {0.68,0.90}) {
      test(frac)
      auto x1 = gm( [f_nominal,frac,integral](double x1) {
        return intervalx2(f_nominal, frac, x1, integral) - x1;
      }, firstx(f_nominal), rtailx(f_nominal,frac,integral) ).first;
      auto x2 = intervalx2(f_nominal, frac, x1, integral);
      cout << setprecision(6);
      cout << "x1 = " << x1 << endl;
      cout << "x2 = " << x2 << endl;
      lbl.DrawLatex(lxmin,ly-=0.05,cat(
        frac*100,"% : (",
        fixed, setprecision(2),
        x1, ',', x2,
        ") [GeV]"
      ).c_str());
    }
  } else {
    lbl.DrawLatex(lxmin,ly-=0.05,cat(
      "68% : (",
      fixed, setprecision(2),
      ltailx(f_nominal,0.16,integral), ',', rtailx(f_nominal,0.16,integral),
      ") [GeV]"
    ).c_str());
    lbl.DrawLatex(lxmin,ly-=0.05,cat(
      "90% : (",
      fixed, setprecision(2),
      ltailx(f_nominal,0.05,integral), ',', rtailx(f_nominal,0.05,integral),
      ") [GeV]"
    ).c_str());
  }

  canv.SaveAs(argv[2]);

  delete fin;
  return 0;
}
