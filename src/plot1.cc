#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

int err_prec(double err, int n=2) noexcept {
  const double log_err = std::log10(err);
  if (log_err < -8) return 0;
  else if (log_err < 0) return n - log_err;
  else if (log_err < n-1) return n-1;
  else return 0;
}

struct fit {
  TFile *file;
  TTree *tree;
  TGraph *hist, *fcn;
  map<string,pair<Double_t,Double_t>> stats;

  fit(const char* fname)
  : file(new TFile(fname,"read"))
  {
    if (file->IsZombie()) throw runtime_error(
      string("Bad file name ")+fname
    );
    hist = dynamic_cast<TGraph*>(file->Get("hist"));
    fcn  = dynamic_cast<TGraph*>(file->Get("fcn"));
    tree = dynamic_cast<TTree*>(file->Get("stats"));

    for (auto* x : *tree->GetListOfBranches())
      tree->SetBranchAddress(x->GetName(), (void*)&stats[x->GetName()]);
    tree->GetEntry(0);
  }
  ~fit() {
    delete file;
  }

  string operator()(const string& name, int prec=-1, bool single=false) const {
    const auto& par = stats.at(name);
    stringstream ss;
    if (prec==-1) {
      prec = err_prec(par.second);
      ss << fixed << setprecision(prec) << par.first << " #pm "
         << setprecision(prec) << par.second;
    } else {
      if (single) {
        ss << fixed << setprecision(prec) << par.first;
      } else {
        ss << fixed << setprecision(prec)
           << '(' << par.first <<", "<< par.second << ')';
      }
    }
    return ss.str();
  }
};

int main(int argc, char** argv)
{
  if (argc!=4) {
    cout << "usage: " << argv[0]
         << " fit_8TeV.root fit_13TeV.root plot.pdf" << endl;
    return 1;
  }

  fit f8(argv[1]), f13(argv[2]);

  TCanvas canv;
  canv.SetMargin(0.06,0.01,0.1,0.02);
  // canv.SetLogy();
  canv.SetTicks();

  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(18);
  lbl.SetNDC();

  f8.hist->SetMarkerSize(1);
  f8.hist->SetMarkerColor(30);
  f8.hist->SetTitle(";#it{m}_{#gamma#gamma} [GeV]");
  f8.hist->Draw("AP");
  f8.hist->GetXaxis()->SetRangeUser(104.5,140.5);
  f8.hist->GetXaxis()->SetTitleOffset(1.2);
  canv.Update();
  f8.fcn ->SetLineColorAlpha(3,0.5);
  f8.fcn ->Draw("Lsame");

  f13.hist->SetMarkerSize(1);
  f13.hist->SetMarkerColor(46);
  f13.hist->Draw("Psame");
  f13.fcn ->SetLineColorAlpha(2,0.5);
  f13.fcn ->Draw("Lsame");

  double lxmin = 0.67;
  double ly = 0.9;

  lbl.DrawLatex(lxmin,ly,"ATLAS")->SetTextFont(73);
  lbl.DrawLatex(lxmin+0.095,ly,"Internal");
  lbl.DrawLatex(lxmin,ly-=0.06,
    "#it{H#rightarrow#gamma#gamma}, #it{m_{H}} = 125 GeV");

  ly-=0.05;
  TLegend leg(lxmin,ly-0.215,lxmin+0.25,ly);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(f8.hist,"8 TeV MC","lep");
  leg.AddEntry(f8.fcn, "8 TeV Model","l");
  leg.AddEntry(f13.hist,"13 TeV MC","lep");
  leg.AddEntry(f13.fcn, "13 TeV Model","l");
  leg.Draw();

  lxmin = 0.12;
  ly = 0.9;

  lbl.DrawLatex(lxmin,ly,"13 TeV model fit parameters:");
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#mu_{CB} = " + f13("crys_mean") + " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#sigma_{CB} = " + f13("sigma_offset") + " [GeV]"
  ).c_str());

  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#mu_{GA} = " + f13("gaus_mean") + " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#kappa_{GA} = " + f13("gaus_kappa") + " [GeV]"
  ).c_str());

  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#phi_{CB} = " + f13("fcb")
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#alpha_{CB} = " + f13("crys_alpha")
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#it{n}_{CB} = " + f13("crys_norm",1,true) + " (fixed)"
  ).c_str());

  ly-=0.02;
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "FWHM = " + f13("fwhm",2,true) + " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#sigma_{68} = " + f13("sigma_68",2) + " [GeV]"
  ).c_str());
  lbl.DrawLatex(lxmin,ly-=0.05,(
    "#sigma_{90} = " + f13("sigma_90",2) + " [GeV]"
  ).c_str());

  canv.SaveAs(argv[3]);

  return 0;
}
