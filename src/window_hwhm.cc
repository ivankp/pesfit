#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooFitResult.h>
//#include <RooCurve.h>

#include "catstr.hh"
#include "seqmap.hh"
#include "structmap.hh"
#include "root_safe_get.hh"
#include "val_err.hh"
#include "workspace.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

structmap(val_err<double>,hist_t,
  (nominal)(scale_down)(scale_up)(res_down)(res_up));

Double_t exp2(const Double_t* x, const Double_t* p) noexcept {
  const Double_t m = (*x-100.)/100.;
  return exp(p[0]*m - p[1]*m*m);
}

int main(int argc, char** argv)
{
  string ifname, ofname, wfname, cfname;
  bool logy, bg, dopull;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->required(),
       "*input root file")
      ("output,o", po::value(&ofname)->required(),
       "*output pdf file name")
      ("workspace,w", po::value(&wfname)->default_value("data/ws.root"),
       "ROOT file with RooWorkspace for CB fits")
      ("config,c", po::value(&cfname),
       "configuration file name")

      ("logy,l", po::bool_switch(&logy),
       "logarithmic Y axis")
      ("bg,b", po::bool_switch(&bg),
       "add exp2 background to signal")
      ("pull,p", po::bool_switch(&dopull),
       "calculate fit pulls")
    ;

    po::positional_options_description pos;
    pos.add("input",1);
    pos.add("output",1);
    pos.add("workspace",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("config")) {
      po::store( po::parse_config_file<char>(
        vm["config"].as<string>().c_str(), desc), vm);
    }
    po::notify(vm);
  } catch (std::exception& e) {
    cerr << "\033[31mArgs: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  TFile *fin = new TFile(ifname.c_str(),"read");
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
    *h_nominal    = get<TH1>(fin,"nominal"),
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

  workspace ws(wfname,bg);

  if (bg) { // add background
    Double_t params[2] = { -3.4472e+00, 4.6032e-01 };
    TF1 fexp2("exp2",exp2,105,140,2);
    fexp2.SetParameters(params);

    TH1 *h_bg = new TH1D("hexp2",";m_{#gamma#gamma} [GeV];N Events",
      h_nominal->GetNbinsX(),
      h_nominal->GetXaxis()->GetXmin(),
      h_nominal->GetXaxis()->GetXmax());
    h_bg->Eval(&fexp2);
    h_bg->Scale(24740./h_bg->Integral());

    auto& nsig_mc = stats["nsig_mc"];
    const double scale_to_nsig_nom
      = (nsig_mc.nominal = 124.84)/h_nominal->Integral();

    h_scale_down->Scale(scale_to_nsig_nom);
    h_scale_up  ->Scale(scale_to_nsig_nom);
    h_res_down  ->Scale(scale_to_nsig_nom);
    h_res_up    ->Scale(scale_to_nsig_nom);

    nsig_mc.scale_down = h_scale_down->Integral();
    nsig_mc.scale_up   = h_scale_up  ->Integral();
    nsig_mc.res_down   = h_res_down  ->Integral();
    nsig_mc.res_up     = h_res_up    ->Integral();

    h_scale_down->Add(h_bg);
    h_scale_up  ->Add(h_bg);
    h_res_down  ->Add(h_bg);
    h_res_up    ->Add(h_bg);
  }

  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  if (logy) canv.SetLogy();
  canv.SaveAs((ofname+'[').c_str());

  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(20);
  lbl.SetNDC();

  auto draw_stat = [&lbl](const TH1* hist, double y, double x){
    auto lblp = lbl.DrawLatex(0.14,y,hist->GetName());
    lblp->SetTextColor(hist->GetLineColor());
    lblp->DrawLatex(0.31,y,cat(x).c_str());
  };

  // SCALE ************************************************

  h_scale_down->SetTitle("Unconstrained fit");
  h_scale_down->Draw();
  h_scale_up->Draw("same");

  if (bg) {
    lbl.DrawLatex(0.62,0.85,"NSig_fit");
    lbl.DrawLatex(0.73,0.85,"pull");
    lbl.DrawLatex(0.84,0.85,"rel_diff");

    const auto& nsig_mc_stat = stats["nsig_mc"];

    auto fit = ws.fit(h_scale_down);
    fit.second->Draw("same");
    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.scale_down) / nsig_mc_stat.scale_down;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.80,"scale_down");
    lblp->SetTextColor(h_scale_down->GetLineColor());
    lblp->DrawLatex(0.62,0.80,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.80,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.80,Form("%.5f",rel_diff));

    fit = ws.fit(h_scale_up);
    fit.second->Draw("same");
    nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    rel_diff = (nsig - nsig_mc_stat.scale_up) / nsig_mc_stat.scale_up;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    lblp = lbl.DrawLatex(0.45,0.75,"scale_up");
    lblp->SetTextColor(h_scale_up->GetLineColor());
    lblp->DrawLatex(0.62,0.75,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.75,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.75,Form("%.5f",rel_diff));
  } else {
    f_scale_down->Draw("same");
    f_scale_up->Draw("same");
    draw_stat(h_scale_down,0.84,stats["mean_offset_bin0"].scale_down);
    draw_stat(h_scale_up  ,0.80,stats["mean_offset_bin0"].scale_up  );
  }

  canv.SaveAs(ofname.c_str());

  if (bg) h_scale_down->SetTitle("Should be exactly the same as previous");
  else h_scale_down->SetTitle("Mean offset set to window mean - 125");
  ws.fixVal("mean_offset_bin0", stats["hist_window_mean"].scale_down - 125.);
  auto fit = ws.fit(h_scale_down);

  test(fit.first->GetName())
  auto& pars = fit.first->floatParsFinal();
  for (int i=0, n=pars.getSize(); i<n; ++i)
    test(pars[i].GetName())

  h_scale_down->Draw();
  fit.second->Draw("same");

  if (bg) {
    lbl.DrawLatex(0.62,0.85,"NSig_fit");
    lbl.DrawLatex(0.73,0.85,"pull");
    lbl.DrawLatex(0.84,0.85,"rel_diff");

    const auto& nsig_mc_stat = stats["nsig_mc"];

    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.scale_down) / nsig_mc_stat.scale_down;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.80,"scale_down");
    lblp->SetTextColor(h_scale_down->GetLineColor());
    lblp->DrawLatex(0.62,0.80,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.80,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.80,Form("%.5f",rel_diff));
  } else draw_stat(h_scale_down,0.84,
    static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("mean_offset_bin0")
    )->getVal());

  ws.fixVal("mean_offset_bin0", stats["hist_window_mean"].scale_up - 125.);
  fit = ws.fit(h_scale_up);
  h_scale_up->Draw("same");
  fit.second->Draw("same");
  if (bg) {
    const auto& nsig_mc_stat = stats["nsig_mc"];

    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.scale_up) / nsig_mc_stat.scale_up;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.75,"scale_up");
    lblp->SetTextColor(h_scale_up->GetLineColor());
    lblp->DrawLatex(0.62,0.75,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.75,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.75,Form("%.5f",rel_diff));
  } else draw_stat(h_scale_up,0.80,
    static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("mean_offset_bin0")
    )->getVal());

  canv.SaveAs(ofname.c_str());

  // RESOLUTION *******************************************

  h_res_down->SetTitle("Unconstrained fit");
  h_res_down->Draw();
  h_res_up->Draw("same");

  if (bg) {
    lbl.DrawLatex(0.62,0.85,"NSig_fit");
    lbl.DrawLatex(0.73,0.85,"pull");
    lbl.DrawLatex(0.84,0.85,"rel_diff");

    const auto& nsig_mc_stat = stats["nsig_mc"];

    auto fit = ws.fit(h_res_down);
    fit.second->Draw("same");
    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.res_down) / nsig_mc_stat.res_down;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.80,"res_down");
    lblp->SetTextColor(h_res_down->GetLineColor());
    lblp->DrawLatex(0.62,0.80,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.80,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.80,Form("%.5f",rel_diff));

    fit = ws.fit(h_res_up);
    fit.second->Draw("same");
    nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    rel_diff = (nsig - nsig_mc_stat.res_up) / nsig_mc_stat.res_up;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    lblp = lbl.DrawLatex(0.45,0.75,"res_up");
    lblp->SetTextColor(h_res_up->GetLineColor());
    lblp->DrawLatex(0.62,0.75,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.75,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.75,Form("%.5f",rel_diff));
  } else {
    f_res_down->Draw("same");
    f_res_up->Draw("same");
    draw_stat(h_res_down,0.84,stats["sigma_offset_bin0"].res_down);
    draw_stat(h_res_up  ,0.80,stats["sigma_offset_bin0"].res_up  );
  }

  canv.SaveAs(ofname.c_str());

  if (bg) h_res_down->SetTitle("Should be exactly the same as previous");
  else h_res_down->SetTitle("Sigma offset set to HWHM");
  ws.fixVal("sigma_offset_bin0", stats["FWHM"].res_down/2.);
  fit = ws.fit(h_res_down);
  h_res_down->Draw();
  fit.second->Draw("same");
  if (bg) {
    lbl.DrawLatex(0.62,0.85,"NSig_fit");
    lbl.DrawLatex(0.73,0.85,"pull");
    lbl.DrawLatex(0.84,0.85,"rel_diff");

    const auto& nsig_mc_stat = stats["nsig_mc"];

    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.res_down) / nsig_mc_stat.res_down;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.80,"res_down");
    lblp->SetTextColor(h_res_down->GetLineColor());
    lblp->DrawLatex(0.62,0.80,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.80,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.80,Form("%.5f",rel_diff));
  } else draw_stat(h_res_down,0.84,
    static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("sigma_offset_bin0")
    )->getVal());

  ws.fixVal("sigma_offset_bin0", stats["FWHM"].res_up/2.);
  fit = ws.fit(h_res_up);
  h_res_up->Draw("same");
  fit.second->Draw("same");
  if (bg) {
    const auto& nsig_mc_stat = stats["nsig_mc"];

    Double_t nsig = static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("NSig_bin0")
    )->getVal();
    Double_t rel_diff = (nsig - nsig_mc_stat.res_up) / nsig_mc_stat.res_up;
    Double_t pull = NAN;
    if (dopull) {
      pull = static_cast<RooRealVar*>(
        fit.first->floatParsFinal().find("Uncert_EnRes_EnRes")
      )->getVal();
    }

    auto lblp = lbl.DrawLatex(0.45,0.75,"res_up");
    lblp->SetTextColor(h_res_up->GetLineColor());
    lblp->DrawLatex(0.62,0.75,Form("%.2f",nsig));
    lblp->DrawLatex(0.73,0.75,Form("%.3f",pull));
    lblp->DrawLatex(0.84,0.75,Form("%.5f",rel_diff));
  } else draw_stat(h_res_up,0.80,
    static_cast<RooRealVar*>(
      fit.first->floatParsFinal().find("sigma_offset_bin0")
    )->getVal());

  canv.SaveAs(ofname.c_str());

  canv.SaveAs((ofname+']').c_str());

  return 0;
}
