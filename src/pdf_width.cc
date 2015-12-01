#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <cmath>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLine.h>
#include <TLatex.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
//#include <RooCategory.h>

#include "root_safe_get.hh"
#include "assert_ext.hh"
#include "golden_min.hh"
#include "catstr.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;
  
struct workspace {
  TFile *file;
  RooWorkspace *ws;
  RooSimultaneous *sim_pdf;
  //RooCategory *cat;

public:
  workspace(const char* fname)
  : file(new TFile(fname,"read")),
    ws(get<RooWorkspace>(file,"mcfit")),
    sim_pdf(static_cast<RooSimultaneous*>(
      ws->obj("mc_sim_pdf_bin0")))
    //cat(static_cast<RooCategory*>(ws->obj("mc_sample")))
  { }
  ~workspace() {
    //delete cat;
    delete sim_pdf;
    delete file;
  }

  inline RooWorkspace* operator->() noexcept { return ws; }
  
  void set_var(const string& name, double x) {
    auto* var = static_cast<RooRealVar*>(ws->var(name.c_str()));
    if (!var) throw runtime_error(
      "No variable "+name+" in workspace "+ws->GetName()
    );
    var->setRange(x,x);
    var->setVal(x);
    cout << name << " = " << var->getVal() << endl;
  }
  
  TH1* m_yy_hist(int nbins) const {
    return sim_pdf->createHistogram("m_yy",nbins);
  }
};

template<typename T>
pair<double,double> get_fwhm(const T& hist) {
  const double bw = hist.GetBinWidth(1);

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
    
  return {hist.GetBinLowEdge(a), hist.GetBinLowEdge(b)};
}

template<typename T>
pair<double,double> get_sigma_tails(const T& hist, double frac) {
  Double_t total = 0.;
  for (Int_t i=1, n=hist.GetSize()-1; i<n; ++i)
    total += hist[i];

  Double_t tail = total*(1.-frac)/2.;
  frac *= total;

  Double_t left = 0.;
  for (Int_t i=1, n=hist.GetSize()-1; i<n; ++i) {
    left += hist[i];
    if (left > tail) {
      left = hist.GetBinLowEdge(i);
      break;
    }
  }

  Double_t right = 0.;
  for (Int_t i=hist.GetSize()-2; i>0; --i) {
    right += hist[i];
    if (right > tail) {
      right = hist.GetBinLowEdge(i+1);
      break;
    }
  }
  
  return {left,right};
}

template<typename T>
pair<double,double> get_sigma_mode(const T& hist, double frac) {
  Double_t total = 0.;
  Int_t mode = 1;
  for (Int_t i=1, n=hist.GetSize()-1; i<n; ++i) {
    total += hist[i];
    if (hist[mode] < hist[i]) mode = i;
  }

  frac *= total;

  Int_t a=mode, b=mode;
  bool l = true;
  Double_t center = hist[mode];
  while (center < frac) {
    center += ( (l = !l) ? hist[--a] : hist[++b] );
  }
  return {hist.GetBinLowEdge(l ? a+1 : a),
          hist.GetBinLowEdge(l ? b+1 : b)};
}

int main(int argc, char** argv)
{
  if (argc<4) {
    cout << "usage: " << argv[0]
         << " workspace.root output.pdf config_files" << endl;
    return 1;
  }
  const string ofname(argv[2]);
  assert_ext(argv[1],{"root"});
  assert_ext(ofname,{"pdf"});

  workspace ws(argv[1]);
  
  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  canv.SetLogy();
  canv.SaveAs((ofname+'[').c_str());
  
  TLine  line;
  TLatex text;
  text.SetTextSize(0.035);

  for (int i=3; i<argc; ++i) {
    cout << '\n' << argv[i] << endl;
    ifstream cfg(argv[i]);
    string var;
    double val;
    while (cfg >> var >> val) {
      ws.set_var(var,val);
    }
  
    TH1* hist = ws.m_yy_hist(10000);
    hist->SetTitle(argv[i]);
    hist->SetStats(false);
  
    hist->Draw();
    canv.Update();
    
    const auto fwhm = get_fwhm(*dynamic_cast<TH1F*>(hist));

    line.SetLineColor(2);
    line.DrawLine(
      fwhm.first,  pow(10.,gPad->GetUymin()),
      fwhm.first,  pow(10.,gPad->GetUymax()));
    line.DrawLine(
      fwhm.second, pow(10.,gPad->GetUymin()),
      fwhm.second, pow(10.,gPad->GetUymax()));
      
    text.DrawLatexNDC(0.15, 0.83,
      cat(setprecision(4),"FWHM = ",
          fwhm.second-fwhm.first," GeV").c_str()
    )->SetTextColor(2);
    
    const auto sigma_68_tails = get_sigma_tails(*dynamic_cast<TH1F*>(hist),0.68);

    line.SetLineColor(3);
    line.DrawLine(
      sigma_68_tails.first,  pow(10.,gPad->GetUymin()),
      sigma_68_tails.first,  pow(10.,gPad->GetUymax()));
    line.DrawLine(
      sigma_68_tails.second, pow(10.,gPad->GetUymin()),
      sigma_68_tails.second, pow(10.,gPad->GetUymax()));
      
    text.DrawLatexNDC(0.15, 0.78,
      cat(setprecision(4),"sigma_68_tails = ",
          (sigma_68_tails.second-sigma_68_tails.first)/2," GeV").c_str()
    )->SetTextColor(3);
    
    const auto sigma_68_mode = get_sigma_mode(*dynamic_cast<TH1F*>(hist),0.68);

    line.SetLineColor(4);
    line.DrawLine(
      sigma_68_mode.first,  pow(10.,gPad->GetUymin()),
      sigma_68_mode.first,  pow(10.,gPad->GetUymax()));
    line.DrawLine(
      sigma_68_mode.second, pow(10.,gPad->GetUymin()),
      sigma_68_mode.second, pow(10.,gPad->GetUymax()));
      
    text.DrawLatexNDC(0.15, 0.73,
      cat(setprecision(4),"sigma_68_mode = ",
          (sigma_68_mode.second-sigma_68_mode.first)/2," GeV").c_str()
    )->SetTextColor(4);
    
    const auto sigma_95_tails = get_sigma_tails(*dynamic_cast<TH1F*>(hist),0.95);

    line.SetLineColor(3);
    line.DrawLine(
      sigma_95_tails.first,  pow(10.,gPad->GetUymin()),
      sigma_95_tails.first,  pow(10.,gPad->GetUymax()));
    line.DrawLine(
      sigma_95_tails.second, pow(10.,gPad->GetUymin()),
      sigma_95_tails.second, pow(10.,gPad->GetUymax()));
      
    text.DrawLatexNDC(0.15, 0.68,
      cat(setprecision(4),"sigma_95_tails = ",
          (sigma_95_tails.second-sigma_95_tails.first)/2," GeV").c_str()
    )->SetTextColor(3);
    
    const auto sigma_95_mode = get_sigma_mode(*dynamic_cast<TH1F*>(hist),0.95);

    line.SetLineColor(4);
    line.DrawLine(
      sigma_95_mode.first,  pow(10.,gPad->GetUymin()),
      sigma_95_mode.first,  pow(10.,gPad->GetUymax()));
    line.DrawLine(
      sigma_95_mode.second, pow(10.,gPad->GetUymin()),
      sigma_95_mode.second, pow(10.,gPad->GetUymax()));
      
    text.DrawLatexNDC(0.15, 0.63,
      cat(setprecision(4),"sigma_95_mode = ",
          (sigma_95_mode.second-sigma_95_mode.first)/2," GeV").c_str()
    )->SetTextColor(4);

    canv.SaveAs(ofname.c_str());
    delete hist;
  }

  canv.SaveAs((ofname+']').c_str());
  
  return 0;
}
