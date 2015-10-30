// Developed by Ivan Pogrebnyak, MSU

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <initializer_list>
#include <memory>
#include <stdexcept>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
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
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooCurve.h>

#include "catstr.hh"
#include "senum.hh"
#include "val_err.hh"
#include "TGraph_fcns.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::pair;
using std::vector;
using std::map;
using std::unordered_set;
using std::initializer_list;
using std::unique_ptr;
using std::runtime_error;
using std::stringstream;
using std::scientific;
using std::setprecision;
using std::fixed;
namespace po = boost::program_options;

#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

#if GCC_VERSION < 40900 // Test for GCC < 4.9.0
#include <boost/regex.hpp>
using boost::regex;
using boost::smatch;
using boost::regex_match;
#define regex_icase boost::regex::icase
#else
#include <regex>
using std::regex;
using std::smatch;
using std::regex_match;
using std::regex_constants::icase;
#define regex_icase std::regex_constants::icase
#endif

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

template<typename T> inline T sq(T x) noexcept { return x*x; }

namespace std {
  template <typename T1, typename T2>
  istream& operator>>(istream& in, pair<T1,T2>& p) {
    string s;
    in >> s;
    size_t sep = s.find(':');
    if (sep==string::npos) throw invalid_argument(
      cat('\"',s,"\": pair values must be delimited by \':\'"));
    stringstream (s.substr(0,sep)) >> p.first;
    stringstream (s.substr(sep+1)) >> p.second;
    return in;
  }
}

template<typename T, typename Key=std::string>
class seqmap: public std::vector<std::pair<Key,T>> {
public:
  T& operator[](const Key& key) {
    auto it = std::find_if(this->begin(),this->end(),
      [&key](const std::pair<Key,T>& p){ return p.first==key; });
    if (it==this->end()) {
      this->emplace_back(key,T());
      return this->back().second;
    } else return it->second;
  }
};

template<typename T>
inline T* get(TDirectory* d, const char* name) {
  TObject *obj = d->Get(name);
  if (!obj) {
    throw runtime_error( cat(
      "No object ",name," in ",d->GetName()
    ) );
  } else if (obj->InheritsFrom(T::Class())) {
    return static_cast<T*>(obj);
  } else {
    throw runtime_error( cat(
      obj->ClassName(), ' ', obj->GetName(),
      " does not inherit from ", T::Class()->GetName()
    ) );
  }
}

// options ------------------
senum(Fit,(none)(gaus)(cb))
Fit::type fit;

string ofname, cfname, wfname;
vector<string> ifname;
vector<Color_t> colors;
Int_t nbins;
pair<double,double> xrange;
bool logy, fix_alpha;
int prec;
vector<pair<string,pair<double,double>>> new_ws_ranges;
// --------------------------

// global -------------------
seqmap<seqmap<val_err<double>>> stats;
TTree *tree;
TCanvas *canv;
TLatex *lbl;
// --------------------------

using FitResult = unique_ptr<RooFitResult>;

class workspace {
  TFile *file;
  RooWorkspace *ws;
  RooSimultaneous *sim_pdf;
  RooCategory *rcat;
  RooRealVar *myy;

public:
  workspace(const string& fname)
  : file(new TFile(fname.c_str(),"read")),
    ws(get<RooWorkspace>(file,"mcfit")),
    sim_pdf(static_cast<RooSimultaneous*>(ws->obj("mc_sim_pdf_bin0"))),
    rcat(static_cast<RooCategory*>(ws->obj("mc_sample"))),
    myy(static_cast<RooRealVar*>(ws->var("m_yy")))
  {
    for (const auto& range : new_ws_ranges)
      ws->var(range.first.c_str())
        ->setRange(range.second.first,range.second.second);
  }

  ~workspace() {
    delete myy;
    delete rcat;
    delete sim_pdf;
    delete file;
  }

  inline RooWorkspace* operator->() noexcept { return ws; }

  pair<FitResult,RooCurve*> fit(TH1* hist) const {
    // Produce a RooDataHist object from the TH1
    RooDataHist *rdh = new RooDataHist(
      "mc_dh","mc_dh",RooArgSet(*myy),hist);

    // This map might look a bit useless,
    // but the PDF is designed to fit several datasets simultaneously.
    // This is done by assigning the PDF categories,
    // here we only have one histogram
    map<string,RooDataHist*> rdhmap;
    rdhmap["mc_125"] = rdh;

    // With the map we can build one combined dataset
    // that has the proper link of the category information.
    RooDataHist crdh("c_mc_dh","c_mc_dh",
      RooArgSet(*myy), RooFit::Index(*rcat), RooFit::Import(rdhmap));

    // Now we are ready to fit! We have a PDF and a RooDataHist
    RooFitResult *res = sim_pdf->fitTo(crdh,
      RooFit::Extended(false),
      RooFit::InitialHesse(true),
      RooFit::SumW2Error(true),
      RooFit::Save(true),
      RooFit::NumCPU(4),
      RooFit::Minimizer("Minuit2"),
      RooFit::Offset(true),
      RooFit::Strategy(2)
    );

    RooPlot *frame = myy->frame();
    crdh.plotOn(frame,
      RooFit::LineColor(12),
      RooFit::Cut("mc_sample == mc_sample::mc_125")
    );
    sim_pdf->plotOn(frame,
      RooFit::LineColor(85),
      RooFit::Slice(*rcat, "mc_125"),
      RooFit::ProjWData(RooArgSet(*rcat), crdh)
    );
    auto *curve = frame->getCurve();
    //frame->Draw("same");
    curve->Draw("same");
    res->Print("v");

    delete rdh;
    return {FitResult(res),curve};
  }
};

inline Double_t get_FWHM(const TH1* hist) noexcept {
  const Double_t half_max = hist->GetMaximum()/2;
  return hist->GetBinCenter(hist->FindLastBinAbove(half_max))
       - hist->GetBinCenter(hist->FindFirstBinAbove(half_max));
}
inline Double_t get_FWHM(const TGraph* gr) noexcept {
  const Double_t half_max = max(gr).second/2;
  return rfindx(gr,half_max) - lfindx(gr,half_max);
}

Double_t mean_window(const TH1* hist, Double_t a, Double_t b) noexcept {
  Int_t ai = hist->FindFixBin(a);
  Int_t bi = hist->FindFixBin(b);
  Double_t mean = 0., stdev = 0., sumw = 0.;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    const Double_t x = hist->GetBinCenter(i);
    sumw += w;
    mean += x*w;
  }
  mean /= sumw;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    const Double_t x = hist->GetBinCenter(i);
    stdev += sq(x-mean)*w;
  }
  stdev = sqrt(stdev/sumw);

  ai = hist->FindFixBin(mean-1.5*stdev);
  bi = hist->FindFixBin(mean+2.0*stdev);
  mean = 0.; sumw = 0.;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    sumw += w;
    mean += hist->GetBinCenter(i) * w;
  }
  return mean/sumw;
}

void make_hist(TH1*& hist, const char* name, const string& proc,
               double scale, const char* branch) {
  const string cmd1(cat(
    branch,"/1000>>hist(",nbins,",",xrange.first,",",xrange.second,")"));
  const string cmd2(
    "HGamEventInfoAuxDyn.crossSectionBRfilterEff"
    "*HGamEventInfoAuxDyn.weight"
    "*(HGamEventInfoAuxDyn.isPassed==1)");
  tree->Draw(cmd1.c_str(),cmd2.c_str());
  cout << endl << name
       << endl << cmd1
       << endl << cmd2 << endl;
  TH1 *temp = get<TH1>(gDirectory,"hist");
  temp->Scale(1000.*scale/temp->GetBinWidth(1));

  stats[name]["xsec_"+proc] = temp->Integral(0,temp->GetNbinsX()+1,"width");

  if (!hist)
    (hist = (TH1*)temp->Clone(name))->SetDirectory(0);
  else hist->Add(temp);
}

vector<FitResult> draw(const initializer_list<TH1*>& hs) {

  vector<FitResult> res;
  if (fit==Fit::cb) res.reserve(hs.size());

  int i=0;
  Color_t color;
  for (auto it=hs.begin(), end=hs.end(); it!=end; ++it) {
    TH1 *hist = *it;
    const char *name = hist->GetName();
    hist->SetStats(false);
    hist->SetLineWidth(2);
    color = colors[(i++) % colors.size()];
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);

    if (it==hs.begin()) {
      hist->SetXTitle("m_{#gamma#gamma} [GeV]");
      hist->SetYTitle("d#sigma/dm_{#gamma#gamma} [fb/GeV]");
      hist->SetTitleOffset(1.3,"Y");
      hist->Draw();
      lbl->DrawLatex(.24,0.88-0.04*i,"Entries");
      lbl->DrawLatex(.32,0.88-0.04*i,"mean");
      lbl->DrawLatex(.40,0.88-0.04*i,"stdev");
    } else hist->Draw("same");

    Double_t mean  = hist->GetMean(),
             stdev = hist->GetStdDev();

    auto &hstat = stats[name];
    hstat["hist_N"] = hist->GetEntries();
    hstat["hist_mean"] = {mean,hist->GetMeanError()};
    hstat["hist_stdev"] = {stdev,hist->GetStdDevError()};

    hstat["hist_window_mean"] = mean_window(hist,120,130);

    switch (fit) { // FITTING +++++++++++++++++++++++++++++++++++++++
      case Fit::none: break;

      case Fit::gaus: {
        auto fit_res = hist->Fit("gaus","S");
        mean  = fit_res->Value(1);
        stdev = fit_res->Value(2);

        auto &hstat = stats[name];
        hstat["gaus_mean"] = {mean,fit_res->Error(1)};
        hstat["gaus_stdev"] = {stdev,fit_res->Error(2)};
      break; }

      case Fit::cb: {
        cout << "\033[32mFitting " << name << "\033[0m" << endl;
        static workspace ws(wfname);
        auto fit_res = ws.fit(hist);

        auto &hstat = stats[name];
        for (const string& varname : {
          "crys_alpha_bin0", "crys_norm_bin0", "fcb_bin0", "gaus_kappa_bin0",
          "gaus_mean_offset_bin0", "mean_offset_bin0", "sigma_offset_bin0"
        }) {
          auto *var = static_cast<RooRealVar*>(
            fit_res.first->floatParsFinal().find(varname.c_str()));
          hstat[varname] = {var->getVal(),var->getError()};
        }
        hstat["FWHM"] = get_FWHM(fit_res.second);

        if (fix_alpha) if (!strcmp(name,"nominal")) {
          auto *alpha = ws->var("crys_alpha_bin0");
          alpha->setRange(alpha->getVal(),alpha->getVal());
        }

/*
        auto Nsig = static_cast<RooRealVar*>(
          fit_res->floatParsFinal().find("NSig_bin0")
        );
        const double Nsig_rat = Nsig->getVal()/Nsig->getError();
        hstat["p0"] = Nsig_rat > 0 ?
                      0.5*TMath::Prob( sq(Nsig_rat) , 1.) :
                      1 - 0.5*TMath::Prob( sq(Nsig_rat) , 1.);
*/

      res.emplace_back(move(fit_res.first));

      break; }
    }

    auto lblp = lbl->DrawLatex(.12,0.84-0.04*i,hist->GetName());
    lblp->SetTextColor(color);
    lblp->DrawLatex(.24,0.84-0.04*i,cat(hist->GetEntries()).c_str());
    lblp->DrawLatex(.32,0.84-0.04*i,
      cat(fixed,setprecision(2),mean).c_str());
    lblp->DrawLatex(.40,0.84-0.04*i,
      cat(fixed,setprecision(2),stdev).c_str());
  }
  canv->SaveAs(ofname.c_str());

  return res;
}

int main(int argc, char** argv)
{
  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->multitoken()->required(),
       "*input root file names")
      ("output,o", po::value(&ofname)->required(),
       "*output pdf file name")
      ("config,c", po::value(&cfname),
       "configuration file name")

      ("logy,l", po::bool_switch(&logy),
       "logarithmic Y axis")
      ("fit,f", po::value(&fit)->default_value(Fit::none),
       cat("fit type: ",Fit::_str_all()).c_str())
      ("workspace,w", po::value(&wfname)->default_value("data/ws.root"),
       "ROOT file with RooWorkspace for CB fits")
      ("fix-alpha", po::bool_switch(&fix_alpha),
       "fix crys_alpha_bin0 parameter after nominal fit")
      ("xrange,x", po::value(&xrange)->default_value({105,140},"105:140"),
       "histograms\' X range")
      ("nbins,n", po::value(&nbins)->default_value(100),
       "histograms\' number of bins")
      ("prec", po::value(&prec)->default_value(-1),
       "summary table precision, -1 prints uncertainty")
      ("colors", po::value(&colors)->multitoken()->
        default_value(decltype(colors)({602,46}), "{602,46}"),
       "histograms\' colors")

      ("ws-setRange", po::value(&new_ws_ranges),
       "call RooWorkspace::setRange()")
    ;

    po::positional_options_description pos;
    pos.add("input",-1);

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

  TH1 *nom=nullptr,
      *scale_down=nullptr, *scale_up=nullptr,
      *res_down=nullptr,   *res_up=nullptr;

  // LOOP over input files
  for (const string& f : ifname) {
    TFile *file = new TFile(f.c_str(),"read");
    if (file->IsZombie()) return 1;
    cout << "Data file: " << f << endl;
    tree = get<TTree>(file,"CollectionTree");

    const size_t slash = f.find('/')+1;
    const double xsecscale = 1./get<TH1>(file,
      ("CutFlow_"+f.substr(slash,f.find('.')-slash)+"_weighted").c_str()
    )->GetBinContent(3);

    // Regex for process identification
    static regex proc_re(".*[\\._]?(gg.|VBF|ttH|WH|ZH)[0-9]*[\\._]?.*",
                         regex_icase);
    smatch proc_match;
    if (!regex_match(f, proc_match, proc_re))
      throw runtime_error(cat("Filename \"",f,"\" does not specify process"));
    const string proc(proc_match.str(1));

    // protect from repeated processes
    static unordered_set<string> procs;
    if (!procs.emplace(proc).second)
      throw runtime_error(cat("File \"",f,"\" repeats process ",proc));

    // Make or add histograms
    make_hist(nom,"nominal",proc,xsecscale,
              "HGamEventInfoAuxDyn.m_yy");
    make_hist(scale_down,"scale_down",proc,xsecscale,
              "HGamEventInfo_EG_SCALE_ALL__1downAuxDyn.m_yy");
    make_hist(scale_up,"scale_up",proc,xsecscale,
              "HGamEventInfo_EG_SCALE_ALL__1upAuxDyn.m_yy");
    make_hist(res_down,"res_down",proc,xsecscale,
              "HGamEventInfo_EG_RESOLUTION_ALL__1downAuxDyn.m_yy");
    make_hist(res_up,"res_up",proc,xsecscale,
              "HGamEventInfo_EG_RESOLUTION_ALL__1upAuxDyn.m_yy");

    delete file;
  }
  cout << endl;

  // ---------------------------------------

  canv = new TCanvas();
  canv->SetMargin(0.1,0.04,0.1,0.1);
  if (logy) canv->SetLogy();
  canv->SaveAs((ofname+'[').c_str());

  lbl = new TLatex();
  lbl->SetTextFont(43);
  lbl->SetTextSize(15);
  lbl->SetNDC();

  TH2 *corr = draw({nom}).front()
    ->correlationHist("Nominal signal fit correlation matrix");

  if (logy) canv->SetLogy(false);
  gStyle->SetPaintTextFormat(".3f");
  canv->SetMargin(0.17,0.12,0.1,0.1);
  corr->SetStats(false);
  corr->SetMarkerSize(1.8);
  corr->Draw("COLZ TEXT");
  canv->SaveAs(ofname.c_str());
  canv->SetMargin(0.1,0.04,0.1,0.1);
  if (logy) canv->SetLogy(true);

  draw({scale_down,scale_up});
  draw({res_down,res_up});

  // Print summary page *********************************************
  {
    canv->Clear();
    vector<unique_ptr<TPaveText>> txt;
    int i, n=6;
    txt.reserve(n);
    for (i=0; i<n; ++i) {
      TPaveText *pt;
      txt.emplace_back(pt = new TPaveText(float(i)/n,0.,float(i+1)/n,1.,"NBNDC"));
      pt->SetFillColor(0);
    }
    i=1;
    int m = 1;
    txt[0]->AddText("");
    for (auto& hist : stats) {
      txt[i]->AddText(hist.first.c_str());
      for (auto& var : hist.second) {
        if (i==1) {
          txt[0]->AddText(var.first.c_str());
          ++m;
        }
        stringstream ss;
        if (var.first.substr(0,var.first.find('_'))=="xsec") {
          ss << setprecision(3) << var.second.val;
        } else if (var.first.substr(var.first.rfind('_')+1)=="N") {
          ss << fixed << setprecision(0) << var.second.val;
        } else {
          if (prec==-1) var.second.print(ss," #pm ");
          else ss << scientific << setprecision(prec) << var.second.val;
        }
        txt[i]->AddText(ss.str().c_str());
      }
      ++i;
    }
    TLine *line = new TLine();
    for (i=0; i<n; ++i) {
      txt[i]->Draw();
      if (i) line->DrawLineNDC(float(i)/n,0.,float(i)/n,1.);
    }
    for (i=1; i<m; ++i) line->DrawLineNDC(0.,float(i)/m,1.,float(i)/m);
    canv->SaveAs(ofname.c_str());
  }
  // ****************************************************************

  if (fit==Fit::cb) {
    canv->Clear();
    double scale      = stats["nominal"   ]["mean_offset_bin0"].val;
    double scale_down = stats["scale_down"]["mean_offset_bin0"].val;
    double scale_up   = stats["scale_up"  ]["mean_offset_bin0"].val;
           scale_down = scale_down - scale;
           scale_up   = scale_up   - scale;
    double scale_sym  = (scale_up-scale_down)/2;

    double win        = stats["nominal"   ]["hist_window_mean"].val;
    double win_down   = stats["scale_down"]["hist_window_mean"].val;
    double win_up     = stats["scale_up"  ]["hist_window_mean"].val;
           win_down   = win_down - win;
           win_up     = win_up   - win;
    double win_sym    = (win_up-win_down)/2;

    double res        = stats["nominal"   ]["sigma_offset_bin0"].val;
    double res_down   = stats["res_down"  ]["sigma_offset_bin0"].val;
    double res_up     = stats["res_up"    ]["sigma_offset_bin0"].val;
           res_down   = res_down - res;
           res_up     = res_up   - res;
    double res_sym    = (res_up-res_down)/2;

    double fwhm       = stats["nominal"   ]["FWHM"].val;
    double fwhm_down  = stats["res_down"  ]["FWHM"].val;
    double fwhm_up    = stats["res_up"    ]["FWHM"].val;
           fwhm_down  = fwhm_down - fwhm;
           fwhm_up    = fwhm_up   - fwhm;
    double fwhm_sym   = (fwhm_up-fwhm_down)/2;

    vector<unique_ptr<TPaveText>> txt;
    txt.reserve(6);
    for (int i=0; i<6; ++i) {
      TPaveText *pt;
      txt.emplace_back(pt = new TPaveText(i/6.,0.,(i+1)/6.,1.,"NBNDC"));
      pt->SetFillColor(0);
    }

    txt[0]->AddText("[GeV]");
    txt[1]->AddText("Scale");
    txt[2]->AddText("Window");
    txt[3]->AddText("Resolution");
    txt[4]->AddText("HWHM");
    txt[5]->AddText("FWHM/FWHM_{nom}");

    txt[0]->AddText("Nominal");
    txt[1]->AddText(Form("%.3f",scale));
    txt[2]->AddText(Form("%.3f",win));
    txt[3]->AddText(Form("%.3f",res));
    txt[4]->AddText(Form("%.3f",fwhm/2));
    txt[5]->AddText("1");

    txt[0]->AddText("Variation");
    txt[1]->AddText(Form("%.3f, +%.3f",scale_down,scale_up));
    txt[2]->AddText(Form("%.3f, +%.3f",win_down,win_up));
    txt[3]->AddText(Form("%.3f, +%.3f",res_down,res_up));
    txt[4]->AddText(Form("%.3f, +%.3f",fwhm_down/2,fwhm_up/2));
    txt[5]->AddText(Form("%.3f, +%.3f",fwhm_down/fwhm,fwhm_up/fwhm));

    txt[0]->AddText("Average Variation");
    txt[1]->AddText(Form("#pm %.3f",scale_sym));
    txt[2]->AddText(Form("#pm %.3f",win_sym));
    txt[3]->AddText(Form("#pm %.3f",res_sym));
    txt[4]->AddText(Form("#pm %.3f",fwhm_sym/2));
    txt[5]->AddText(Form("#pm %.3f",fwhm_sym/fwhm));

    TLine *line = new TLine();

    for (int i=0; i<6; ++i) {
      txt[i]->Draw();
      if (i) line->DrawLineNDC(i/6.,0.,i/6.,1.);
    }
    for (int i=1; i<4; ++i) line->DrawLineNDC(0.,i/4.,1.,i/4.);
    canv->SaveAs(ofname.c_str());
  }

  canv->SaveAs((ofname+']').c_str());
  delete canv;
  delete lbl;

  return 0;
}
