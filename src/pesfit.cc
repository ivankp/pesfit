#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <memory>
#include <regex>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>

#include "catstr.hh"
#include "senum.hh"

using namespace std;
namespace po = boost::program_options;

template<typename T> using usmap = unordered_map<string,T>;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

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
bool print_stats, logy;
// --------------------------

// global -------------------
usmap<usmap<double>> stats;
TTree *tree;
TCanvas *canv;
TLatex *lbl;
// --------------------------

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
  { }
  
  ~workspace() {
    delete myy;
    delete rcat;
    delete sim_pdf;
    delete file;
  }
  
  unique_ptr<RooFitResult> fit(TH1* hist) const {
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
    frame->Draw("same");
    
    delete rdh;
    return unique_ptr<RooFitResult>(res);
  }
};

void make_hist(TH1*& hist, const char* name, const string& proc,
               double scale, const char* branch) {
  const string cmd1(cat(branch,"/1000>>hist(",nbins,",105,160)"));
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

void draw(const initializer_list<TH1*>& hs) {
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
      lbl->DrawLatex(.67,0.88-0.04*i,"Entries");
      lbl->DrawLatex(.77,0.88-0.04*i,"mean");
      lbl->DrawLatex(.87,0.88-0.04*i,"stdev");
    } else hist->Draw("same");

    Double_t mean  = hist->GetMean(),
             stdev = hist->GetStdDev();

    stats[name]["hist_mean"] = mean;
    stats[name]["hist_mean_err"] = hist->GetMeanError();
    stats[name]["hist_stdev"] = stdev;
    stats[name]["hist_stdev_err"] = hist->GetStdDevError();

    switch (fit) { // FITTING +++++++++++++++++++++++++++++++++++++++
      case Fit::none: break;

      case Fit::gaus: {
        auto fr = hist->Fit("gaus","S");
        mean  = fr->Value(1);
        stdev = fr->Value(2);

        stats[name]["gaus_mean"] = mean;
        stats[name]["gaus_mean_err"] = fr->Error(1);
        stats[name]["gaus_stdev"] = stdev;
        stats[name]["gaus_stdev_err"] = fr->Error(2);
      break; }

      case Fit::cb: {
        cout << "\033[32mFitting " << name << "\033[0m" << endl;
        static workspace ws(wfname);
        ws.fit(hist)->Print("v");
      break; }
    }

    auto lblp = lbl->DrawLatex(.52,0.84-0.04*i,hist->GetName());
    lblp->SetTextColor(color);
    lblp->DrawLatex(.67,0.84-0.04*i,cat(hist->GetEntries()).c_str());
    lblp->DrawLatex(.77,0.84-0.04*i,
      cat(fixed,setprecision(2),mean).c_str());
    lblp->DrawLatex(.87,0.84-0.04*i,
      cat(fixed,setprecision(2),stdev).c_str());
  }
  canv->SaveAs(ofname.c_str());
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
      ("stats,s", po::bool_switch(&print_stats),
       "print interesting values for histograms")
      ("nbins,n", po::value(&nbins)->default_value(100),
       "histograms\' number of bins")
      ("colors", po::value(&colors)->multitoken()->
        default_value(decltype(colors)({602,46}), "{602,46}"),
       "histograms\' colors")
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
  } catch (exception& e) {
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
                         regex_constants::icase);
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
  lbl->SetTextSize(20);
  lbl->SetNDC();

  draw({nom});
  draw({scale_down,scale_up});
  draw({res_down,res_up});

  canv->SaveAs((ofname+']').c_str());
  delete canv;
  delete lbl;
  
  if (print_stats) {
    for (auto& stat : stats) {
      cout << endl << stat.first << endl;
      for (auto& var : stat.second) {
        cout <<"  "<< var.first <<" = "<< var.second << endl;
      }
    }
  }

  return 0;
}
