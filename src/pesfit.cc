#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
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

string ofname, cfname;
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
  vector<unique_ptr<TLatex>> lblp;
  lblp.reserve(15);
  auto latex = [&lblp](TLatex* lbl, Double_t x, Double_t y, const string& text) {
    auto *p = lbl->DrawLatex(x,y,text.c_str());
    lblp.emplace_back(p);
    return p;
  };
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
      latex(lbl,.67,0.88-0.04*i,"Entries");
      latex(lbl,.77,0.88-0.04*i,"mean");
      latex(lbl,.87,0.88-0.04*i,"stdev");
    } else hist->Draw("same");

    Double_t mean  = hist->GetMean(),
             stdev = hist->GetStdDev();

    stats[name]["hist_mean"] = mean;
    stats[name]["hist_mean_err"] = hist->GetMeanError();
    stats[name]["hist_stdev"] = stdev;
    stats[name]["hist_stdev_err"] = hist->GetStdDevError();

    switch (fit) {
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
        cerr << "cb fit not implemented yet" << endl;
        break; }
    }

    auto lblp = latex(lbl,.52,0.84-0.04*i,hist->GetName());
    lblp->SetTextColor(color);
    latex(lblp,.67,0.84-0.04*i,cat(hist->GetEntries()));
    latex(lblp,.77,0.84-0.04*i,
      cat(fixed,setprecision(2),mean));
    latex(lblp,.87,0.84-0.04*i,
      cat(fixed,setprecision(2),stdev));
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
