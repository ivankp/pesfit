#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <unordered_set>
#include <memory>
#include <stdexcept>
#include <cmath>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>

#include "regex.hh"
#include "catstr.hh"
#include "root_safe_get.hh"
#include "TGraph_fcns.hh"
#include "binned.hh"
#include "workspace.hh"
#include "golden_min.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#ifndef Old8TeVFile
#define Old8TeVFile 0
#endif

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args ) {
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

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

class branches_manager {
  TTree *tree;
  vector<string> names;
public:
  branches_manager(TTree* tree): tree(tree) { }
  template<typename T>
  void operator()(string&& branch, T* x) {
    tree->SetBranchAddress(branch.c_str(), x);
    names.emplace_back(branch);
  }
  ~branches_manager() {
    tree->SetBranchStatus("*",0);
    for (const auto& name : names)
      tree->SetBranchStatus(name.c_str(),1);
  }
};

int main(int argc, char** argv)
{
  vector<string> ifname;
  string ofname, wfname, cfname;
  bool logy, sumw2;
  int nbins;
  pair<double,double> xrange;
  vector<Color_t> colors;
  vector<pair<string,pair<double,double>>> new_ws_ranges;
  double sigma_frac;
  pair<int,pair<double,double>> vert;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->multitoken()->required(),
       "*input root file names")
      ("output,o", po::value(&ofname)->required(),
       "*output pdf or root file name")
      ("workspace,w", po::value(&wfname)->required(),
       "ROOT file with RooWorkspace for CB fits")
      ("config,c", po::value(&cfname),
       "configuration file name")

      ("sigma-frac,s", po::value(&sigma_frac)->default_value(0.68,"0.68"),
       "confidence interval fraction")
      ("vert,v", po::value(&vert)->required(),
       "num vertices binning")
      ("nbins,b", po::value(&nbins)->default_value(100),
       "histograms\' number of bins")
      ("xrange,x", po::value(&xrange)->default_value({105,140},"105:140"),
       "histograms\' X range")
      ("logy,l", po::bool_switch(&logy),
       "logarithmic Y axis")
      ("sumw2", po::bool_switch(&sumw2),
       "histograms' Sumw2")
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

    const string ofext = ofname.substr(ofname.rfind('.')+1);
    if (ofext!="pdf") throw runtime_error(
      "Output file extension "+ofext+" is not pdf"
    );

  } catch (exception& e) {
    cerr << "\033[31mArgs: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  // Book histograms ************************************************
  constexpr auto hist_types = {"selected"};
  binned<vector<pair<TH1*,unique_ptr<TH1>>>>
    hmap(vert.first,vert.second.first,vert.second.second);

  for (size_t i=0, n=hmap.nbins()+1; i<=n; ++i) {
    auto& hs = hmap.at(i);
    hs.reserve(hist_types.size());
    for (const auto *hist_type : hist_types) {
      string name;
      if (i==0) name = cat(
        hist_type,"_m_yy_nvert<",hmap.right_edge(i)
      ); else if (i==n) name = cat(
        hist_type,"_m_yy_nvert>",hmap.left_edge(i)
      ); else name = cat(
        hist_type,"_m_yy_nvert[",hmap.left_edge(i),',',hmap.right_edge(i),')'
      );
      hs.emplace_back(
        new TH1D( name.c_str(),
          ";m_{#gamma#gamma} [GeV];d#sigma/dm_{#gamma#gamma} [fb/GeV]",
          nbins,xrange.first,xrange.second),
        make_unique<TH1D>("tmp_hist","",nbins,xrange.first,xrange.second)
      );
      hs.back().second->SetDirectory(0);
      hs.back().first ->Sumw2(sumw2);
      hs.back().second->Sumw2(sumw2);
    }
  }

  cout << "Histograms:" << endl;
  for (const auto& hs : hmap)
    for (const auto& h : hs)
      cout << "  " << h.first->GetName() << endl;
  cout << endl;

  // LOOP over input files ******************************************
  for (const string& f : ifname) {
    TFile *file = new TFile(f.c_str(),"read");
    if (file->IsZombie()) return 1;
    cout << "Data file: " << f << endl;

    #if !(Old8TeVFile)
    TTree *tree = get<TTree>(file,"CollectionTree");
    
    const size_t slash = f.rfind('/')+1;
    const double xsecscale = 1e3/get<TH1>(file,
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

    #else
    TTree *tree = get<TTree>(file,"Hgg_tree");
    const double xsecscale = 1.;
    #endif

    // Branch variables
    #if !(Old8TeVFile)
    static array<pair<Float_t,Int_t>,hist_types.size()> var;
    static Float_t crossSectionBRfilterEff;//, weight;
    static Char_t isPassed;
    {
      branches_manager br(tree);

      br("HGamEventInfoAuxDyn.m_yy", &var.at(0).first);
      br("HGamEventInfoAuxDyn.numberOfPrimaryVertices", &var.at(0).second);

      br("HGamEventInfoAuxDyn.crossSectionBRfilterEff",
         &crossSectionBRfilterEff);
      // br("HGamEventInfoAuxDyn.weight", &weight);
      br("HGamEventInfoAuxDyn.isPassed", &isPassed);
    }
    #else
    static array<pair<Double_t,Int_t>,hist_types.size()> var;
    //static Double_t weight;
    static Int_t truth_mH;
    {
      branches_manager br(tree);

      br("m_yy", &var.at(0).first);
      br("NPV",  &var.at(0).second);
      br("Higgs_truth_mass", &truth_mH);
      //br("weight", &weight);
    }
    #endif

    // LOOP over tree entries
    for (Long64_t nent=tree->GetEntries(), ent=0; ent<nent; ++ent) {
      tree->GetEntry(ent);

      for (size_t i=0; i<var.size(); ++i) {
        #if !(Old8TeVFile)
        if (isPassed==1)
          hmap[var[i].second][i].second->Fill(
            var[i].first/1e3,
            crossSectionBRfilterEff//*weight
          );
        #else
        if (truth_mH==125)
          hmap[var[i].second][i].second->Fill(
            var[i].first/1e3 //, weight
          );
        #endif
      }
    }

    // Add histograms
    for (auto it=hmap.begin(), end=hmap.end(true); it!=end; ++it) {
      for (const auto& h : *it) {
        h.second->Scale(xsecscale/h.second->GetBinWidth(1));
        h.first->Add(h.second.get());
        h.second->Reset();
      }
    }

    delete file;
  }

  // Fit functions **************************************************
  workspace ws(wfname);
  golden_min gm;

  vector<array<tuple<TGraph*,double,double>,hist_types.size()>> fits;
  fits.reserve(hmap.nbins());

  for (const auto& hs : hmap) {
    size_t i=0;
    fits.emplace_back();
    for (const auto& h : hs) {
      TGraph* fit_gr = ws.fit(h.first).second;
      const Double_t integral = integrate(fit_gr);
      // minimize sigma
      Double_t x1 = gm( [fit_gr,sigma_frac,integral](double x1) {
        return intervalx2(fit_gr, sigma_frac, x1, integral) - x1;
      }, firstx(fit_gr), rtailx(fit_gr,sigma_frac,integral) ).first;
      Double_t x2 = intervalx2(fit_gr, sigma_frac, x1, integral);

      fit_gr->SetName(cat(h.first->GetName(),"_fit").c_str());
      fits.back()[i++] = make_tuple(fit_gr,x1,x2);
    }
  }

  // Draw histograms ************************************************
  TCanvas canv;
  canv.SetMargin(0.08,0.04,0.1,0.02);
  canv.SetTicks();

  TLatex lbl;
  lbl.SetTextFont(43);
  lbl.SetTextSize(18);
  lbl.SetNDC();

  TLine line;

  canv.SaveAs((ofname+'[').c_str());

  // Summary plot
  array<TH1*,hist_types.size()> hists; hists.fill(nullptr);
  for (const auto& fit : fits) {
    static int b=0;
    int h=0;
    for (const auto *hist_type : hist_types) {
      if (b==0) hists[h] = new TH1D(
        hist_type,
        cat(";Number primary vertices;"
            "#sigma_{",sigma_frac*100,"} [GeV]").c_str(),
        hmap.nbins(),hmap.get_bins().data());
      hists[h]->SetBinContent(b+1,(get<2>(fit[h])-get<1>(fit[h]))/2.);
      ++h;
    }
    ++b;
  }

  TLegend leg(0.12,0.67,0.45,0.76);
  array<const char*,hist_types.size()> leg_lbl {
    "Selected vertex"
  };

  for (auto* h : hists) {
    static int i=0;
    Color_t color = colors[i % colors.size()];
    h->SetStats(false);
    h->SetLineWidth(2);
    h->GetYaxis()->SetTitleOffset(1.05);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    // h->SetMarkerStyle(20+i);
    h->Draw(i ? "same" : "");
    if (i==0) {
      leg.AddEntry(h,leg_lbl[i]);
      double lxmin = 0.12;
      double ly = 0.9;
      lbl.DrawLatex(lxmin,ly,"ATLAS")->SetTextFont(73);
      lbl.DrawLatex(lxmin+0.095,ly,"Internal");
      lbl.DrawLatex(lxmin,ly-=0.06,"#it{#sqrt{s}} = 13 TeV");
      lbl.DrawLatex(lxmin,ly-=0.06,
        "#it{H#rightarrow#gamma#gamma}, #it{m_{H}} = 125 GeV");
    }
    ++i;
  }

  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.Draw();

  canv.SaveAs(ofname.c_str());

  // Histograms in vertex bins
  if (logy) canv.SetLogy();

  for (size_t b=1, n=hmap.nbins(); b<=n+1; ++b) {
    Color_t color;
    int i=0;
    for (const auto& hp : hmap.at(b)) {
      TH1* h = hp.first;
      h->SetStats(false);
      h->SetLineWidth(2);
      h->GetYaxis()->SetTitleOffset(1.05);
      h->SetLineColor(color = colors[i % colors.size()]);
      h->Draw(i ? "same" : "");
      const auto& fit = fits[b-1][i];
      if (b<=n) {
        get<0>(fit)->Draw("same");
        canv.Update();
        double y1 = canv.GetUymin();
        double y2 = canv.GetUymax();
        if (logy) {
          y1 = pow(10.,y1);
          y2 = pow(10.,y2);
        }
      	line.DrawLine(get<1>(fit), y1, get<1>(fit), y2)->SetLineColor(color);
      	line.DrawLine(get<2>(fit), y1, get<2>(fit), y2)->SetLineColor(color);

        canv.Update();
        TLine line(canv.GetUxmax(),canv.GetUymin(),canv.GetUxmin(),canv.GetUymax());
        line.Draw();
      }

      lbl.DrawLatex(0.12,0.9-0.05*i,h->GetName())->SetTextColor(color);
      lbl.DrawLatex(0.80,0.9-0.05*i,
        cat(h->GetEntries()).c_str())->SetTextColor(color);
      ++i;
    }
    canv.SaveAs(ofname.c_str());
  }

  canv.SaveAs((ofname+']').c_str());

  return 0;
}
