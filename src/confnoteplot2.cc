#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

#include <boost/program_options.hpp>

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
#include "workspace.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

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

int main(int argc, char** argv)
{
  vector<string> ifname;
  string ofname, wfname, cfname;
  bool logy;
  int nbins;
  pair<double,double> xrange;
  vector<Color_t> colors;
  vector<pair<string,pair<double,double>>> new_ws_ranges;

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

      ("logy,l", po::bool_switch(&logy),
       "logarithmic Y axis")
      ("nbins,n", po::value(&nbins)->default_value(100),
       "histograms\' number of bins")
      ("xrange,x", po::value(&xrange)->default_value({105,140},"105:140"),
       "histograms\' X range")
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
  hmap({0,5,10,15,20,25,30});

  for (size_t i=1, n=hmap.nbins(); i<=n; ++i) {
    static auto& hs = hmap.at(i);
    hs.reserve(hist_types.size());
    for (const auto *hist_type : hist_types) {
      hs.emplace_back(
        new TH1D( cat(
          hist_type,"_m_yy_nvert[",hmap.left_edge(i),',',hmap.right_edge(i),')'
        ).c_str(),
        ";m_{#gamma#gamma} [GeV];d#sigma/dm_{#gamma#gamma} [fb/GeV]",
        nbins,xrange.first,xrange.second),
        make_unique<TH1D>("tmph","",nbins,xrange.first,xrange.second)
      );
      hs.back().second->SetDirectory(0);
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
    TTree *tree = get<TTree>(file,"CollectionTree");

    const size_t slash = f.rfind('/')+1;
    const double xsecscale = 1000./get<TH1>(file,
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

    // Branch variables
    static vector<pair<Float_t,Int_t>> var(hist_types.size());
    tree->SetBranchAddress(
      "HGamEventInfoAuxDyn.m_yy",&var.at(0).first);
    tree->SetBranchAddress(
      "HGamEventInfoAuxDyn.numberOfPrimaryVertices",&var.at(0).second);

    static Float_t crossSectionBRfilterEff, weight;
    static Char_t isPassed;
    tree->SetBranchAddress(
      "HGamEventInfoAuxDyn.crossSectionBRfilterEff",&crossSectionBRfilterEff);
    tree->SetBranchAddress(
      "HGamEventInfoAuxDyn.weight",&weight);
    tree->SetBranchAddress(
      "HGamEventInfoAuxDyn.isPassed",&isPassed);

    // LOOP over tree entries
    for (Long64_t nent=tree->GetEntries(), ent=0; ent<nent; ++ent) {
      tree->GetEntry(ent);

      for (size_t i=0, n=hist_types.size(); i<n; ++i) {
        auto& hist = hmap[var[i].second][i];
        test( var[i].first )
        test( var[i].second )
        test( hist.first->GetName() )
        hist.second->Fill(
          var[i].first,
          crossSectionBRfilterEff*weight*isPassed
        );
      }
      if (ent==10) break;
    }

    // Add histograms
    for (const auto& hs : hmap) {
      for (const auto& h : hs) {
        h.second->Scale(xsecscale/h.second->GetBinWidth(1));
        h.first->Add(h.second.get());
        h.second->Reset();
      }
    }

    delete file;
  }

  // Fit functions **************************************************
  workspace ws(wfname);

  // Draw histograms ************************************************
  TCanvas canv;

  canv.SaveAs((ofname+'[').c_str());

  for (const auto& hs : hmap) {
    int i=0;
    Color_t color;
    for (const auto& h : hs) {
      h.first->SetStats(false);
      h.first->SetLineWidth(2);
      h.first->SetLineColor(color = colors[(i++) % colors.size()]);
      h.first->Draw(i ? "" : "same");
    }
    canv.SaveAs(ofname.c_str());
  }

  canv.SaveAs((ofname+']').c_str());

  return 0;
}
