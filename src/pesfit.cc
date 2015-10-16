#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>

#include "catstr.hh"

using namespace std;
namespace po = boost::program_options;

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
string ofname, cfname;
vector<string> ifname;
vector<Color_t> colors;
Int_t nbins;
// --------------------------

template<typename T>
inline void SetColor(T* obj, int i) noexcept {
  const Color_t color = colors[i%colors.size()];
  obj->SetLineColor(color);
  obj->SetMarkerColor(color);
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
  
  TH1 *hist, *nom, *scale_down, *scale_up, *res_down, *res_up;

  for (const string& f : ifname) {
    static bool first = true;

    TFile *file = new TFile(f.c_str(),"read");
    if (file->IsZombie()) return 1;
    cout << f << endl;
    TTree *tree = get<TTree>(file,"CollectionTree");
    
    const double xsecscale =
      1./get<TH1>(file,"CutFlow_PowhegPy8_ggH125_small_weighted")->GetBinContent(3);
    // regex to find the right hist
    // convert to differential cross section
    // fitt gaussian and crystal ball
    
    // Is HGamEventInfoAuxDyn.crossSectionBRfilterEff name?
    
    // --------------------------------------------------------------
    
    tree->Draw(
      cat("HGamEventInfoAuxDyn.m_yy/1000>>hist(",nbins,",105,160)").c_str(),
      "HGamEventInfoAuxDyn.weight*(HGamEventInfoAuxDyn.isPassed==1)"
    );
    hist = get<TH1>(gDirectory,"hist");
    hist->Scale(xsecscale);
    if (first)
      (nom = (TH1*)hist->Clone("nom"))->SetDirectory(0);
    else nom->Add(hist);
    
    // --------------------------------------------------------------
    
    tree->Draw(
      cat("HGamEventInfo_EG_SCALE_ALL__1downAuxDyn.m_yy/1000>>hist(",nbins,",105,160)").c_str(),
      "HGamEventInfoAuxDyn.weight*(HGamEventInfoAuxDyn.isPassed==1)"
    );
    hist = get<TH1>(gDirectory,"hist");
    hist->Scale(xsecscale);
    if (first)
      (scale_down = (TH1*)hist->Clone("scale_down"))->SetDirectory(0);
    else scale_down->Add(hist);
    
    tree->Draw(
      cat("HGamEventInfo_EG_SCALE_ALL__1upAuxDyn.m_yy/1000>>hist(",nbins,",105,160)").c_str(),
      "HGamEventInfoAuxDyn.weight*(HGamEventInfoAuxDyn.isPassed==1)"
    );
    hist = get<TH1>(gDirectory,"hist");
    hist->Scale(xsecscale);
    if (first)
      (scale_up = (TH1*)hist->Clone("scale_up"))->SetDirectory(0);
    else scale_up->Add(hist);
    
    // --------------------------------------------------------------
    
    tree->Draw(
      cat("HGamEventInfo_EG_RESOLUTION_ALL__1downAuxDyn.m_yy/1000>>hist(",nbins,",105,160)").c_str(),
      "HGamEventInfoAuxDyn.weight*(HGamEventInfoAuxDyn.isPassed==1)"
    );
    hist = get<TH1>(gDirectory,"hist");
    hist->Scale(xsecscale);
    if (first)
      (res_down = (TH1*)hist->Clone("res_down"))->SetDirectory(0);
    else res_down->Add(hist);
    
    tree->Draw(
      cat("HGamEventInfo_EG_RESOLUTION_ALL__1upAuxDyn.m_yy/1000>>hist(",nbins,",105,160)").c_str(),
      "HGamEventInfoAuxDyn.weight*(HGamEventInfoAuxDyn.isPassed==1)"
    );
    hist = get<TH1>(gDirectory,"hist");
    hist->Scale(xsecscale);
    if (first)
      (res_up = (TH1*)hist->Clone("res_up"))->SetDirectory(0);
    else res_up->Add(hist);
    
    // --------------------------------------------------------------
    
    delete file;
    first = false;
  }

  // ---------------------------------------

  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  canv.SaveAs((ofname+'[').c_str());

  // ---------------------------------------
    
  nom->SetStats(false);
  nom->SetLineWidth(2);
  SetColor(nom,0);
  nom->Draw();
  canv.SaveAs(ofname.c_str());

  // ---------------------------------------
    
  scale_down->SetStats(false);
  scale_down->SetLineWidth(2);
  SetColor(scale_down,0);
  scale_down->Draw();
    
  scale_up->SetStats(false);
  scale_up->SetLineWidth(2);
  SetColor(scale_up,1);
  scale_up->Draw("same");
  canv.SaveAs(ofname.c_str());

  // ---------------------------------------
    
  res_down->SetStats(false);
  res_down->SetLineWidth(2);
  SetColor(res_down,0);
  res_down->Draw();
    
  res_up->SetStats(false);
  res_up->SetLineWidth(2);
  SetColor(res_up,1);
  res_up->Draw("same");
  canv.SaveAs(ofname.c_str());

  // ---------------------------------------

  canv.SaveAs((ofname+']').c_str());

  return 0;
}
