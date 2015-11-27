#include <iostream>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>

#include "assert_ext.hh"
#include "branches.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

namespace std {
  template <typename T1, typename T2>
  istream& operator>>(istream& in, pair<T1,T2>& p) {
    string s;
    in >> s;
    size_t sep = s.find(':');
    if (sep==string::npos) throw invalid_argument(
      '\"'+s+"\": pair values must be delimited by \':\'");
    stringstream (s.substr(0,sep)) >> p.first;
    stringstream (s.substr(sep+1)) >> p.second;
    return in;
  }
}

int main(int argc, char** argv)
{
  vector<string> ifname;
  string ofname, cfname;
  pair<double,double> xrange;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->multitoken()->required(),
       "*input root file names")
      ("output,o", po::value(&ofname)->required(),
       "*output pdf or root file name")
      ("config,c", po::value(&cfname),
       "configuration file name")

      ("xrange,x", po::value(&xrange)->default_value({105,140},"105:140"),
       "m_yy data range")
    ;

    po::positional_options_description pos;
    pos.add("output",1);
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

    assert_ext(ifname,{"root"});
    assert_ext(ofname,{"root"});

  } catch (std::exception& e) {
    cerr << "\033[31mArgs: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  xrange.first  *= 1e3;
  xrange.second *= 1e3;

  TFile *fout = new TFile(ofname.c_str(),"recreate");
  TTree *tree = new TTree("mc","");

  Double_t myy, weight;
  Int_t npv;

  tree->Branch("m_yy",   &myy,    "m_yy/D"  );
  tree->Branch("npv",    &npv,    "npv/I"   );
  tree->Branch("weight", &weight, "weight/D");

  for (const auto& fname : ifname) {
    TFile *f = new TFile(fname.c_str(),"read");
    if (f->IsZombie()) return 1;
    cout << f->GetName() << endl;
    TTree *t = nullptr;
    Double_t weight_scale = 0.;

    if ((t = dynamic_cast<TTree*>(f->Get("CollectionTree")))) {

      Float_t fmyy;
      Float_t eff;
      // Float_t _weight;
      Char_t isPassed;

      branches(t,
        "HGamEventInfoAuxDyn.m_yy", &fmyy,
        "HGamEventInfoAuxDyn.crossSectionBRfilterEff", &eff,
        // "HGamEventInfoAuxDyn.weight", &weight,
        "HGamEventInfoAuxDyn.isPassed", &isPassed
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (isPassed==1 && fmyy>xrange.first && fmyy<xrange.second) {
          weight_scale += eff/* *weight */;
        }
      }
      cout << "total_weight = " << weight_scale << endl;
      weight_scale = (1./weight_scale)/ifname.size();

      branches(t,
        "HGamEventInfoAuxDyn.m_yy", &fmyy,
        "HGamEventInfoAuxDyn.numberOfPrimaryVertices", &npv,
        "HGamEventInfoAuxDyn.crossSectionBRfilterEff", &eff,
        // "HGamEventInfoAuxDyn.weight", &_weight,
        "HGamEventInfoAuxDyn.isPassed", &isPassed
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (isPassed==1 && fmyy>xrange.first && fmyy<xrange.second) {
          myy = fmyy/1e3;
          weight = eff * weight_scale;
          tree->Fill();
        }
      }

    } else if ((t = dynamic_cast<TTree*>(f->Get("Hgg_tree")))) {

      // Double_t _weight;
      Int_t truth_mH;

      branches(t,
        "m_yy", &myy,
        // "weight", &_weight,
        "Higgs_truth_mass", &truth_mH
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (truth_mH==125 && myy>xrange.first && myy<xrange.second) {
          // weight_scale += weight;
          weight_scale += 1.;
        }
      }
      cout << "total_weight = " << weight_scale << endl;
      weight_scale = (1./weight_scale)/ifname.size();

      branches(t,
        "m_yy", &myy,
        "NPV", &npv,
        "Higgs_truth_mass", &truth_mH
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (truth_mH==125 && myy>xrange.first && myy<xrange.second) {
          myy /= 1e3;
          weight = weight_scale;
          tree->Fill();
        }
      }

    } else {
      cerr << "Neither CollectionTree nor Hgg_tree found in file "
           << f->GetName() << endl;
      return 1;
    }

    delete f;
  }

  fout->Write(0,TObject::kOverwrite);
  delete fout;

  return 0;
}
