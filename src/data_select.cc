#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "assert_ext.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

template<typename T>
inline void branches_impl(TTree* tree, const char* name, T* add) {
  tree->SetBranchAddress(name, add);
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus(name,1);
}

template<typename T, typename... TT>
inline void branches_impl(
  TTree* tree, const char* name, T* add, TT... bb
) {
  tree->SetBranchAddress(name, add);
  branches(tree,bb...);
  tree->SetBranchStatus(name,1);
}

template<typename... TT>
inline void branches(TTree* tree, TT... bb) {
  tree->SetBranchStatus("*",1);
  branches_impl(tree,bb...);
}

int main(int argc, char** argv)
{
  if (argc<3) {
    cout << "usage: " << argv[0]
         << " output.root input.root [input2.root ...]" << endl;
    return 1;
  }
  assert_ext(&argv[1],argc-1,{"root"});
  const int NF = argc - 2;

  TFile *fout = new TFile(argv[1],"recreate");
  TTree *tree = new TTree("mc","");

  Double_t myy, weight;
  Int_t npv;

  tree->Branch("m_yy",   &myy,    "m_yy/D"  );
  tree->Branch("npv",    &npv,    "npv/I"   );
  tree->Branch("weight", &weight, "weight/D");

  for (int i=2; i<argc; ++i) {
    TFile *f = new TFile(argv[i],"read");
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

        if (isPassed==1 && fmyy>105000. && fmyy<140000.) {
          weight_scale += eff/* *weight */;
        }
      }
      cout << "total_weight = " << weight_scale << endl;
      weight_scale = (1./weight_scale)/NF;

      branches(t,
        "HGamEventInfoAuxDyn.m_yy", &fmyy,
        "HGamEventInfoAuxDyn.numberOfPrimaryVertices", &npv,
        "HGamEventInfoAuxDyn.crossSectionBRfilterEff", &eff,
        // "HGamEventInfoAuxDyn.weight", &_weight,
        "HGamEventInfoAuxDyn.isPassed", &isPassed
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (isPassed==1 && fmyy>105000. && fmyy<140000.) {
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

        if (truth_mH==125 && myy>105000. && myy<140000.) {
          // weight_scale += weight;
          weight_scale += 1.;
        }
      }
      cout << "total_weight = " << weight_scale << endl;
      weight_scale = (1./weight_scale)/NF;

      branches(t,
        "m_yy", &myy,
        "NPV", &npv,
        "Higgs_truth_mass", &truth_mH
      );

      for (Long64_t nent=t->GetEntries(), ent=0; ent<nent; ++ent) {
        t->GetEntry(ent);

        if (truth_mH==125 && myy>105000. && myy<140000.) {
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
