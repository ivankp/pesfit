#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <thread>
#include <cstring>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooCurve.h>
#include <RooHist.h>

#include "root_safe_get.hh"
#include "assert_ext.hh"
#include "catstr.hh"
#include "TGraph_fcns.hh"
#include "golden_min.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

struct workspace {
  TFile *file;
  RooWorkspace *ws;
  RooSimultaneous *sim_pdf;
  RooCategory *cat;
  RooRealVar *myy;

public:
  workspace(const char* fname)
  : file(new TFile(fname,"read")),
    ws(get<RooWorkspace>(file,"mcfit")),
    sim_pdf(static_cast<RooSimultaneous*>(
      ws->obj("mc_sim_pdf_bin0"))),
    cat(static_cast<RooCategory*>(ws->obj("mc_sample"))),
    myy(static_cast<RooRealVar*>(ws->var("m_yy")))
  { }
  ~workspace()  {
    delete myy;
    delete cat;
    delete sim_pdf;
    delete file;
  }

  inline RooWorkspace* operator->() noexcept { return ws; }
};

class mc_tree_ptr {
  TFile *file;
  TTree *tree;
public:
  mc_tree_ptr(const char* fname)
  : file(new TFile(fname,"read")),
    tree(get<TTree>(file,"mc"))
  {
    cout << file->GetName() << endl;
  }
  ~mc_tree_ptr() { delete file; }
  inline TTree* operator->() noexcept { return tree; }
  inline operator TTree*() const noexcept { return tree; }
};

int main(int argc, char** argv)
{
  if (argc!=4 && argc!=5 && argc!=6) {
    cout << "usage: " << argv[0]
         << " workspace.root input.root output.(root|pdf)"
            " name [-w] [-l]" << endl;
    return 1;
  }
  assert_ext(argv+1,2,{"root"});
  assert_ext(argv[3],{"root","pdf"});

  bool w = false, logy = false;
  string name;
  for (int i=4; i<argc; ++i) {
    if (argv[i][0]!='-') name = argv[i];
    else if (!strcmp(argv[i],"-w")) w = true;
    else if (!strcmp(argv[i],"-l")) logy = true;
  }

  workspace ws(argv[1]);
  mc_tree_ptr tree(argv[2]);
  const string ofname(argv[3]);

  ws.myy->setBins(70);

  RooRealVar wvar("weight","weight",0,2);

  // Build merged RooDataSet
  RooDataSet set("m_yy_set", "m_yy_set",
    w ? RooArgSet(*ws.myy,wvar) : RooArgSet(*ws.myy),
    RooFit::Index(*ws.cat),
    RooFit::Import(map<string,RooDataSet*> {
      {"mc_125", new RooDataSet(
        name.c_str(), name.c_str(),
        tree,
        w ? RooArgSet(*ws.myy,wvar) : RooArgSet(*ws.myy),
        0,
        w ? "weight" : 0
      )}
    }),
    w ? RooFit::WeightVar(wvar) : RooCmdArg()
  );

  cout << endl << endl;
  ws.cat->Print("v");
  cout << endl << endl;
  set.Print("v");
  cout << endl << endl;

  // -------------------------------------------------
  // Perform the Fit
  RooFitResult* fit = ws.sim_pdf->fitTo(set,
    RooFit::Extended(false),
    RooFit::SumW2Error(w),
    RooFit::Save(true),
    RooFit::NumCPU(std::thread::hardware_concurrency()),
    RooFit::Strategy(2),
    RooFit::Minimizer("Minuit2"),
    RooFit::Offset(true)
  );

  fit->Print("v");

  RooPlot *frame = ws.myy->frame(RooFit::Title(name.size() ? name.c_str() : argv[2]));
  // Draw Monte Carlo histogram
  set.plotOn(frame,
    RooFit::MarkerColor(50),
    RooFit::MarkerSize(0.5)
  );
  // Draw fitted function
  ws.sim_pdf->plotOn(frame,
    RooFit::LineColor(85),
    RooFit::Slice(*ws.cat, "mc_125"),
    RooFit::ProjWData(RooArgSet(*ws.cat), set),
    RooFit::Precision(1e-5)
  );
  
  TFile *file = (ofname.substr(ofname.rfind('.')+1)=="root")
    ? new TFile(ofname.c_str(),"recreate") : nullptr;

  TCanvas *canv = new TCanvas();
  if (!file) {
    canv->SetMargin(0.1,0.04,0.1,0.1);
    if (logy) canv->SetLogy();
  }

  frame->Draw();

  if (file) {
    file = new TFile(ofname.c_str(),"recreate");
    TTree *tree = new TTree("stats","CB+GA fit parameters");
    
    auto *hist = new TGraphAsymmErrors(*frame->getHist("h_m_yy_set"));
    hist->SetTitle("MC pseudo-data");
    hist->Write("hist");
    auto *fcn = new TGraph(*frame->getCurve());
    fcn->SetTitle("Fitted function");
    fcn->Write("fcn");
    
    vector<pair<string,pair<Double_t,Double_t>>> pars {
      {"crys_alpha_bin0",{}},
      {"crys_norm_bin0",{}},
      {"fcb_bin0",{}},
      {"gaus_kappa_bin0",{}},
      {"gaus_mean_offset_bin0",{}},
      {"mean_offset_bin0",{}},
      {"sigma_offset_bin0",{}},
      {"sigma_68",{}}
    };
    auto par = pars.begin();
    for (auto end=pars.end()-1; par<end; ++par) {
      auto *var = static_cast<RooRealVar*>(
        fit->floatParsFinal().find(par->first.c_str()));
      par->second = {var->getVal(),var->getError()};
    }
    
    const Double_t integral = integrate(fcn);
    const double frac = 0.68;
    golden_min gm;
    auto x1 = gm( [fcn,frac,integral](double x1) {
      return intervalx2(fcn, frac, x1, integral) - x1;
    }, firstx(fcn), rtailx(fcn,frac,integral) ).first;
    auto x2 = intervalx2(fcn, frac, x1, integral);
    par->second = {x1,x2};
    
    for (const auto& par : pars)
      tree->Branch(par.first.c_str(), (void*)&par.second,
                   (par.first+"[2]/D").c_str());

    tree->Fill();
    file->Write(0,TObject::kOverwrite);
    delete file;
  } else {
    canv->SaveAs(ofname.c_str());
  }
  delete canv;

  delete fit;
  return 0;
}
