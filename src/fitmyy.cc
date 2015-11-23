#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooCurve.h>

#include "root_safe_get.hh"
#include "assert_ext.hh"
#include "named.hh"
#include "catstr.hh"

using namespace std;

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
  if (argc<4) {
    cout << "usage: " << argv[0]
         << " output.(root|pdf) workspace.root"
            " name:input.root [name:input2.root ...]" << endl;
    return 1;
  }
  assert_ext(argv[1],{"root","pdf"});
  assert_ext(&argv[2],argc-2,{"root"});

  workspace ws(argv[2]);

  vector<named<mc_tree_ptr>> trees;
  trees.reserve(argc-3);
  for (int i=3; i<argc; ++i) {
    named<string> fname;
    stringstream(argv[i]) >> fname;
    trees.emplace_back(move(fname.name),fname.x.c_str());
  }

  // -------------------------------------------------
  // Fetch samples and
  // Merge data sets into a single combined data set
  RooRealVar wvar("weight","weight",0,2);
  map<string,RooDataSet*> sets;
  for (auto& tree : trees) {
    sets.emplace( tree.name, new RooDataSet(
      tree.cname(), tree.cname(),
      tree.x, RooArgSet(*ws.myy,wvar), 0, "weight"
    ) );
  }
  // Build merged RooDataSet
  RooDataSet set("m_yy_set", "m_yy_set",
    RooArgSet(*ws.myy,wvar),
    RooFit::Index(*ws.cat),
    RooFit::Import(sets),
    RooFit::WeightVar(wvar));

  // -------------------------------------------------
  // Perform the Fit
  RooFitResult* fit = ws.sim_pdf->fitTo(set,
    RooFit::Extended(false),
    RooFit::SumW2Error(true),
    RooFit::Save(true),
    RooFit::NumCPU(4),
    RooFit::Strategy(2)
  );

  fit->Print("v");

  TCanvas canv;
  canv.SetMargin(0.1,0.04,0.1,0.1);
  canv.SetLogy();
  // canv.SaveAs(cat(argv[1],'[').c_str());

  for (auto& tree : trees) {
    static int i=0;
    RooPlot *frame = ws.myy->frame(RooFit::Title(" "/*tree.cname()*/));
    // Draw Monte Carlo histogram
    set.plotOn(frame,
      RooFit::MarkerColor(40+10*i),
      RooFit::MarkerSize(0.5),
      RooFit::Cut(cat("mc_sample == mc_sample::",tree.name).c_str())
    );
    // Draw fitted function
    // !!!!!!!!!!!!!!!!!!!!!!!!!! I get a segfault here!!!
    // ws.sim_pdf->plotOn(frame,
    //   RooFit::LineColor(85),
    //   RooFit::Slice(*ws.cat, tree.cname()),
    //   RooFit::ProjWData(RooArgSet(*ws.cat), set)
    //   // RooFit::Precision(1e-5)
    // );
    // !!!!!!!!!!!!!!!!!!!!!!!!!!

    // auto *curve = frame->getCurve();
    frame->Draw(i++ ? "same" : "");
  }
  canv.SaveAs(argv[1]);
  // canv.SaveAs(cat(argv[1],']').c_str());

  delete fit;
  return 0;
}
