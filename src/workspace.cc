#include "workspace.hh"

#include <map>

#include <TFile.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooCurve.h>

#include "root_safe_get.hh"

workspace::workspace(const std::string& fname, bool bg)
: bg(bg),
  file(new TFile(fname.c_str(),"read")),
  ws(get<RooWorkspace>(file,bg ? "datafit" : "mcfit")),
  sim_pdf(static_cast<RooSimultaneous*>(
    ws->obj(bg ? "sig_bkg_sim_pdf" : "mc_sim_pdf_bin0"))),
  rcat(static_cast<RooCategory*>(ws->obj(bg ? "sample" : "mc_sample"))),
  myy(static_cast<RooRealVar*>(ws->var("m_yy")))
{ }

workspace::~workspace() {
  delete myy;
  delete rcat;
  delete sim_pdf;
  delete file;
}

void workspace::setRange(const char* name, Double_t min, Double_t max) {
  ws->var(name)->setRange(min,max);
}

void workspace::fixVal(const char* name, Double_t val) {
  auto *var = ws->var(name);
  var->setRange(val,val);
  var->setVal(val);
}

auto workspace::fit(TH1* hist) const -> std::pair<FitResult,TGraph*> {
  // Produce a RooDataHist object from the TH1
  RooDataHist rdh("dh","dh",RooArgSet(*myy),hist);

  // This map might look a bit useless,
  // but the PDF is designed to fit several datasets simultaneously.
  // This is done by assigning the PDF categories,
  // here we only have one histogram
  std::map<std::string,RooDataHist*> rdhmap {
    {bg ? "data_bin0" : "mc_125", &rdh}
  };

  // With the map we can build one combined dataset
  // that has the proper link of the category information.
  RooDataHist crdh("c_dh","c_dh",
    RooArgSet(*myy), RooFit::Index(*rcat), RooFit::Import(rdhmap));

  // Now we are ready to fit! We have a PDF and a RooDataHist
  RooFitResult *res = sim_pdf->fitTo(crdh,
    RooFit::Extended(bg),
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
    RooFit::Cut(bg ? "sample == sample::data_bin0"
                   : "mc_sample == mc_sample::mc_125")
  );
  sim_pdf->plotOn(frame,
    RooFit::LineColor(85),
    RooFit::Slice(*rcat, bg ? "data_bin0" : "mc_125"),
    RooFit::ProjWData(RooArgSet(*rcat), crdh)
  );
  auto *curve = frame->getCurve();
  // frame->Draw("same");
  // curve->Draw("same");
  res->Print("v");

  return {FitResult(res),curve};
}
