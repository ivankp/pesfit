#ifndef pesfit_workspace_hh
#define pesfit_workspace_hh

#include <string>
#include <utility>
#include <memory>

#include <TH1.h>

#include <RooFitResult.h>

class TFile;
class TGraph;

class RooWorkspace;
class RooSimultaneous;
class RooCategory;
class RooRealVar;

using FitResult = std::unique_ptr<RooFitResult>;

class workspace {
  bool bg;
  TFile *file;
  RooWorkspace *ws;
  RooSimultaneous *sim_pdf;
  RooCategory *rcat;
  RooRealVar *myy;

public:
  workspace(const std::string& fname, bool bg=false);
  ~workspace();

  inline RooWorkspace* operator->() noexcept { return ws; }

  void setRange(const char* name, Double_t min, Double_t max);
  void fixVal(const char* name, Double_t val);

  std::pair<FitResult,TGraph*> fit(TH1* hist) const;
};

#endif
