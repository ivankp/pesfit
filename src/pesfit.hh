#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <initializer_list>
#include <memory>
#include <stdexcept>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooCurve.h>

#include "regex.hh"
#include "catstr.hh"
#include "senum.hh"
#include "seqmap.hh"
#include "structmap.hh"
#include "val_err.hh"
#include "TGraph_fcns.hh"
#include "root_safe_get.hh"
#include "workspace.hh"
#include "window_mean.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::pair;
using std::vector;
using std::map;
using std::unordered_set;
using std::initializer_list;
using std::unique_ptr;
using std::runtime_error;
using std::stringstream;
using std::scientific;
using std::setprecision;
using std::fixed;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;
