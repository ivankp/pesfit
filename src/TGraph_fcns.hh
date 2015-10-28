// Developed by Ivan Pogrebnyak, MSU

#ifndef TGraph_fcns_hh
#define TGraph_fcns_hh

#include <utility>
#include <TGraph.h>

std::pair<Int_t,Double_t> max(const TGraph* gr) noexcept;
Double_t lfindx(const TGraph* gr, Double_t y) noexcept;
Double_t rfindx(const TGraph* gr, Double_t y) noexcept;

#endif
