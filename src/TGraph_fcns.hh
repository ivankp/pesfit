// Developed by Ivan Pogrebnyak, MSU

#ifndef TGraph_fcns_hh
#define TGraph_fcns_hh

#include <utility>
#include <TGraph.h>

inline Double_t firstx(const TGraph* gr) noexcept { return gr->GetX()[0]; }
inline Double_t  lastx(const TGraph* gr) noexcept { return gr->GetX()[gr->GetN()-1]; }

std::pair<Int_t,Double_t> max(const TGraph* gr) noexcept;
Double_t lfindx(const TGraph* gr, Double_t y) noexcept;
Double_t rfindx(const TGraph* gr, Double_t y) noexcept;

Double_t integrate(const TGraph* gr, Int_t ixi=0, Int_t ixf=-1) noexcept;
Double_t ltailx(const TGraph* gr, Double_t frac, Double_t totalint=0.) noexcept;
Double_t rtailx(const TGraph* gr, Double_t frac, Double_t totalint=0.) noexcept;

Double_t intervalx2(const TGraph* gr, Double_t frac, Double_t x1, Double_t totalint=0.) noexcept;

Double_t fwhm(const TGraph* gr) noexcept;

#endif
