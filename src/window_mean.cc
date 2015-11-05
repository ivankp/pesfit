#include "window_mean.hh"

template<typename T> inline T sq(T x) noexcept { return x*x; }

Double_t window_mean(const TH1* hist, Double_t a, Double_t b) noexcept {
  Int_t ai = hist->FindFixBin(a);
  Int_t bi = hist->FindFixBin(b);
  Double_t mean = 0., stdev = 0., sumw = 0.;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    const Double_t x = hist->GetBinCenter(i);
    sumw += w;
    mean += x*w;
  }
  mean /= sumw;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    const Double_t x = hist->GetBinCenter(i);
    stdev += sq(x-mean)*w;
  }
  stdev = sqrt(stdev/sumw);

  ai = hist->FindFixBin(mean-1.5*stdev);
  bi = hist->FindFixBin(mean+2.0*stdev);
  mean = 0.; sumw = 0.;
  for (int i=ai; i<=bi; ++i) {
    const Double_t w = hist->GetBinContent(i);
    sumw += w;
    mean += hist->GetBinCenter(i) * w;
  }
  return mean/sumw;
}
