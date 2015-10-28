// Developed by Ivan Pogrebnyak, MSU

#include "TGraph_fcns.hh"

#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;

pair<Int_t,Double_t> max(const TGraph* gr) noexcept {
  Int_t maxi = 0, n = gr->GetN();
  if (n==0) {
    stringstream ss;
    ss << "Graph " << gr->GetName() << " has 0 points";
    throw runtime_error(ss.str());
  }
  Double_t max, x, y;
  gr->GetPoint(0,x,max);
  for (int i=1; i<n; ++i) {
    gr->GetPoint(i,x,y);
    if (y>max) {
      max  = y;
      maxi = i;
    }
  }
  return {maxi,max};
}

Double_t lfindx(const TGraph* gr, Double_t y) noexcept {
  Double_t x2, y2;
  int i = 0, n = gr->GetN();
  for (; i<n; ++i) {
    gr->GetPoint(i,x2,y2);
    if (y2>=y) break;
  }
  if (i==0 || i==n) {
    stringstream ss;
    ss << "Graph " << gr->GetName()
       << " never attains value " << y;
    throw runtime_error(ss.str());
  }
  Double_t x1, y1;
  gr->GetPoint(i-1,x1,y1);
  const Double_t a = (y2-y1)/(x2-x1);
  const Double_t b = y2 - a*x2;
  return (y-b)/a;
}

Double_t rfindx(const TGraph* gr, Double_t y) noexcept {
  Double_t x1, y1;
  int n = gr->GetN(), i = n-1;
  for (; i>=0; --i) {
    gr->GetPoint(i,x1,y1);
    if (y1>=y) break;
  }
  if (i==-1 || i==n-1) {
    stringstream ss;
    ss << "Graph " << gr->GetName()
       << " never attains value " << y;
    throw runtime_error(ss.str());
  }
  Double_t x2, y2;
  gr->GetPoint(i+1,x2,y2);
  const Double_t a = (y2-y1)/(x2-x1);
  const Double_t b = y2 - a*x2;
  return (y-b)/a;
}

