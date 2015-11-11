// Developed by Ivan Pogrebnyak, MSU

#include "TGraph_fcns.hh"

#include <string>
#include <sstream>
#include <stdexcept>

#include <iostream>

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

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

Double_t fint(const TGraph* gr, Int_t ixi, Int_t ixf) noexcept {
  const int n = gr->GetN();
  if (n<2 || n<ixf+1) {
    stringstream ss;
    ss << "Graph " << gr->GetName() << " has only " <<n<< " points";
    throw runtime_error(ss.str());
  }
  if ((ixf-ixi)<1) {
    stringstream ss;
    ss << (ixf-ixi+1) << " points selected in graph " << gr->GetName();
    throw runtime_error(ss.str());
  }
  Double_t integral = 0.;
  Double_t x1, x2, y1, y2;
  gr->GetPoint(ixi++, x1, y1);
  for (; ixi<=ixf; ++ixi) {
    gr->GetPoint(ixi, x2, y2);
    integral += (x2-x1)*(y2+y1);
    x1 = x2;
    y1 = y2;
  }
  return integral/2;
}

Double_t ltailx(const TGraph* gr, Double_t frac) noexcept {
  const int n = gr->GetN();
  const Double_t intfrac = fint(gr,0,n-1)*2*frac;
  Double_t integral = 0.;
  Double_t x1, x2, y1, y2;
  gr->GetPoint(0, x1, y1);
  for (int i=1; i<n; ++i) {
    gr->GetPoint(i, x2, y2);
    Double_t part = (x2-x1)*(y2+y1);
    if (integral+part > intfrac) {
      Double_t a = (y2-y1)/(x2-x1);
      Double_t b = y2 - a*x2;
      Double_t c = (integral - intfrac) - (b + y1)*x1;
      b += (y1 - a*x1);
      x2 = (sqrt(b*b-4.*a*c)-b)/(2.*a);
      break;
    } else {
      integral += part;
      x1 = x2;
      y1 = y2;
    }
  }
  return x2;
}

Double_t rtailx(const TGraph* gr, Double_t frac) noexcept {
  const int n = gr->GetN();
  const Double_t intfrac = fint(gr,0,n-1)*2*frac;
  Double_t integral = 0.;
  Double_t x1, x2, y1, y2;
  gr->GetPoint(n-1, x2, y2);
  for (int i=n-2; i>=0; --i) {
    gr->GetPoint(i, x1, y1);
    Double_t part = (x2-x1)*(y2+y1);
    if (integral+part > intfrac) {
      Double_t a = (y1-y2)/(x2-x1);
      Double_t b = y2 + a*x2;
      Double_t c = (integral - intfrac) + (y2 + b)*x2;
      b += (y2 + a*x2);
      x1 = (-sqrt(b*b-4.*a*c)+b)/(2.*a);
      break;
    } else {
      integral += part;
      x2 = x1;
      y2 = y1;
    }
  }
  return x1;
}

