// Developed by Ivan Pogrebnyak, MSU

#include "TGraph_fcns.hh"

#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include <iostream>

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using namespace std;

pair<Int_t,Double_t> max(const TGraph* gr) noexcept {
  Int_t maxi = 0, n = gr->GetN();
  if (n==0) {
    stringstream ss;
    ss << "Graph " << gr->GetName() << " has 0 points";
    throw out_of_range(ss.str());
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
    throw out_of_range(ss.str());
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
    throw out_of_range(ss.str());
  }
  Double_t x2, y2;
  gr->GetPoint(i+1,x2,y2);
  const Double_t a = (y2-y1)/(x2-x1);
  const Double_t b = y2 - a*x2;
  return (y-b)/a;
}

Double_t integrate(const TGraph* gr, Int_t ixi/*=0*/, Int_t ixf/*=-1*/) noexcept {
  const auto n = gr->GetN();
  if (ixf==-1) ixf = n-1;
  if (n<2 || n<ixf+1) {
    stringstream ss;
    ss << "Graph " << gr->GetName() << " has only " <<n<< " points";
    throw out_of_range(ss.str());
  }
  if ((ixf-ixi)<1) {
    stringstream ss;
    ss << (ixf-ixi+1) << " points selected in graph " << gr->GetName();
    throw out_of_range(ss.str());
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

Double_t ltailx(const TGraph* gr, Double_t frac, Double_t totalint/*=0.*/) noexcept {
  const auto n = gr->GetN();
  if (totalint==0.) totalint = integrate(gr);
  totalint *= frac*2;
  Double_t integral = 0.;
  Double_t x1, x2, y1, y2;
  gr->GetPoint(0, x1, y1);
  for (auto i=1; i<n; ++i) {
    gr->GetPoint(i, x2, y2);
    Double_t part = (x2-x1)*(y2+y1);
    if (integral+part > totalint) {
      Double_t a = (y2-y1)/(x2-x1);
      Double_t b = y2 - a*x2;
      Double_t c = (integral - totalint) - (b + y1)*x1;
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

Double_t rtailx(const TGraph* gr, Double_t frac, Double_t totalint/*=0.*/) noexcept {
  const auto n = gr->GetN();
  if (totalint==0.) totalint = integrate(gr);
  totalint *= frac*2;
  Double_t integral = 0.;
  Double_t x1, x2, y1, y2;
  gr->GetPoint(n-1, x2, y2);
  for (auto i=n-2; i>=0; --i) {
    gr->GetPoint(i, x1, y1);
    Double_t part = (x2-x1)*(y2+y1);
    if (integral+part > totalint) {
      Double_t a = (y1-y2)/(x2-x1);
      Double_t b = y2 + a*x2;
      Double_t c = (integral - totalint) + (y2 + b)*x2;
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

Double_t intervalx2(const TGraph* gr, Double_t frac, Double_t x1, Double_t totalint/*=0.*/) noexcept {
  const auto n = gr->GetN();
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  if (x1>*(x+n-1)) {
    stringstream ss;
    ss << x1 << " > " << *(x+n-1);
    throw out_of_range(ss.str());
  }
  auto k = upper_bound(x, x+n, x1) - x;
  if (k==0) {
    stringstream ss;
    ss << x1 << " < " << *x;
    throw out_of_range(ss.str());
  }

  if (totalint==0.) totalint = integrate(gr);
  totalint *= frac*2;

  Double_t a = (y[k]-y[k-1])/(x[k]-x[k-1]);
  Double_t b = y[k] - a*x[k];
  Double_t integral = (x[k]-x1)*(y[k]+(a*x1+b));
  Double_t x2 = x1;

  bool enough = false;
  for (++k; k<n; ++k) {
    Double_t part = (x[k]-x[k-1])*(y[k]+y[k-1]);
    if (integral+part > totalint) {
      a = (y[k]-y[k-1])/(x[k]-x[k-1]);
      b = y[k] - a*x[k];
      Double_t c = (integral - totalint) - (b + y[k-1])*x[k-1];
      b += (y[k-1] - a*x[k-1]);
      x2 = (sqrt(b*b-4.*a*c)-b)/(2.*a);
      enough = true;
      break;
    } else integral += part;
  }

  if (!enough) {
    stringstream ss;
    ss << "Integral from x1=" << x1 << " is less then "
       << frac << " of the total";
    throw out_of_range(ss.str());
  }

  return x2;
}

Double_t fwhm(const TGraph* gr) noexcept {
  const Double_t half_max = max(gr).second/2;
  return rfindx(gr,half_max) - lfindx(gr,half_max);
}
