// Developed by Ivan Pogrebnyak, MSU
// Adopted from Numerical Recipes, The art of scientific computing,
// 3rd ed., p.459

#ifndef golden_min_hh
#define golden_min_hh

#include <utility>
#include <cmath>

#include <iostream>
#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

class golden_min {
  const double tol;

public:
  golden_min(double tol=3e-8): tol(tol) { }
  
  template<typename Fcn>
  std::pair<double,double> operator()(Fcn fcn, double xmin, double xmax) {
    static constexpr double R=0.61803399, C=1.-R;
    double ax=xmin, cx=xmax, bx=(ax+cx)/2., x1, x2, x0=ax, x3=cx;
    if (cx < ax) std::swap(ax,cx);
    if (abs(cx-bx) > abs(bx-ax)) {
      x1=bx;
      x2=bx+C*(cx-bx);
    } else {
      x2=bx;
      x1=bx-C*(bx-ax);
    }
    double f1=fcn(x1);
    double f2=fcn(x2);
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
      if (f2 < f1) {
        shift(R*x2+C*x3,x2,x1,x0);
        shift(fcn(x2),f2,f1);
      } else {
        shift(R*x1+C*x0,x1,x2,x3);
        shift(fcn(x1),f1,f2);
      }
    }
    if (f1 < f2) return {x1,f1};
    else return {x2,f2};
  }
  
  template<typename A, typename B>
  static inline void shift(const B& b, A& a) noexcept {
    a = b;
  }
  template<typename C, typename B, typename... AA>
  static inline void shift(const C& c, B& b, AA&... aa) noexcept {
    shift(b,aa...);
    b = c;
  }

};

#endif
