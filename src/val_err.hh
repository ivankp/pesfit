// Developed by Ivan Pogrebnyak, MSU

#ifndef val_err_hh
#define val_err_hh

#include <ostream>
#include <iomanip>
#include <utility>

// number of precision digits
int err_prec(double err, int n=2) noexcept {
  const double log_err = std::log10(err);
  if (log_err < -8) return 0;
  else if (log_err < 0) return n - log_err;
  else if (log_err < n-1) return n-1;
  else return 0;
}

template<typename V> struct val_err {
private:
  template<typename T> inline T sq(T x) noexcept { return x*x; }

public:
  typedef V type;
  type val, err;

  val_err(): val(0), err(0) { }
  template<typename T1>
  val_err(const T1& val): val(val), err(0) { }
  template<typename T1, typename T2>
  val_err(const T1& val, const T2& err): val(val), err(err) { }
  template<typename T1, typename T2>
  val_err(const std::pair<T1,T2>& x): val(x.first), err(x.second) { }

  // Cast to value type
  inline operator type() const noexcept { return val; }

  // Addition
  template <typename T>
  inline val_err<type> operator+(const T& x) const noexcept {
    return {val+x, err};
  }
  template <typename T>
  inline val_err<type> operator+(const val_err<T>& x) const noexcept {
    return {val + x.val, sqrt( sq(err) + sq(x.err) )};
  }

  // Difference
  template <typename T>
  inline val_err<type> operator-(const T& x) const noexcept {
    return {val-x, err};
  }
  template <typename T>
  inline val_err<type> operator-(const val_err<T>& x) const noexcept {
    return {val - x.val, sqrt( sq(err) + sq(x.err) )};
  }

  // Product
  template <typename T>
  inline val_err<type> operator*(const T& x) const noexcept {
    return {val*x, err*x};
  }
  template <typename T>
  inline val_err<type> operator*(const val_err<T>& x) const noexcept {
    const auto a = val * x.val;
    return {a, a*sqrt( sq(err/val) + sq(x.err/x.val) )};
  }

  // Quotient
  template <typename T>
  inline val_err<type> operator/(const T& x) const noexcept {
    return {val/x, err/x};
  }
  template <typename T>
  inline val_err<type> operator/(const val_err<T>& x) const noexcept {
    const auto a = val / x.val;
    return {a, a*sqrt( sq(err/val) + sq(x.err/x.val) )};
  }

  // Print
  std::ostream& print(std::ostream& out, const char* pm=" ± ") const {
    const int prec = err_prec(err);
    const auto flags = out.flags();
    out << std::fixed << std::setprecision(prec) << val << pm
        << std::setprecision(prec) << err;
    out.flags(flags);
    return out;
  }
};

// output stream operator
template<typename T>
std::ostream& operator<<(std::ostream& out, const val_err<T>& x) {
  return x.print(out);
}

// output stream operator
template<typename T>
std::ostream& operator<<(std::ostream& out, const std::pair<val_err<T>,const char*>& x) {
  return x.first.print(out,x.second);
}

#endif
