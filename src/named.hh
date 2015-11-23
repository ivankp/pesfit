#ifndef named_hh
#define named_hh

#include <string>
#include <utility>
#include <iostream>
#include <cctype>

template<typename T>
struct named {
  std::string name;
  T x;

  named(): name() { }
  template<typename Name, typename... TT> named(Name&& name, TT&&... xx)
  : name(std::forward<Name>(name)), x(std::forward<TT>(xx)...) { }

  template<typename T1, typename T2> named(const std::pair<T1, T2>& p)
  : name(p.first), x(p.second) { }
  template<typename T1, typename T2> named(std::pair<T1, T2>&& p)
  : name(std::move(p.first)), x(std::move(p.second)) { }

  named(const named& p) = default;
  named(named&& p) = default;

  inline operator T&() { return x; }
  inline operator const T&() const { return x; }
  inline const char* cname() { return name.c_str(); }
};

template <typename T>
std::istream& operator>>(std::istream& in, named<T>& x) {
  // Skip over the leading whitespace
  while (!isgraph(in.peek())) in.get();
  std::getline(in, x.name, ':');
  in >> x.x;
  return in;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, named<T>& x) {
  out << x.name <<':'<< x.x;
  return out;
}

#endif
