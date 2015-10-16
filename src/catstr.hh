// Developed by Ivan Pogrebnyak, MSU

#ifndef catstr_hh
#define catstr_hh

#include <string>
#include <sstream>
#include <utility>

template<typename T>
inline void cat_impl(std::stringstream& ss, T&& t) {
  ss << t;
}

template<typename T, typename... TT>
inline void cat_impl(std::stringstream& ss, T&& t, TT&&... tt) {
  ss << t;
  cat_impl(ss,std::forward<TT>(tt)...);
}

template<typename... TT>
inline std::string cat(TT&&... tt) {
  std::stringstream ss;
  cat_impl(ss,std::forward<TT>(tt)...);
  return ss.str();
}

#endif
